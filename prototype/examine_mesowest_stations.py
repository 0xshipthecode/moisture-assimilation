#
#
#  This script attempted to find a nice spatial correlation structure in the observation
#  of fuel moisture (fm vs. equilibrium residualds) using moisture equilibria computed
#  from relative humidity and temperature measured directly at the site of the station.
#  However, comparisons show that the fuel moisture from the sensor is very different
#  from the computed equilibrium.  This indicates the residuals are unusable and indeed
#  there is a very weak spatial structure.
#
#
from time_series_utilities import match_sample_times, build_observation_data
from spatial_model_utilities import render_spatial_field, great_circle_distance
                                    
from wrf_model_data import WRFModelData
from mean_field_model import MeanFieldModel
from observation_stations import MesoWestStation, Observation
from statistics import compute_ols_estimator
from diagnostics import init_diagnostics, diagnostics

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
import pytz
import os
import string


station_data_dir = "../real_data/colorado_stations/"
wrf_data_file = "wrfout_sel_1km.nc"

def compute_model_equi_fm(H_Percent, T_Kelvin):
    """
    Compute the equilibrium temperature given by the
    moisture model from the relative humidity measured at the station
    and the temperature measured at the station.
    """
    H = H_Percent
    T = T_Kelvin

    Ed = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    Ew = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
                
    Ed *= 0.01
    Ew *= 0.01

    return 0.5 * (Ed + Ew)


if __name__ == '__main__':
    
    init_diagnostics('results/examine_mesowest_station_data.log')

    # load stations and match then to grid points
    # load station data from files
    with open(os.path.join(station_data_dir, 'station_infos'), 'r') as f:
        si_list = f.read().split('\n')

    si_list = filter(lambda x: len(x) > 0, map(string.strip, si_list))

    # load wrf model data
    wrf_data = WRFModelData(os.path.join(station_data_dir, wrf_data_file),
                            fields = ['T2', 'Q2', 'PSFC', 'EFMS', 'FMC_EQUI'],
                            tz_name = 'GMT')
    lon, lat, tm = wrf_data.get_lons(), wrf_data.get_lats(), wrf_data.get_times()
    efms, wrf_equi = wrf_data['EFMS'], wrf_data['FMC_EQUI']
    T2 = wrf_data['T2']

    print('Loaded %d times from WRF output (%s - %s).' % (len(tm), str(tm[0]), str(tm[-1])))

    # for each station id, load the station
    stations = []
    for sinfo in si_list:
        code = sinfo.split(',')[0]
        mws = MesoWestStation(sinfo, wrf_data)
        for suffix in [ '_1', '_2', '_3', '_4', '_5', '_6', '_7' ]:
            mws.load_station_data(os.path.join(station_data_dir, '%s%s.xls' % (code, suffix)))
        stations.append(mws)

    print('Loaded %d stations.' % len(stations))

    # check stations for nans
    stations = filter(MesoWestStation.data_ok, stations)
    print('Have %d stations with complete data.' % len(stations))

    # setup measurement variances
    for st in stations:
        st.set_measurement_variance('fm10', 0.2)
        st.set_measurement_variance('air_temp', 0.5)
        st.set_measurement_variance('rh', 10)

    # compute the mean of the equilibria
    Ed, Ew = wrf_data.get_moisture_equilibria()
    E = 0.5 * (Ed + Ew)

    mfm = MeanFieldModel()

    # construct basemap for rendering
    llcrnrlon = min(map(lambda x: x.lon, stations))
    llcrnrlat = min(map(lambda x: x.lat, stations))
    urcrnrlon = max(map(lambda x: x.lon, stations))
    urcrnrlat = max(map(lambda x: x.lat, stations))

    m = Basemap(llcrnrlon=llcrnrlon,
                llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon,
                urcrnrlat=urcrnrlat,
                projection = 'mill')

    # show the equilibrium field and render position of stations on top
    render_spatial_field(m, lon, lat, E[0,:,:], 'Equilibrium moisture')
    plt.figure()
    m.drawparallels(np.arange(llcrnrlat,urcrnrlat,1.0))
    m.drawmeridians(np.arange(llcrnrlon,urcrnrlon,1.0))
    for s in stations:
        slon, slat = s.get_position()
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')
        plt.text(x, y, s.get_name(), color = 'white')

    x, y = m(lon.ravel(), lat.ravel())
    plt.plot(x, y, 'k+', markersize = 4)
    plt.savefig('results/station_localization.png')

    # find common observation times for the station and for the WRF model
    tms = stations[0].get_obs_times()
    mtm, wrf_tndx, _ = match_sample_times(tm, tms)

    # part C - fit the mean field model
    obs_data = build_observation_data(stations, 'fm10', wrf_data, mtm)
    Nt = len(obs_data)
    gammas = []
    residuals = {}
    ot = []
    for t in range(Nt):
        if tm[t] in obs_data:
            ot.append(t)
            obs_at_t = np.array(obs_data[tm[t]])
            mfm.fit_to_data(E[t, :, :], obs_at_t)
            gammas.append(mfm.gamma)
            mu = mfm.predict_field(E[t,:,:])
            
            for obs in obs_at_t:
                i, j = obs.get_nearest_grid_point()
                sn = obs.get_station().get_name()
                rs = residuals[sn] if sn in residuals else []
                rs.append(obs.get_value() - E[t, i, j] * gammas[-1])
                residuals[sn] = rs
    
    gammas = np.array(gammas)
    ot = np.array(ot)

    # plot the fitteg gammas
    plt.figure()
    plt.plot(gammas)
    plt.title('MFM/Gamma values vs observation time')
    plt.savefig('results/mfm_gammas.png')
 
    # for each station, plot the equilibrium vs fm10 and accumulate residuals
    residuals = {}
    fm10_all = np.zeros((len(mtm), len(stations)))
    T_all = np.zeros_like(fm10_all)
    wfm10 = np.zeros_like(fm10_all)
    wT = np.zeros_like(fm10_all)
    for s, sndx in zip(stations, range(len(stations))):
        
        # compute equilibrium for station data
        H = np.array([o.get_value() for o in s.get_observations_for_times('rh', mtm)])
        T = np.array([o.get_value() for o in s.get_observations_for_times('air_temp', mtm)]) + 273.15
        equis = compute_model_equi_fm(H, T)

        # obtain fm10 station measurements 
        fm10so = [o.get_value() for o in s.get_observations_for_times('fm10', mtm)]
        fm10s = np.array(fm10so)
        Tso = [o.get_value() for o in s.get_observations_for_times('air_temp', mtm)]
        fm10_all[:, sndx] = fm10s
        T_all[:, sndx] = Tso 

        # obtain equilibrium at nearest grid point
        i, j = s.get_nearest_grid_point()
        fm10w = np.array(efms[wrf_tndx, 1, i, j])
        wfm10[:, sndx] = fm10w
        wT[:, sndx] = np.array(T2[wrf_tndx, i, j])
        equiw = np.array(wrf_equi[wrf_tndx, 1, i, j])
        equiw[equiw > 1.0] = 0.0
        plt.figure()
        plt.plot(mtm, equis, 'r+-', mtm, fm10s, 'rx-', mtm, equiw, 'g+-',mtm, fm10w, 'gx-')
        plt.title('fm10 and equi in WRF and in station %s' % s.get_name())
        plt.legend(['st. equi', 'st. fm10', 'WRF equi', 'WRF fm10'])
        plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%m'))
        for l in plt.gca().get_xticklabels():
            l.set_rotation(90)

        plt.savefig('results/%s_fm10_vs_equi.png' % s.get_name())
        
        residuals[s.get_name()] = fm10s - fm10w
 
    # **********************************************************************************
    # compute COVARIANCE between station residuals and plot this vs. distance
    Ns = len(stations)
    ss = [s.get_name() for s in stations]
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    WD = np.zeros((Ns,Ns))
    E = np.zeros((Ns,Ns))
    for i in range(Ns):
        r1, (lon1, lat1) = residuals[ss[i]], stations[i].get_position()
        i1,j1 = stations[i].get_nearest_grid_point()
        for j in range(Ns):
            r2, (lon2, lat2) = residuals[ss[j]], stations[j].get_position()
            i2,j2 = stations[j].get_nearest_grid_point()
            cc = np.cov(r1, r2)
            C[i,j] = cc[0, 1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)
            WD[i,j] = great_circle_distance(lon[i1,j1], lat[i1,j1], lon[i2,j2], lat[i2,j2])
            E[i,j] = np.abs(stations[i].get_elevation() - stations[j].get_elevation()) / 1000.0

    f = plt.figure(figsize = (16,16))
    f.subplots_adjust(hspace = 0.5, wspace = 0.5)
    ax = plt.subplot(221)
    plt.imshow(C, interpolation = 'nearest')
    plt.title('Covariance [-]')
    plt.colorbar()
    ax.set_xticks(np.arange(len(ss)))
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticks(np.arange(len(ss)))
    ax.set_yticklabels(ss)
    ax = plt.subplot(222)
    plt.imshow(D, interpolation = 'nearest')
    plt.title('Distance [km]')
    ax.set_xticks(np.arange(len(ss)))
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticks(np.arange(len(ss)))
    ax.set_yticklabels(ss)
    plt.colorbar()
    plt.subplot(223)
    colors = 'rgbcmyk'
    for r in range(D.shape[0]):
        ndx = np.setdiff1d(np.arange(Ns), [r])
        plt.plot(D[r,ndx], C[r,ndx], '%so' % colors[r%7])
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Covariance [-]')
    plt.title('Aggregate plot of covar vs. dist')

    plt.subplot(224)
    iu1 = np.triu_indices(Ns, 1)
    Dt = D[iu1]
    Ct = C[iu1]
    plt.plot(Dt, Ct, 'ro')
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Covariance [-]')
    plt.title('Aggregate plot of covar vs. distance')
    plt.savefig('results/station_residuals_covariance.png')

    with open('results/distance_matrix.txt', 'w') as f:
        for i in range(D.shape[0]):
            f.write(', '.join(map(str, D[i,:])))
            f.write('\n')
        f.close()

    with open('results/wdistance_matrix.txt', 'w') as f:
        for i in range(WD.shape[0]):
            f.write(', '.join(map(str, WD[i,:])))
            f.write('\n')
        f.close()

    with open('results/fm10_data.txt', 'w') as f:
        for i in range(len(mtm)):
            f.write(', '.join(map(str, fm10_all[i,:])))
            f.write('\n')
        f.close()
    
    with open('results/T_data.txt', 'w') as f:
        for i in range(len(mtm)):
            f.write(', '.join(map(str, T_all[i,:])))
            f.write('\n')
        f.close()

    with open('results/wfm10_data.txt', 'w') as f:
        for i in range(len(mtm)):
            f.write(','.join(map(str, wfm10[i,:])))
            f.write('\n')
        f.close()

    with open('results/wT_data.txt', 'w') as f:
        for i in range(len(mtm)):
            f.write(', '.join(map(str, wT[i,:])))
            f.write('\n')
        f.close()

    with open('results/elevations.txt', 'w') as f:
        for i in range(len(E)):
            f.write(', '.join(map(str, E[i,:])))
            f.write('\n')
        f.close()

            

    # plot the covariance vs. distance from the viewpoint of each station
    for i in range(len(stations)):
        n = stations[i].get_name()
        ci = C[i,:].copy()
        di = D[i,:].copy()
        ndx = np.argsort(di)
        ci, di = ci[ndx], di[ndx]
        plt.figure()
        plt.plot(di, ci, 'o-')
        plt.title('Distance vs. cov for station %s' % n)
        plt.xlabel('Distance [km]')
        plt.ylabel('Covariance [-]')
        plt.savefig('results/%s_cov_vs_dist.png' % n)


    # **********************************************************************************
    # compute CORRELATION COEFFICIENT between station residuals and plot this vs. distance
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    E = np.zeros((Ns,Ns))
    for i in range(Ns):
        r1, (lon1, lat1) = residuals[ss[i]], stations[i].get_position()
        for j in range(Ns):
            s2name = stations[j].get_name()
            r2, (lon2, lat2) = residuals[ss[j]], stations[j].get_position()
            cc = np.corrcoef(r1, r2)
            C[i,j] = cc[0, 1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)
            E[i,j] = np.abs(stations[i].get_elevation() - stations[j].get_elevation()) / 1000.0

    f = plt.figure(figsize = (16,16))
    f.subplots_adjust(hspace = 0.5, wspace = 0.5)
    ax = plt.subplot(221)
    plt.imshow(C, interpolation = 'nearest')
    plt.clim([0.0, 1.0])
    plt.title('Correlation coefficient [-]')
    plt.colorbar()
    ax.set_xticks(np.arange(len(ss)))
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticks(np.arange(len(ss)))
    ax.set_yticklabels(ss)
    ax = plt.subplot(222)
    plt.imshow(D, interpolation = 'nearest')
    plt.title('Distance [km]')
    ax.set_xticks(np.arange(len(ss)))
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticks(np.arange(len(ss)))
    ax.set_yticklabels(ss)
    plt.colorbar()
    plt.subplot(223)
    colors = 'rgbcmyk'
    for r in range(D.shape[0]):
        ndx = np.setdiff1d(np.arange(Ns), [r])
        plt.plot(D[r,ndx], C[r,ndx], '%so' % colors[r%7])
    plt.ylim([0.0, 1.0])
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation coefficient [-]')
    plt.title('Aggregate plot of cc vs. dist')

    plt.subplot(224)
    iu1 = np.triu_indices(Ns, 1)
    Et = E[iu1]
    Dt = D[iu1]
    Ct = C[iu1]
    beta, Rsq = compute_ols_estimator(np.hstack([Dt[:,np.newaxis], np.ones_like(Dt[:,np.newaxis])]), Ct[:,np.newaxis])
    print beta, Rsq
    olsDT = beta[1] + beta[0] * Dt
    plt.plot(Dt, Ct, 'ro')
    plt.plot(Dt, olsDT, 'k-')
    plt.ylim([0.0, 1.0])
    plt.legend(['Data', 'OLS fit'], loc = 'lower left', prop = { 'size' : 10 } )
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation coefficient [-]')
    plt.title('Aggregate plot of cc vs. distance')
    plt.savefig('results/station_residuals_correlation.png')
    
    # plot the correlation coefficient vs. distance for each station
    for i in range(len(stations)):
        n = stations[i].get_name()
        ci = C[i,:].copy()
        di = D[i,:].copy()
        ndx = np.argsort(di)
        ci, di = ci[ndx], di[ndx]
        plt.figure()
        plt.plot(di, ci, 'o-')
        plt.title('Distance vs. corr. coef. for station %s' % n)
        plt.xlabel('Distance [km]')
        plt.ylabel('Covariance [-]')
        plt.savefig('results/%s_cc_vs_dist.png' % n)
    
