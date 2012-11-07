
from time_series_utilities import match_time_series, build_observation_data
from spatial_model_utilities import render_spatial_field, great_circle_distance
                                    
from wrf_model_data import WRFModelData
from observation_stations import Station, Observation
from mean_field_model import MeanFieldModel

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
import pytz
import os


station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Moisture",
                  "Alpine_Moisture",
                  "Ranchita_Moisture",
                  "Camp_Elliot_Moisture"
                ]
                  
station_data_dir = "../real_data/witch_creek/"


def plot_stations_vs_model_ts(stations, field_name, field, wrf_data):
    
    # extract station time series and corresponding grid point time series
    f = plt.figure()
    f.subplots_adjust(hspace = 1.2)
    sp = 1
    for s in stations:
        i, j = s.get_nearest_grid_point()
        model_tm = wrf_data.get_times()
        mindx, obs_list = s.get_observations_for_times(field_name, model_tm)
        obs_ts = [obs.get_value() for obs in obs_list]
        common_tm = [model_tm[t] for t in mindx]
        ax = plt.subplot(4, 2, sp)
        ax.xaxis.set_major_formatter(DateFormatter('%H:%m')) 
        plt.plot(common_tm, obs_ts, 'ro-', common_tm, field[:, i, j], 'gx-', linewidth = 2)
        plt.title('%s vs. model %s' % (s.get_name(), field_name))
        for l in ax.get_xticklabels():
            l.set_rotation(90)
        sp += 1



if __name__ == '__main__':

    # load the smallest domain
    wrf_data = WRFModelData('../real_data/witch_creek/realfire03_d04_20071021.nc', tz_name = 'US/Pacific')

    # read in vars
    lon, lat = wrf_data.get_lons(), wrf_data.get_lats() 
    P, Q, T, rain = wrf_data['PSFC'], wrf_data['Q2'], wrf_data['T2'], wrf_data['RAINNC']
    tm = wrf_data.get_times()
    Nt = len(tm)
        
    # load stations and match then to grid points
    # load station data from files
    tz = pytz.timezone('US/Pacific')
    stations = [Station(os.path.join(station_data_dir, s), tz, wrf_data) for s in station_list]
    
    # construct basemap for rendering
    domain_rng = wrf_data.get_domain_extent()
    m = Basemap(llcrnrlon=domain_rng[0],llcrnrlat=domain_rng[1],
                urcrnrlon=domain_rng[2],urcrnrlat=domain_rng[3],
                projection = 'mill')

    # compute the equilibrium moisture on grid points (for all times t)
    Ed, Ew = wrf_data.get_moisture_equilibria()
    E = 0.5 * (Ed + Ew)
    
    mfm = MeanFieldModel()

    # show the equilibrium field and render position of stations on top
    render_spatial_field(m, lon, lat, E[0,:,:], 'Equilibrium moisture')
    for s in stations:
        slon, slat = s.get_position()
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')
        plt.text(x, y, s.get_name())

    x, y = m(lon.ravel(), lat.ravel())
    plt.plot(x, y, 'k+', markersize = 4)
        
    # part C - fit the mean field model
    obs_data = build_observation_data(stations, 'fuel_moisture', wrf_data)
    gammas = []
    residuals = {}
    ot = []
    for t in range(Nt):
        if tm[t] in obs_data:
            ot.append(t)
            obs_at_t = obs_data[tm[t]] 
            mfm.fit_to_data(E[t, :, :], obs_at_t)
            gammas.append(mfm.gamma)
            mu = mfm.predict_field()
            
            for obs in obs_at_t:
                i, j = obs.get_nearest_grid_point()
                sn = obs.get_station().get_name()
                rs = residuals[sn] if sn in residuals else []
                rs.append(obs.get_value() - E[t, i, j])
                residuals[sn] = rs
    
    gammas = np.array(gammas)
    ot = np.array(ot)

    # part B, compare values at the same times
    for st_field_name, field in [ ('air_temp', T[ot,:,:] - 273.15), ('fuel_moisture', E[ot,:,:])]:
        plot_stations_vs_model_ts(stations, st_field_name, field, wrf_data)

    # plot fitted model
    plot_stations_vs_model_ts(stations, 'fuel_moisture', gammas[:, np.newaxis, np.newaxis] * E[ot,:,:], wrf_data)
    
    # **********************************************************************************
    # compute COVARIANCE between station residuals and plot this vs. distance
    Ns = len(stations)
    ss = [s.get_name() for s in stations]
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    for i in range(Ns):
        r1, (lon1, lat1) = residuals[ss[i]], stations[i].get_position()
        for j in range(Ns):
            s2name = stations[j].get_name()
            r2, (lon2, lat2) = residuals[ss[j]], stations[j].get_position()
            cc = np.cov(r1, r2)
            C[i,j] = cc[0, 1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)

    f = plt.figure(figsize = (8,8))
    f.subplots_adjust(hspace = 0.5, wspace = 0.5)
    ax = plt.subplot(221)
    plt.imshow(C, interpolation = 'nearest')
    plt.title('Covariance [-]')
    plt.colorbar()
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticklabels(ss)
    ax = plt.subplot(222)
    plt.imshow(D, interpolation = 'nearest')
    plt.title('Distance [km]')
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticklabels(ss)
    plt.colorbar()
    plt.subplot(223)
    colors = 'rgbcmyk'
    for r in range(D.shape[0]):
        ndx = np.setdiff1d(np.arange(Ns), [r])
        plt.plot(D[r,ndx], C[r,ndx], '%so' % colors[r])
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

    f = open('distance_vs_covariance.txt', 'w')
    for i in range(Dt.shape[0]):
        f.write('%g,%g\n' % (Dt[i], Ct[i]))
    f.close() 


    # **********************************************************************************
    # compute CORRELATION COEFFICIENT between station residuals and plot this vs. distance
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    for i in range(Ns):
        r1, (lon1, lat1) = residuals[ss[i]], stations[i].get_position()
        for j in range(Ns):
            s2name = stations[j].get_name()
            r2, (lon2, lat2) = residuals[ss[j]], stations[j].get_position()
            cc = np.corrcoef(r1, r2)
            C[i,j] = cc[0, 1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)

    f = plt.figure(figsize = (8,8))
    f.subplots_adjust(hspace = 0.5, wspace = 0.5)
    ax = plt.subplot(221)
    plt.imshow(C, interpolation = 'nearest')
    plt.clim([0.0, 1.0])
    plt.title('Correlation coefficient [-]')
    plt.colorbar()
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticklabels(ss)
    ax = plt.subplot(222)
    plt.imshow(D, interpolation = 'nearest')
    plt.title('Distance [km]')
    ax.set_xticklabels(ss, rotation = 90)
    ax.set_yticklabels(ss)
    plt.colorbar()
    plt.subplot(223)
    colors = 'rgbcmyk'
    for r in range(D.shape[0]):
        ndx = np.setdiff1d(np.arange(Ns), [r])
        plt.plot(D[r,ndx], C[r,ndx], '%so' % colors[r])
    plt.ylim([0.0, 1.0])
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation coefficient [-]')
    plt.title('Aggregate plot of cc vs. dist')

    plt.subplot(224)
    iu1 = np.triu_indices(Ns, 1)
    Dt = D[iu1]
    Ct = C[iu1]
    plt.plot(Dt, Ct, 'ro')
    plt.plot(Dt, 0.8166 - 0.0052 * Dt, 'r-')
    plt.plot(Dt, 0.8565 - 0.0063 * Dt, 'g-')
    plt.plot(Dt, 0.8729 - 0.0068 * Dt, 'b-')
    plt.ylim([0.0, 1.0])
    plt.legend(['Data', 'OLS', 'Huber', 'Bisquare'], loc = 'lower left', prop = { 'size' : 10 } )
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation coefficient [-]')
    plt.title('Aggregate plot of cc vs. distance')
    
    f = open('distance_vs_covariance.txt', 'w')
    for i in range(Dt.shape[0]):
        f.write('%g,%g\n' % (Dt[i], Ct[i]))
    f.close() 
    
    plt.show()
    