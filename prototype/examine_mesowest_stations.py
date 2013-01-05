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
from time_series_utilities import match_time_series, build_observation_data
from spatial_model_utilities import render_spatial_field, great_circle_distance
                                    
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


def compute_model_equi_fm(station):
    """
    Compute the equilibrium temperature given by the
    moisture model from the relative humidity measured at the station
    and the temperature measured at the station.
    """
    H = np.array(station.get_observations('rh'))
    T = np.array(station.get_observations('air_temp')) + 273.15

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

    # for each station id, load the station
    stations = []
    for sinfo in si_list:
        code = sinfo.split(',')[0]
        mws = MesoWestStation(sinfo)
        for suffix in [ '_1', '_2', '_3', '_4', '_5', '_6', '_7' ]:
            mws.load_station_data(os.path.join(station_data_dir, '%s%s.xls' % (code, suffix)))
        stations.append(mws)

    print('Loaded %d stations.' % len(stations))

    # check stations for nans
    stations = filter(MesoWestStation.data_ok, stations)
    print('Have %d stations with complete data.' % len(stations))

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
#    render_spatial_field(m, lon, lat, E[0,:,:], 'Equilibrium moisture')
    plt.figure()
    m.drawparallels(np.arange(llcrnrlat,urcrnrlat,1.0))
    m.drawmeridians(np.arange(llcrnrlon,urcrnrlon,1.0))
    for s in stations:
        slon, slat = s.get_position()
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')
        plt.text(x, y, s.get_name(), color = 'white')

    # for each station, plot the equilibrium vs fm10 and accumulate residuals
    residuals = {}
    for s in stations:
        equi = compute_model_equi_fm(s)
        fm10 = np.array(s.get_observations('fm10')) / 100.0
        residuals[s.get_name()] = fm10 - equi
        plt.figure()
        plt.plot(equi, 'bo-')
        plt.plot(fm10, 'ro-')
        plt.legend(['equi', 'fm10'])
        plt.title('Model equi vs. fm10 for %s' % s.get_name())
        plt.savefig('results/%s_fm10_vs_equi.png' % s.get_name())
 
    # **********************************************************************************
    # compute COVARIANCE between station residuals and plot this vs. distance
    Ns = len(stations)
    ss = [s.get_name() for s in stations]
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    E = np.zeros((Ns,Ns))
    for i in range(Ns):
        r1, (lon1, lat1) = residuals[ss[i]], stations[i].get_position()
        for j in range(Ns):
            r2, (lon2, lat2) = residuals[ss[j]], stations[j].get_position()
            cc = np.cov(r1, r2)
            C[i,j] = cc[0, 1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)
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

#    f = open('distance_vs_covariance.txt', 'w')
#    for i in range(Dt.shape[0]):
#        f.write('%g,%g\n' % (Dt[i], Ct[i]))
#    f.close() 


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
    
    
