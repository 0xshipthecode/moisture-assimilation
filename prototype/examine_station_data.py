
from time_series_utilities import match_time_series, build_observation_data
from spatial_model_utilities import render_spatial_field, great_circle_distance
                                    
from wrf_model_data import WRFModelData
from observation_stations import MesoWestStation
from diagnostics import init_diagnostics, diagnostics

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
import pytz
import os
import string


if __name__ == '__main__':
    
    init_diagnostics('results/examine_station_data.log')

    # hack to use configuration object without option passing
    cfg = { 'station_data_dir' : '../real_data/colorado_stations/',
            'station_list_file' : 'clean_stations',
            'wrf_data_file' : '../real_data/colorado_stations/wrfout_sel_5km_10mhist.nc',
            'assimilation_window' : 3600,
            'output_dir' : 'results',
            'max_dist' : 200.0,
            'bin_width' : 20
          }

    # load the smallest domain
    wrf_data = WRFModelData(cfg['wrf_data_file'], ['T2', 'PSFC', 'Q2', 'HGT'])

    # read in vars
    lon, lat = wrf_data.get_lons(), wrf_data.get_lats() 
    hgt = wrf_data['HGT']
    tm = wrf_data.get_gmt_times()
    Nt = len(tm)

    print("Loaded %d timestamps from WRF." % Nt)
        
    # load station data from files
    with open(os.path.join(cfg['station_data_dir'], cfg['station_list_file']), 'r') as f:
        si_list = f.read().split('\n')

    si_list = filter(lambda x: len(x) > 0 and x[0] != '#', map(string.strip, si_list))

    # for each station id, load the station
    stations = []
    for code in si_list:
        mws = MesoWestStation(code)
        mws.load_station_info(os.path.join(cfg["station_data_dir"], "%s.info" % code))
        mws.register_to_grid(wrf_data)
        mws.load_station_data(os.path.join(cfg["station_data_dir"], "%s.obs" % code))
        stations.append(mws)

    print('Loaded %d stations.' % len(stations))
    
    # construct basemap for rendering
    domain_rng = wrf_data.get_domain_extent()
    m = Basemap(llcrnrlon=domain_rng[0],llcrnrlat=domain_rng[1],
                urcrnrlon=domain_rng[2],urcrnrlat=domain_rng[3],
                projection = 'mill')

    # show the equilibrium field and render position of stations on top
    render_spatial_field(m, lon, lat, hgt[0,:,:], 'Station positions over topography')
    for s in stations:
        slon, slat = s.get_position()
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')
        plt.text(x, y, s.get_name(), color = 'white')

    x, y = m(lon.ravel(), lat.ravel())
    plt.plot(x, y, 'k+', markersize = 4)
    plt.savefig(os.path.join(cfg['output_dir'], 'station_positions.png'))
        
    # part C - compute the variogram estimator for each assimilation window
    obs_data = build_observation_data(stations, 'FM')
    print("Found a total of %d observations." % len(obs_data))
    assim_win = cfg['assimilation_window']
    max_dist = cfg['max_dist']
    bin_width = cfg['bin_width']
    plt.figure(figsize = (10, 5))

    all_dists = []
    all_sqdiffs = []
    for t in range(0, Nt, 10):

        # find number of observations valid now (depends on assimilation window)
        t_now = tm[t]
        obs_at_t = []
        for k, v in obs_data.iteritems():
            if abs((k-t_now).total_seconds()) < assim_win:
                obs_at_t.extend(v)

        if len(obs_at_t) > 0:
            # construct pairwise dataset
            dists = []
            sqdiffs = []
            for i in obs_at_t:
                pi = i.get_position()
                oi = i.get_value()
                for j in obs_at_t:
                    pj = j.get_position()
                    oj = j.get_value()
                    if i != j:
                        dist = great_circle_distance(pi[0], pi[1], pj[0], pj[1])
                        if dist < max_dist:
                            dists.append(dist)
                            sqdiffs.append(0.5 * (oi - oj)**2)
                            all_dists.append(dist)
                            all_sqdiffs.append(0.5 * (oi - oj)**2)

            # compute an empirical variogram
            bsqdiffs = []
            for lim in np.arange(bin_width, max_dist + bin_width, bin_width):
                ndx = [i for i in range(len(dists)) if (lim-bin_width <= dists[i]) and (dists[i] < lim)]
                bsqdiffs.append(np.mean([sqdiffs[i] for i in ndx]))

            # construct a plot at this time
            plt.clf()
            plt.plot(dists, sqdiffs, 'ro', markersize = 5)
            plt.plot(np.arange(bin_width, max_dist+bin_width, bin_width) - bin_width/2.0, bsqdiffs, 'b-', linewidth = 2.0)
            plt.xlabel('Distance [km]')
            plt.ylabel('Squared obs. diff [-]')
            plt.title('Variogram at time %s' % str(t_now))
            plt.savefig(os.path.join(cfg['output_dir'], 'variogram_est_%03d.png' % t))


    # compute an empirical variogram
    bsqdiffs = []
    for lim in np.arange(bin_width, max_dist + bin_width, bin_width):
        ndx = [i for i in range(len(all_dists)) if (lim-bin_width <= all_dists[i]) and (all_dists[i] < lim)]
        bsqdiffs.append(np.mean([all_sqdiffs[i] for i in ndx]))

    # construct a plot at this time
    plt.clf()
    plt.plot(all_dists, all_sqdiffs, 'ro', markersize = 5)
    plt.plot(np.arange(bin_width, max_dist+bin_width, bin_width) - bin_width/2.0, bsqdiffs, 'b-', linewidth = 2.0)
    plt.xlabel('Distance [km]')
    plt.ylabel('Squared obs. diff [-]')
    plt.title('Variogram for all times')
    plt.savefig(os.path.join(cfg['output_dir'], 'variogram_est_all.png'))

    a = plt.axis()
    plt.axis([a[0], a[1], a[2], 0.01])
    plt.savefig(os.path.join(cfg['output_dir'], 'variogram_est_all_censored.png'))
