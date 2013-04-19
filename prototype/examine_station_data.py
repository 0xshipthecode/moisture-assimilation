
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


def plot_variogram(dists, sqdiffs, bins, title, savename):
    """
    Plots the various (halved) squared differences as a point plot and overlays
    a curve, which is the empirical variogram estimate.
    """
    bin_width = bins[1] - bins[0]

    # compute an empirical variogram
    bsqdiffs = []
    for lim in bins:
        ndx = [i for i in range(len(dists)) if (lim - bin_width <= dists[i]) and (dists[i] < lim)]
        bsqdiffs.append(np.mean([sqdiffs[i] for i in ndx]))

    # construct a plot at this time
    plt.clf()
    plt.plot(dists, sqdiffs, 'ro', markersize = 5)
    plt.plot(np.arange(bin_width, max_dist+bin_width, bin_width) - bin_width/2.0, bsqdiffs, 'b-', linewidth = 2.0)
    plt.xlabel('Distance [km]')
    plt.ylabel('Squared obs. diff [-]')
    plt.title('Variogram at time %s' % str(t_now))
    plt.savefig(savename)



if __name__ == '__main__':
    
    init_diagnostics('results/examine_station_data.log')

    what = [ 'variograms', 'distributions' ]

    # hack to use configuration object without option passing
    cfg = { 'station_info_dir' : '../real_data/colorado_stations',
            'station_obs_dir' : '../real_data/colorado_stations/20120601',
            'station_list_file' : '../real_data/colorado_stations/clean_stations',
            'wrf_data_file' : '../real_data/colorado_stations/wrf20120601_sel_5km_10mhist.nc',
            'assimilation_window' : 3600,
            'output_dir' : 'results',
            'max_dist' : 100.0,
            'bin_width' : 10,
            'standardize' : True
          }

    # cfg = { 'station_info_dir' : '../real_data/witch_creek',
    #         'station_obs_dir' : '../real_data/witch_creek',
    #         'station_list_file' : '../real_data/witch_creek/station_list',
    #         'wrf_data_file' : '../real_data/witch_creek/wrf20071021_witchcreek_all.nc',
    #         'assimilation_window' : 3600,
    #         'output_dir' : 'results_wc',
    #         'max_dist' : 200.0,
    #         'bin_width' : 20,
    #         'standardize' : False
    #       }


    # load the smallest domain
    wrf_data = WRFModelData(cfg['wrf_data_file'], ['T2', 'PSFC', 'Q2', 'HGT'])

    # read in vars
    lon, lat = wrf_data.get_lons(), wrf_data.get_lats() 
    hgt = wrf_data['HGT']
    tm = wrf_data.get_gmt_times()
    Nt = len(tm)

    print("Loaded %d timestamps from WRF." % Nt)
        
    # load station data from files
    with open(cfg['station_list_file'], 'r') as f:
        si_list = f.read().split('\n')

    si_list = filter(lambda x: len(x) > 0 and x[0] != '#', map(string.strip, si_list))

    # for each station id, load the station
    stations = []
    for code in si_list:
        mws = MesoWestStation(code)
        mws.load_station_info(os.path.join(cfg["station_info_dir"], "%s.info" % code))
        mws.register_to_grid(wrf_data)
        mws.load_station_data(os.path.join(cfg["station_obs_dir"], "%s.obs" % code))
        stations.append(mws)

    print('Loaded %d stations.' % len(stations))
    
    # construct basemap for rendering
    domain_rng = wrf_data.get_domain_extent()
    m = Basemap(llcrnrlon=domain_rng[0],llcrnrlat=domain_rng[1],
                urcrnrlon=domain_rng[2],urcrnrlat=domain_rng[3],
                projection = 'mill')

    # show the equilibrium field and render position of stations on top
    plt.figure(figsize = (12,9))
    render_spatial_field(m, lon, lat, hgt[0,:,:], 'Station positions over topography')
    for s in stations:
        slon, slat = s.get_position()
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 4, markerfacecolor = 'magenta')
        plt.text(x, y, s.get_id(), color = 'black')


    # compute station means and standard deviations
    station_stats = {}
    for s in stations:
        v = [o.get_value() for o in s.get_observations("FM")]
        stdev = np.std(v)
        if np.isnan(stdev) or abs(stdev) <  1e-6:
            stdev = 1.0
        station_stats[s.get_id()] = (np.mean(v), stdev)

#    x, y = m(lon.ravel(), lat.ravel())
#    plt.plot(x, y, 'k+', markersize = 4)
    plt.savefig(os.path.join(cfg['output_dir'], 'station_positions.png'))

    if 'variograms' in what:
        
        # part A - compute the variogram estimator for each assimilation window
        obs_data = build_observation_data(stations, 'FM')
        print("Found a total of %d observations." % len(obs_data))
        assim_win = cfg['assimilation_window']
        max_dist = cfg['max_dist']
        bin_width = cfg['bin_width']
        plt.figure(figsize = (10, 5))

        all_dists = []
        all_sqdiffs = []
        bins = np.arange(bin_width, max_dist + bin_width, bin_width)
        for t in range(0, Nt, 10):

            # find number of observations valid now (depends on assimilation window)
            t_now = tm[t]
            obs_at_t = []
            for k, v in obs_data.iteritems():
                if abs((k-t_now).total_seconds()) < assim_win / 2.0:
                    obs_at_t.extend(v)

            if len(obs_at_t) > 0:
                # construct pairwise dataset
                dists = []
                sqdiffs = []
                for i in range(len(obs_at_t)):
                    sid = obs_at_t[i].get_station().get_id()
                    pi = obs_at_t[i].get_position()
                    stats = station_stats[sid]
                    oi = obs_at_t[i].get_value()
                    if cfg['standardize']:
                        oi = (oi - stats[0]) / stats[1]
                    for j in range(i+1, len(obs_at_t)):
                        pj = obs_at_t[j].get_position()
                        oj = obs_at_t[j].get_value()
                        dist = great_circle_distance(pi[0], pi[1], pj[0], pj[1])
                        if dist < max_dist:
                            dists.append(dist)
                            sqdiffs.append(0.5 * (oi - oj)**2)
                            all_dists.append(dist)
                            all_sqdiffs.append(0.5 * (oi - oj)**2)

                fname = 'std_variogram_est_%03d.png' % t if cfg['standardize'] else 'variogram_est_%03d.png' % t
                plot_variogram(dists, sqdiffs, bins, 'Variogram at time %s' % str(t_now), 
                               os.path.join(cfg['output_dir'], fname))

        fname = 'std_variogram_est_all.png' if cfg['standardize'] else 'variogram_est_all.png'
        plot_variogram(all_dists, all_sqdiffs, bins, 'Variogram (all observations)',
                       os.path.join(cfg['output_dir'], fname))

        # hack to plot the lower part of the variogram
        a = plt.axis()
        plt.axis([a[0], a[1], a[2], 0.1 * max(all_sqdiffs)])
        fname = 'std_variogram_est_all_censored.png' if cfg['standardize'] else 'variogram_est_all_censored.png'
        plt.savefig(os.path.join(cfg['output_dir'], fname))

    
    # part B (obtain the distribution of FM observations from each station and plot it)
    if 'distributions' in what:
        plt.figure(figsize = (8, 4))
        for s in stations:
            obs = [o.get_value() for o in s.get_observations("FM")]
            val_freq = {}
            for v in obs:
                a = val_freq.get(v, 0)
                val_freq[v] = a+1

            vals = np.array(sorted(val_freq.keys()))
            freqs = [val_freq[v] for v in vals]
        
            plt.clf()
            plt.bar(vals - 0.01/4, freqs, width = 0.01 / 2)
            plt.title("Distribution of observations at %s" % s.get_id())
            plt.savefig(os.path.join(cfg['output_dir'], 'distr_%s.png' % s.get_id()))
