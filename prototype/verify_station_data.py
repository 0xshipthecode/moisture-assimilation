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


def check_zeros(mws):
    """
    Check if the data are not all zero.
    """
    d = np.array(mws.get_observations_raw('fm10'))
    return np.sum(d==0) < 5


station_data_dir = "../real_data/colorado_stations/"

if __name__ == '__main__':
    
    init_diagnostics('results/verify_mesowest_station_data.log')

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

    # check stations for zeros and for non-changing data
    stations = filter(check_zeros, stations)
    print('Have %d stations with non-zeros in data.' % len(stations))

    # setup measurement variances
    print('Stats for stations that have passed the tests.')
    for st in stations:
        d = st.get_observations_raw('fm10')
        print('%s: fm10 mean (%g)  std(%g)' % (st.get_name(), np.mean(d), np.std(d)))

    # write the clean data to file
    with open('clean_stations', 'w') as f:
        for st in stations:
            n = st.get_name()
            for si in si_list:
                if si.startswith(n):
                    f.write(si)
                    f.write('\n')

