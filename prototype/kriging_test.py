

from spatial_model_utilities import render_spatial_field

from kriging_methods import simple_kriging_data_to_model
from time_series_utilities import build_observation_data
from wrf_model_data import WRFModelData
from mean_field_model import MeanFieldModel
from diagnostics import init_diagnostics, diagnostics

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pytz
from observation_stations import Station
import os


station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Moisture",
                  "Alpine_Moisture",
                  "Ranchita_Moisture",
                  "Camp_Elliot_Moisture",
                ]

                  
station_data_dir = "../real_data/witch_creek/"

    
def run_module():

    # configure diagnostics        
    init_diagnostics("results/kriging_test_diagnostics.txt")
    diagnostics().configure_tag("skdm_obs_res", True, True, True)
    diagnostics().configure_tag("skdm_obs_res_mean", True, True, True)
        
    wrf_data = WRFModelData('../real_data/witch_creek/realfire03_d04_20071022.nc')
    
    # read in vars
    lat, lon = wrf_data.get_lats(), wrf_data.get_lons()
    tm = wrf_data.get_times()
    rain = wrf_data['RAINNC']
    Ed, Ew = wrf_data.get_moisture_equilibria()
    
    # obtain sizes
    Nt = rain.shape[0]
    dom_shape = lat.shape
    locs = np.prod(dom_shape)
    
    # load station data, match to grid points and build observation records
    # load station data from files
    tz = pytz.timezone('US/Pacific')
    stations = [Station(os.path.join(station_data_dir, s), tz, wrf_data) for s in station_list]
    obs_data = build_observation_data(stations, 'fuel_moisture', wrf_data) 
    
    # construct initial vector
    mfm = MeanFieldModel()
    
    # set up parameters
    mod_res_std = np.ones_like(Ed[0,:,:]) * 0.05
    obs_res_std = np.ones((len(stations),)) * 0.1
    
    # construct a basemap representation of the area
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    plt.figure(figsize = (10, 6))
    
    # run model
    ndx = 1
    for t in range(1, Nt):
        model_time = wrf_data.get_times()[t]
        E = 0.5 * (Ed[t,:,:] + Ew[t,:,:])

        # if we have an observation somewhere in time, run kriging
        if model_time in obs_data:
            print("Time: %s, step: %d" % (str(model_time), t))
            
            mfm.fit_to_data(E, obs_data[model_time])
            Efit = mfm.predict_field(E)

            # krige data to observations
            K, V = simple_kriging_data_to_model(obs_data[model_time], obs_res_std, Efit, wrf_data, mod_res_std, t)
                
            plt.clf()
            plt.subplot(2,2,1)
            render_spatial_field(m, lon, lat, Efit, 'Equilibrium')
            plt.clim([0.0, 0.2])
            plt.colorbar()

            plt.subplot(2,2,2)
            render_spatial_field(m, lon, lat, K, 'Kriging field')
            plt.clim([0.0, 0.2])
            plt.colorbar()

            plt.subplot(2,2,3)
            render_spatial_field(m, lon, lat, V, 'Kriging variance')
            plt.clim([0.0, np.max(V)])
            plt.colorbar()
            
            plt.subplot(2,2,4)
            render_spatial_field(m, lon, lat, K - Efit, 'Kriging vs. mean field residuals')
#            plt.clim([0.0, np.max()])
            plt.colorbar()
            
            plt.savefig('model_outputs/kriging_test_t%03d.png' % (ndx))
            ndx += 1 
        
                
    
if __name__ == '__main__':
#    profile.run('run_module(); print', 'spatial_model.stats')
#    
#    stats = pstats.Stats('spatial_model.stats')
#    stats.strip_dirs()
#    stats.sort_stats('cumulative')
#    stats.print_stats()

    run_module()
