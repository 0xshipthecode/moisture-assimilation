# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:36 2012

@author: martin
"""

from spatial_model_utilities import render_spatial_field_fast, great_circle_distance
from time_series_utilities import build_observation_data

from kriging_methods import simple_kriging_data_to_model
from wrf_model_data import WRFModelData
#from cell_model import CellMoistureModel
from cell_model_opt import CellMoistureModel
from mean_field_model import MeanFieldModel
from observation_stations import Station
from diagnostics import init_diagnostics, diagnostics
from online_variance_estimator import OnlineVarianceEstimator

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pytz
import cPickle


station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Moisture",
                  "Alpine_Moisture",
                  "Ranchita_Moisture",
                  "Camp_Elliot_Moisture"
                ]


station_data_dir = "../real_data/witch_creek/"


def run_module():
    
    # read in configuration file to execute run
    print("Reading configuration from [%s]" % sys.argv[1])
    
    with open(sys.argv[1]) as f:
        cfg = eval(f.read())
    
    # ensure output path exists
    if not os.path.isdir(cfg['output_dir']): 
        os.mkdir(cfg['output_dir'])
        
    # configure diagnostics        
    init_diagnostics(os.path.join(cfg['output_dir'], 'moisture_model_v1_diagnostics.txt'))
    diagnostics().configure_tag("skdm_obs_res", False, True, True)
    diagnostics().configure_tag("skdm_cov_cond", False, True, True)

    diagnostics().configure_tag("assim_mV", False, False, True)
    diagnostics().configure_tag("assim_K", False, False, True)
    diagnostics().configure_tag("assim_data", False, False, True)
    diagnostics().configure_tag("assim_mresV", False, False, True)

    diagnostics().configure_tag("kriging_variance", False, False, True)
    diagnostics().configure_tag("kriging_obs_res_var", False, False, True)
    
    wrf_data = WRFModelData(cfg['input_file'], tz_name = 'US/Pacific')
    
    # read in vars
    lat, lon = wrf_data.get_lats(), wrf_data.get_lons()
    tm = wrf_data.get_times()
    rain = wrf_data['RAINNC']
    Ed, Ew = wrf_data.get_moisture_equilibria()
    
    # find maximum moisture overall to set up visualization
#    maxE = max(np.max(Ed), np.max(Ew)) * 1.2
    maxE = 0.3
    
    # obtain sizes
    Nt = rain.shape[0]
    dom_shape = lat.shape
    
    # load station data from files
    tz = pytz.timezone('US/Pacific')
    stations = [Station(os.path.join(station_data_dir, s), tz, wrf_data) for s in station_list]
    for s in stations:
        s.set_measurement_variance('fm10', 0.1)
    
    # build the observation data structure indexed by time
    obs_data_fm10 = build_observation_data(stations, 'fm10', wrf_data)
    
    # construct initial conditions
    E = 0.5 * (Ed[1,:,:] + Ew[1,:,:])
    
    # set up parameters
    Q = np.eye(9) * 0.0001
    P0 = np.eye(9) * 0.001
    dt = 10.0 * 60
    K = np.zeros_like(E)
    V = np.zeros_like(E)
    mV = np.zeros_like(E)
    Kg = np.zeros_like(E)
    predicted_field = np.zeros_like(E)
    mresV = np.zeros_like(E)
    Kf_fn = np.zeros_like(E)
    Vf_fn = np.zeros_like(E)
    mid = np.zeros_like(E)
    Kg = np.zeros((dom_shape[0], dom_shape[1], 1))
    
    # moisture state and observation residual variance estimators
    mod_re = OnlineVarianceEstimator(np.zeros_like(E), np.ones_like(E) * 0.03, 1)
    obs_re = OnlineVarianceEstimator(np.zeros((len(stations),)), np.ones(len(stations),) * 0.1, 1)
    
    # initialize the mean field model (default fit is 1.0 of equilibrium before new information comes in)
    mfm = MeanFieldModel(cfg['lock_gamma'])

    # construct model grid using standard fuel parameters
    Tk = np.array([1.0, 10.0, 100.0]) * 3600
    models = np.zeros(dom_shape, dtype = np.object)
    models_na = np.zeros_like(models)
    for pos in np.ndindex(dom_shape): 
        models[pos] = CellMoistureModel((lat[pos], lon[pos]), 3, E[pos], Tk, P0 = P0)
        models_na[pos] = CellMoistureModel((lat[pos], lon[pos]), 3, E[pos], Tk, P0 = P0)

    m = None

    plt.figure(figsize = (12, 8))
    
    # run model
    for t in range(1, Nt):
        model_time = wrf_data.get_times()[t]
        print("Time: %s, step: %d" % (str(model_time), t))

        # pre-compute equilibrium moisture to save a lot of time
        E = 0.5 * (Ed[t,:,:] + Ew[t,:,:])
        
        # run the model update
        for pos in np.ndindex(dom_shape):
            i, j = pos
            models[pos].advance_model(Ed[t, i, j], Ew[t, i, j], rain[t, i, j], dt, Q)
            models_na[pos].advance_model(Ed[t, i, j], Ew[t, i, j], rain[t, i, j], dt, Q)
            
        # prepare visualization data        
        f = np.zeros((dom_shape[0], dom_shape[1], 3))
        f_na = np.zeros((dom_shape[0], dom_shape[1], 3))
        for p in np.ndindex(dom_shape):
            f[p[0], p[1], :] = models[p].get_state()[:3]
            f_na[p[0], p[1], :] = models_na[p].get_state()[:3]
            mV[pos] = models[p].get_state_covar()[1,1]
	    mid[p] = models[p].get_model_ids()[1]

        # run Kriging on each observed fuel type
        Kf = []
        Vf = []
        fn = []
        for obs_data, fuel_ndx in [ (obs_data_fm10, 1) ]:

            if model_time in obs_data:

                # fit the current estimation of the moisture field to the data 
                base_field = f[:,:,fuel_ndx]
                mfm.fit_to_data(base_field, obs_data[model_time])
                
                # find differences (residuals) between observed measurements and nearest grid points
                # use this to update observation residual standard deviation 
                obs_vals = np.array([o.get_value() for o in obs_data[model_time]])
                ngp_vals = np.array([base_field[o.get_nearest_grid_point()] for o in obs_data[model_time]])
                mod_vals = np.array([f[:,:,fuel_ndx][o.get_nearest_grid_point()] for o in obs_data[model_time]])
                mod_na_vals = np.array([f_na[:,:,fuel_ndx][o.get_nearest_grid_point()] for o in obs_data[model_time]])
                obs_re.update_with(obs_vals - ngp_vals)
                diagnostics().push("kriging_obs_res_var", (t, np.mean(obs_re.get_variance())))
            
                # predict the moisture field using observed fuel type
                predicted_field = mfm.predict_field(base_field)

                # update the model residual estimator and get current best estimate of variance
                mod_re.update_with(f[:,:,fuel_ndx] - predicted_field)
                mresV = mod_re.get_variance()

                # krige data to observations
                Kf_fn, Vf_fn = simple_kriging_data_to_model(obs_data[model_time], obs_re.get_variance() ** 0.5,
                                                            predicted_field, wrf_data, mresV ** 0.5, t)

                krig_vals = np.array([Kf_fn[o.get_nearest_grid_point()] for o in obs_data[model_time]])                
                diagnostics().push("assim_data", (t, fuel_ndx, obs_vals, ngp_vals, krig_vals, mod_vals, mod_na_vals))

                
                # append to storage for kriged fields in this time instant
                Kf.append(Kf_fn)
                Vf.append(Vf_fn)
                fn.append(fuel_ndx)

        # if there were any observations, run the kalman update step
        if len(fn) > 0:
            Nobs = len(fn)
            Kg = np.zeros((dom_shape[0], dom_shape[1], Nobs))
            # run the kalman update in each model independently
            # gather the standard deviations of the moisture fuel after the Kalman update
            for pos in np.ndindex(dom_shape):
                O = np.zeros((Nobs,))
                V = np.zeros((Nobs, Nobs))
                
                # construct observations for this position
                for i in range(Nobs):
                    O[i] = Kf[i][pos]
                    V[i,i] = Vf[i][pos]
                
                # execute the Kalman update 
                Kg[pos[0], pos[1], :] = models[pos].kalman_update(O, V, fn)
                
            # push the mean Kalman gain
            diagnostics().push("assim_K", (t, np.mean(Kg[:,:,0])))
            

        # prepare visualization data        
        f = np.zeros((dom_shape[0], dom_shape[1], 3))
        for p in np.ndindex(dom_shape):
            f[p[0], p[1], :] = models[p].get_state()[:3]
            
        plt.clf()
        plt.subplot(3,3,1)
        render_spatial_field_fast(m, lon, lat, f[:,:,0], '1-hr fuel')
        plt.clim([0.0, maxE])
        plt.colorbar()
        plt.subplot(3,3,2)
        render_spatial_field_fast(m, lon, lat, f[:,:,1], '10-hr fuel')
        plt.clim([0.0, maxE])        
        plt.colorbar()
        plt.subplot(3,3,3)
#        render_spatial_field_fast(m, lon, lat, f[:,:,2], '100-hr fuel')
	render_spatial_field_fast(m, lon, lat, mid, 'Model ids')
#        plt.clim([0.0, maxE])        
	plt.clim([0.0, 5.0])
        plt.colorbar()
        plt.subplot(3,3,4)
        render_spatial_field_fast(m, lon, lat, predicted_field, 'Mean field Fit')
        plt.clim([0.0, maxE])        
        plt.colorbar()
        plt.subplot(3,3,5)
        render_spatial_field_fast(m, lon, lat, Kg[:,:,0], 'Kalman gain for 10-hr fuel')       
        plt.clim([0.0, 1.0])        
        plt.colorbar()
        plt.subplot(3,3,6)
#        render_spatial_field_fast(m, lon, lat, mV, '10-hr fuel variance')
#        plt.clim([0.0, np.max(mV)]) 
        render_spatial_field_fast(m, lon, lat, f_na[:,:,1], '10hr fuel - no assim')
        plt.clim([0.0, maxE])
        plt.colorbar()
        plt.subplot(3,3,7)
        render_spatial_field_fast(m, lon, lat, Kf_fn, 'Kriged observations')
        plt.clim([0.0, maxE])
        plt.colorbar()
        plt.subplot(3,3,8)
        render_spatial_field_fast(m, lon, lat, Vf_fn, 'Kriging variance')
        plt.clim([0.0, np.max(Vf_fn)])
        plt.colorbar()
        plt.subplot(3,3,9)
        render_spatial_field_fast(m, lon, lat, mresV, 'Model res. variance')
        plt.clim([0.0, np.max(mresV)])
        plt.colorbar()
        
        plt.savefig(os.path.join(cfg['output_dir'], 'moisture_model_t%03d.png' % t))

        # push new diagnostic outputs
        diagnostics().push("assim_K", (t, np.mean(Kg[:,:,0])))
        diagnostics().push("assim_mV", (t, np.mean(mV)))
        diagnostics().push("assim_mresV", (t, np.mean(mresV)))
        diagnostics().push("kriging_variance", (t, np.mean(Vf_fn)))

        
    # store the gamma coefficients
    with open(os.path.join(cfg['output_dir'], 'gamma.txt'), 'w') as f:
        f.write(str(diagnostics().pull('mfm_gamma')))
        
    # make a plot of gammas
    plt.figure()
    plt.subplot(211)
    plt.plot(diagnostics().pull('mfm_gamma'))
    plt.title('Mean field model - gamma')
    plt.subplot(212)
    plt.plot(diagnostics().pull('skdm_cov_cond'))
    plt.title('Condition number of covariance matrix')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_gamma.png'))
    
    # make a plot for each substation
    plt.figure()
    for i in range(len(stations)):
        plt.clf()
        # get data for the i-th station
        t_i = [ o[0] for o in diagnostics().pull("assim_data") ]
        obs_i = [ o[2][i] for o in diagnostics().pull("assim_data") ]
        bf_i = [ o[3][i] for o in diagnostics().pull("assim_data") ]
        krig_i = [ o[4][i] for o in diagnostics().pull("assim_data") ]
        mod_i = [ o[5][i] for o in diagnostics().pull("assim_data") ]
        mod_na_i = [ o[6][i] for o in diagnostics().pull("assim_data") ]
        mx = max(max(obs_i), max(mod_i), max(krig_i), max(mod_i))
        plt.plot(t_i, obs_i, 'ro')
        plt.plot(t_i, bf_i, 'gx-')
        plt.plot(t_i, krig_i, 'bo-')
        plt.plot(t_i, mod_i, 'kx-', linewidth = 1.5)
        plt.plot(t_i, mod_na_i, 'mx-')
        plt.ylim([0.0, 1.1 * mx])
        plt.legend(['Obs.', 'Base field', 'Kriged', 'Model', 'NoAssim'])
        plt.title('Station observations fit to model and kriging field')
        plt.savefig(os.path.join(cfg['output_dir'], 'station%02d.png' % (i+1)))
        
    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("assim_K")],
             [d[1] for d in diagnostics().pull("assim_K")], 'ro-')
    plt.title('Average Kalman gain')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_kalman_gain.png'))
    
    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("assim_mV")],
             [d[1] for d in diagnostics().pull("assim_mV")], 'ro-')
    plt.title('Average model variance')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_fm10_model_variance.png'))

    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("assim_mresV")],
             [d[1] for d in diagnostics().pull("assim_mresV")], 'ro-')
    plt.title('Average fm10 residual variance')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_fm10_model_residual_variance.png'))

    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("kriging_variance")],
             [d[1] for d in diagnostics().pull("kriging_variance")], 'ro-')
    plt.title('Kriging field variance')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_kriging_variance.png'))

    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("kriging_obs_res_var")],
             [d[1] for d in diagnostics().pull("kriging_obs_res_var")], 'ro-')
    plt.title('Observation residual variance')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_observation_residual_variance.png'))
    
    diagnostics().dump_store(os.path.join(cfg['output_dir'], 'diagnostics.bin'))
    
    # as a last step encode all the frames as video
    os.system("cd %s; avconv -qscale 1 -r 20 -b 9600 -i moisture_model_t%%03d.png video.mp4" % cfg['output_dir'])


if __name__ == '__main__':
#    profile.run('run_module(); print', 'spatial_model.stats')
#    
#    stats = pstats.Stats('spatial_model.stats')
#    stats.strip_dirs()
#    stats.sort_stats('cumulative')
#    stats.print_stats()

    run_module()
