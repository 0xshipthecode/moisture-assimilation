# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:36 2012

@author: martin
"""

from spatial_model_utilities import render_spatial_field_fast, great_circle_distance
from time_series_utilities import build_observation_data

from kriging_methods import universal_kriging_data_to_model, trend_surface_model_kriging

from wrf_model_data import WRFModelData
from cell_model_opt import CellMoistureModel
from mean_field_model import MeanFieldModel
from observation_stations import MesoWestStation
from diagnostics import init_diagnostics, diagnostics
from online_variance_estimator import OnlineVarianceEstimator

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pytz
import cPickle
import string


def plot_model_snapshot(cfg, tm, t, fuel_ndx, obs, krig, mod, mod_na):
    """
    Plot the model values at observations points, observations, kriging results etc.
    """
    plt.figure()

    plt.plot(mod_na, 'go', markersize = 5)
    plt.plot(mod, 'bo', markersize = 5)
    leg = [ 'Model+Assim', 'Model',]
    mx = max(max(mod), max(mod_na), 0.5)

    for (v,l,c) in [ (obs, 'Obs.', 'ro'), (krig, 'Kriged', 'mx') ]:
        # adjust plot depending on whether observations are available
        if v is not None:
            mx = max(max(v), mx)
            leg.append(l)
            plt.plot(v, c, markersize = 5)

    plt.ylim([0.0, 1.1 * mx])
    plt.legend(leg)
    plt.title('Model behavior for time %s' % str(tm[t]))
    plt.savefig(os.path.join(cfg['output_dir'], 'model_snapshot_f%d_t%03d.png' % (fuel_ndx, t)))


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

    # Error covariance matrix condition number in kriging
    diagnostics().configure_tag("skdm_cov_cond", False, True, True)

    # Assimilation parameters
    diagnostics().configure_tag("assim_K0", False, True, True)
    diagnostics().configure_tag("assim_K1", True, True, True)
    diagnostics().configure_tag("assim_data", False, False, True)

    diagnostics().configure_tag("kriging_variance", False, False, True)

    diagnostics().configure_tag("obs_residual_var", True, True, True)
    diagnostics().configure_tag("fm10_model_residual_var", True, True, True)
    diagnostics().configure_tag("fm10_model_var", True, True, True)
    diagnostics().configure_tag("fm10_kriging_var", True, True, True)

    ### Load and preprocess WRF model data

    # load WRF data
    wrf_data = WRFModelData(cfg['input_file'], tz_name = 'US/Mountain')
    
    # read in spatial and temporal extent of WRF variables
    lat, lon = wrf_data.get_lats(), wrf_data.get_lons()
    tm = wrf_data.get_gmt_times()
    Nt = cfg['Nt'] if cfg['Nt'] is not None else len(tm)
    dom_shape = lat.shape

    # interpolate rain in the same way as equilibria are
    rain = wrf_data['RAINNC']
    rain = 0.5 * (rain[:-1,:,:] + rain[1:,:,:])

    # moisture equilibria are now computed from averaged Q,P,T at beginning and end of period
    Ed, Ew = wrf_data.get_moisture_equilibria()

    ### Load observation data from the stations

    # load station data from files
    with open(os.path.join(cfg['station_data_dir'], cfg['station_list_file']), 'r') as f:
        si_list = f.read().split('\n')

    si_list = filter(lambda x: len(x) > 0, map(string.strip, si_list))

    # for each station id, load the station
    stations = []
    for sinfo in si_list:
        code = sinfo.split(',')[0]
        mws = MesoWestStation(sinfo, wrf_data)
        for suffix in [ '_1', '_2', '_3', '_4', '_5', '_6', '_7' ]:
            mws.load_station_data(os.path.join(cfg['station_data_dir'], '%s%s.xls' % (code, suffix)))
        stations.append(mws)

    print('Loaded %d stations.' % len(stations))
    
    # check stations for nans
    stations = filter(MesoWestStation.data_ok, stations)
    print('Have %d stations with complete data.' % len(stations))

    # set the measurement variance of the stations
    for s in stations:
        s.set_measurement_variance('fm10', cfg['fm10_meas_var'])

    # build the observation data
    obs_data_fm10 = build_observation_data(stations, 'fm10', wrf_data, tm)

    ### Initialize model and visualization

    # find maximum moisture overall to set up visualization
    maxE = 0.5
    
    # construct initial conditions
    E = 0.5 * (Ed[0,:,:] + Ew[0,:,:])
    
    # set up parameters
    Q = np.eye(9) * cfg['Q']
    P0 = np.eye(9) * cfg['P0']
    dt = (tm[1] - tm[0]).seconds
    print("INFO: Computed timestep from WRF is is %g seconds." % dt)
    K = np.zeros_like(E)
    V = np.zeros_like(E)
    mV = np.zeros_like(E)
    predicted_field = np.zeros_like(E)
    mresV = np.zeros_like(E)
    Kf_fn = np.zeros_like(E)
    Vf_fn = np.zeros_like(E)
    mid = np.zeros_like(E)
    Kg = np.zeros((dom_shape[0], dom_shape[1], 9))
    cV12 = np.zeros_like(E)
    
    # moisture state and observation residual variance estimators
    mod_re = OnlineVarianceEstimator(np.zeros_like(E), np.ones_like(E) * 0.05, 1)
    obs_re = OnlineVarianceEstimator(np.zeros((len(stations),)), np.ones(len(stations),) * 0.05, 1)
    
    # initialize the mean field model (default fit is 1.0 of equilibrium before new information comes in)
    mfm = MeanFieldModel(cfg['lock_gamma'])

    # construct model grid using standard fuel parameters
    Tk = np.array([1.0, 10.0, 100.0]) * 3600
    models = np.zeros(dom_shape, dtype = np.object)
    models_na = np.zeros_like(models)
    for p in np.ndindex(dom_shape): 
        models[p] = CellMoistureModel((lat[p], lon[p]), 3, E[p], Tk, P0 = P0)
        models_na[p] = CellMoistureModel((lat[p], lon[p]), 3, E[p], Tk, P0 = P0)

    m = None
    plt.figure(figsize = (12, 8))
    
    ###  Run model for each WRF timestep and assimilate data when available
    for t in range(1, Nt):
        model_time = tm[t]
        print("INFO: time: %s, step: %d" % (str(model_time), t))

        # run the model update
        for p in np.ndindex(dom_shape):
            i, j = p
            models[p].advance_model(Ed[t-1, i, j], Ew[t-1, i, j], rain[t-1, i, j], dt, Q)
            models_na[p].advance_model(Ed[t-1, i, j], Ew[t-1, i, j], rain[t-1, i, j], dt, Q)
            
        # prepare visualization data
        f = np.zeros((dom_shape[0], dom_shape[1], 3))
        f_na = np.zeros((dom_shape[0], dom_shape[1], 3))
        for p in np.ndindex(dom_shape):
            f[p[0], p[1], :] = models[p].get_state()[:3]
            f_na[p[0], p[1], :] = models_na[p].get_state()[:3]
            P = models[p].get_state_covar()
            cV12[p] = P[0,1]
            mV[p] = P[1,1]
            mid[p] = models[p].get_model_ids()[1]

        diagnostics().push("fm10_model_var", np.mean(mV))

        # run Kriging on each observed fuel type
        Kf = []
        Vf = []
        fn = []
        for obs_data, fuel_ndx in [ (obs_data_fm10, 1) ]:


            # run the kriging subsystem and the Kalman update only if we have observations
            if model_time in obs_data:

                # fit the current estimation of the moisture field to the data 
                base_field = f[:,:,fuel_ndx]
                mfm.fit_to_data(base_field, obs_data[model_time])
                
                # find differences (residuals) between observed measurements and nearest grid points
                # use this to update observation residual standard deviation 
                obs_vals = np.array([o.get_value() for o in obs_data[model_time]])
                mod_vals = np.array([base_field[o.get_nearest_grid_point()] for o in obs_data[model_time]])
                mod_na_vals = np.array([f_na[:,:,fuel_ndx][o.get_nearest_grid_point()] for o in obs_data[model_time]])
                obs_re.update_with(obs_vals - mod_vals)
                diagnostics().push("obs_residual_var", (t, np.mean(obs_re.get_variance())))
            
                # predict the moisture field using observed fuel type
                predicted_field = mfm.predict_field(base_field)

                # update the model residual estimator and get current best estimate of variance
                mod_re.update_with(f[:,:,fuel_ndx] - predicted_field)
                mresV = mod_re.get_variance()
                diagnostics().push("fm10_model_residual_var", (t, np.mean(mresV)))

                # krige observations to grid points
                Kf_fn, Vf_fn = trend_surface_model_kriging(obs_data[model_time], wrf_data, predicted_field)

                krig_vals = np.array([Kf_fn[o.get_nearest_grid_point()] for o in obs_data[model_time]])                
                diagnostics().push("assim_data", (t, fuel_ndx, obs_vals, krig_vals, mod_vals, mod_na_vals))
                plot_model_snapshot(cfg, tm, t, fuel_ndx, obs_vals, krig_vals, mod_vals, mod_na_vals)

                diagnostics().push("fm10_kriging_var", np.mean(Vf_fn))

                # append to storage for kriged fields in this time instant
                Kf.append(Kf_fn)
                Vf.append(Vf_fn)
                fn.append(fuel_ndx)

#            else:
#
                # plot the behavior of the model outside of observations
#                plot_model_snapshot(cfg, tm, t, fuel_ndx, None, None, mod_vals, mod_na_vals)
                

        # if there were any observations, run the kalman update step
        if len(fn) > 0:
            Nobs = len(fn)
            # run the kalman update in each model independently
            # gather the standard deviations of the moisture fuel after the Kalman update
            for p in np.ndindex(dom_shape):
                O = np.zeros((Nobs,))
                V = np.zeros((Nobs, Nobs))
                
                # construct observations for this position
                for i in range(Nobs):
                    O[i] = Kf[i][p]
                    V[i,i] = Vf[i][p]
                
                # execute the Kalman update
                Kp = models[p].kalman_update(O, V, fn)
                Kg[p[0], p[1], :] = Kp[:, 0]

            # push new diagnostic outputs
            diagnostics().push("assim_K0", (t, np.mean(Kg[:,:,0])))
            diagnostics().push("assim_K1", (t, np.mean(Kg[:,:,1])))
            diagnostics().push("kriging_variance", (t, np.mean(Vf_fn)))

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
        render_spatial_field_fast(m, lon, lat, f_na[:,:,1], '10hr fuel - no assim')
        plt.clim([0.0, maxE])
        plt.colorbar()
        plt.subplot(3,3,4)
        render_spatial_field_fast(m, lon, lat, Kg[:,:,0], 'Kalman gain for 1-hr fuel')  
        plt.clim([0.0, 3.0])        
        plt.colorbar()
        plt.subplot(3,3,5)
        render_spatial_field_fast(m, lon, lat, Kg[:,:,1], 'Kalman gain for 10-hr fuel')       
        plt.clim([0.0, 1.0])        
        plt.colorbar()
        plt.subplot(3,3,6)
	render_spatial_field_fast(m, lon, lat, Kf_fn, 'Kriging field')
	plt.clim([0.0, maxE])
	plt.colorbar()
        plt.subplot(3,3,7)
        render_spatial_field_fast(m, lon, lat, mid, 'Model ids')
        plt.clim([0.0, 5.0])
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
    D = diagnostics().pull("assim_data")
    for i in range(len(stations)):
        plt.clf()
        # get data for the i-th station
        t_i = [ o[0] for o in D]
        obs_i = [ o[2][i] for o in D]
        krig_i = [ o[3][i] for o in D]
        mod_i = [ o[4][i] for o in D]
        mod_na_i = [ o[5][i] for o in D]
        mx = max(max(obs_i), max(mod_i), max(krig_i), max(mod_i))
        plt.plot(t_i, obs_i, 'ro')
        plt.plot(t_i, krig_i, 'bo-')
        plt.plot(t_i, mod_i, 'kx-', linewidth = 1.5)
        plt.plot(t_i, mod_na_i, 'mx-')
        plt.ylim([0.0, 1.1 * mx])
        plt.legend(['Obs.', 'Kriged', 'Model', 'NoAssim'])
        plt.title('Station observations fit to model and kriging field')
        plt.savefig(os.path.join(cfg['output_dir'], 'station%02d.png' % (i+1)))
        
    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("assim_K1")],
             [d[1] for d in diagnostics().pull("assim_K1")], 'ro-')
    plt.title('Average Kalman gain')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_kalman_gain_10hr.png'))
    
    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("assim_K0")],
             [d[1] for d in diagnostics().pull("assim_K0")], 'ro-')
    plt.title('Average Kalman gain')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_kalman_gain_1hr.png'))

    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("fm10_model_var")],
             [d[1] for d in diagnostics().pull("fm10_model_var")], 'ro-')
    plt.title('Average fm10 model variance')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_fm10_model_variance.png'))

    plt.figure()
    plt.plot([d[0] for d in diagnostics().pull("fm10_model_residual_var")],
             [d[1] for d in diagnostics().pull("fm10_model_residual_var")], 'ro-')
    plt.title('Average fm10 model residual variance')
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
    
    plt.figure()
    plt.plot(diagnostics().pull("mfm_mape"), 'ro-', linewidth = 2)
    plt.title('Mean absolute prediction error of station data')
    plt.savefig(os.path.join(cfg['output_dir'], 'plot_station_mape.png'))

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
