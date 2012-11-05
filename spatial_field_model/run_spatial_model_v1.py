# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:36 2012

@author: martin
"""

from spatial_model_utilities import render_spatial_field, equilibrium_moisture, load_stations_from_files, \
                                    match_stations_to_gridpoints, match_sample_times, great_circle_distance

from wrf_model_data import WRFModelData
from cell_model import CellMoistureModel

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os

#import profile
#import pstats

import multiprocessing


station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
#                  "Mt Laguna_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Moisture",
#                  "Palomar_Fuel_Moisture",
                  "Alpine_Moisture",
#                  "Valley_Center_Moisture",
                  "Ranchita_Moisture",
                  "Camp_Elliot_Moisture",
#                  "Pine Hills_Moisture" 
                ]

                  
station_data_dir = "../real_data/witch_creek/"


def construct_correlation_matrix(gridndx, mlons, mlats):
    """
    Construct a distance-based correlation matrix between residuals at given longitudes
    and lattitudes.
    """
    N = len(gridndx)
    D = np.zeros((N,N))
    
    # compute distances in km between locations
    for (i,j), i1 in zip(gridndx, range(N)):
        lon1, lat1 = mlons[i,j], mlats[i,j]
        for (k,l), i2 in zip(gridndx, range(N)):
            if i1 != i2:
                lon2, lat2 = mlons[k,l], mlats[k,l]
                D[i1,i2] = great_circle_distance(lon1, lat1, lon2, lat2)
            
    # estimate correlation coeff
    return np.maximum(np.eye(N), 0.8565 - 0.0063 * D)


def simple_kriging_data_to_model(obs_data, W, mS, t):
    """
    Simple kriging of data points to model points.  The kriging results in
    the matrix K, which contains mean of the kriged observations and the
    matrix V, which contains the kriging variance. 
    
        synopsis: K, V = simple_kriging_data_to_model(obs_data, W, t)
        
    """
    mlons, mlats = W.get_lons(), W.get_lats()
    P, Q, T = W['PSFC'][t,:,:], W['Q2'][t,:,:], W['T2'][t,:,:] 
    K = np.zeros_like(mlons)
    V = np.zeros_like(mlons)
    Nobs = len(obs_data)
    fm_obs = np.zeros((Nobs,))
    fm_stds = np.zeros((Nobs,))
    station_lonlat = []
        
    # accumulate the indices of the nearest grid points
    ndx = 0
    gridndx = []
    for mr in obs_data.values():
        fm_obsi, grid_pos, lonlat, fmres_std = mr['fm_obs'], mr['nearest_grid_point'], mr['lonlat'], mr['fm_std']
        gridndx.append(grid_pos)
        fm_obs[ndx] = fm_obsi
        fm_stds[ndx] = fmres_std
        station_lonlat.append(lonlat)
        ndx += 1
        
    # compute nominal state for grid points
    Ed, Ew = equilibrium_moisture(P, Q, T)
    mu_mod = 0.899022 * 0.5 * (Ed + Ew)

    # compute nominal state for station data
    mu_obs = np.zeros((Nobs,))
    for g, i in zip(gridndx, range(Nobs)):
        mu_obs[i] = mu_mod[g]

    # compute observation residuals (using model fit from examine_station_data)
    res_obs = fm_obs - mu_obs
    
    # construct the covariance matrix and invert it
    C = construct_correlation_matrix(gridndx, mlons, mlats)
    oS = np.diag(fm_stds)
    Sigma = np.dot(np.dot(oS, C), oS)
    SigInv = np.linalg.inv(Sigma)
    
    # run the kriging estimator for each model grid point
    K = np.zeros_like(mlats)
    cov = np.zeros_like(mu_obs)
    for p in np.ndindex(K.shape):
        # compute the covariance array anew for each grid point
        for k in range(Nobs):
            lon, lat = station_lonlat[k]
            cc = max(0.8565 - 0.0063 * great_circle_distance(mlons[p], mlats[p], lon, lat), 0.0)
            cov[k] = mS[p] * cc * fm_stds[k]
        csi = np.dot(cov, SigInv)
        K[p] = mu_mod[p] + np.dot(csi, res_obs) 
        V[p] = mS[p]**2 - np.dot(csi, cov)
    
    return K, V


def build_observation_data(stations, W):
    """
    Repackage the matched time series into a time-indexed structure which gives details on the observed data and active observation stations.
    
        synopsis: obs_data = build_observation_data(stations, fm_ts)
        
    """
    Ns = len(stations)
    
    # iterate over time instants and accumulate them into observation packets
    obs_data = {}
    for sname, s in stations.iteritems():
        i, j = s['nearest_grid_point']
        mtm, ndx1, _ = match_sample_times(W.get_times(), sorted(s['fuel_moisture'].keys()))
        Ed, Ew = equilibrium_moisture(W['PSFC'][ndx1, i, j], W['Q2'][ndx1, i, j], W['T2'][ndx1, i, j])
        fm_equi = 0.899022 * 0.5 * (Ed + Ew)
        fm_st = [ s['fuel_moisture'][t] for t in mtm ]
        fm_std = np.std(fm_st - fm_equi) # estimate the standard deviation of the residuals
        for tm, obs in  zip(mtm, fm_st):
            obs_i = obs_data[tm] if tm in obs_data else {}
            obs_s = { 'fm_obs' : obs,
                      'fm_std' : fm_std,
                      'nearest_grid_point' : (i,j),
                      'lonlat' : (s['lon'], s['lat']),
                      'name' : sname }
            obs_i[sname] = obs_s
            obs_data[tm] = obs_i
            
    return obs_data
    
    
def advance_cell_model(x):
    model, pos, Ed, Ew, r, dt, Qij = x
    model.advance_model(Ed, Ew, r, dt, Qij)
    return (pos, model)

def run_module():
        
    W = WRFModelData('../real_data/witch_creek/realfire03_d04_20071022.nc')
    
    # read in vars
    lat, lon = W.get_lats(), W.get_lons()
    tm = W.get_times()
    rain = W['RAINNC']
    Q2 = W['Q2']
    T2 = W['T2']
    P = W['PSFC']
    
    # obtain sizes
    Nt = rain.shape[0]
    dom_shape = lat.shape
    locs = np.prod(dom_shape)
    
    # load station data
    stations = load_stations_from_files(station_data_dir, station_list, 'US/Pacific')
    match_stations_to_gridpoints(stations, lon, lat)
    
    # manipulate observation data into a time indexed structure
    obs_data = build_observation_data(stations, W) 
    
    # construct initial vector
    Ed, Ew = equilibrium_moisture(P[1,:,:], Q2[1,:,:], T2[1,:,:])
    E = 0.5 * (Ed + Ew)
    
    # set up parameters
    mS = np.zeros_like(E)
    Qij = np.eye(9) * 0.001
    dt = 10.0 * 60
    K = np.zeros_like(E)
    V = np.zeros_like(E)
    Kg = np.zeros_like(E)

    # construct model grid using standard fuel parameters
    Tk = np.array([1.0, 10.0, 100.0]) * 3600
    models = np.zeros(dom_shape, dtype = np.object)
    for pos in np.ndindex(dom_shape): 
        models[pos] = CellMoistureModel((lat[pos], lon[pos]), 3, E[pos], Tk, P0 = Qij)
    
    
    # construct a basemap representation of the area
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    plt.figure(figsize = (10, 6))
    
    # run model
    for t in range(1, Nt):
        model_time = W.get_times()[t]
        print("Time: %s, step: %d" % (str(model_time), t))

        # pre-compute equilibrium moisture to save a lot of time
        Ed, Ew = equilibrium_moisture(P[t,:,:], Q2[t,:,:], T2[t,:,:])
        
        # run the model update (use multiprocessing)
        for pos in np.ndindex(dom_shape):
            i, j = pos
            models[pos].advance_model(Ed[i, j], Ew[i, j], rain[t, i, j], dt, Qij)
            
        # gather the standard deviations of the moisture fuel of the observed type from the model
        for pos in np.ndindex(dom_shape):
            mS[pos] = models[pos].P[1,1]**0.5

        # if we have an observation somewhere in time
        doing_kalman_update = False
        if model_time in obs_data:
            doing_kalman_update = True
            
            # krige data to observations
            K, V = simple_kriging_data_to_model(obs_data[model_time], W, mS, t)

            # run the kalman update in each model
            # gather the standard deviations of the moisture fuel after the Kalman update
            for pos in np.ndindex(dom_shape):
                Kg[pos] = models[pos].kalman_update(K[pos], V[pos], 1)
                mS[pos] = models[pos].P[1,1]

        # prepare visualization data        
        f = np.zeros((dom_shape[0], dom_shape[1], 3))
        for p in np.ndindex(dom_shape):
            f[p[0], p[1], :] = models[p].get_state()[:3]
        E = 0.5 * (Ed + Ew)
            
        plt.clf()
        plt.subplot(3,3,1)
        render_spatial_field(m, lon, lat, f[:,:,0], 'Fast fuel')
        plt.clim([0.0, 0.2])
        plt.colorbar()
        plt.subplot(3,3,2)
        render_spatial_field(m, lon, lat, f[:,:,1], 'Mid fuel')
        plt.clim([0.0, 0.2])        
        plt.colorbar()
        plt.subplot(3,3,3)
        render_spatial_field(m, lon, lat, f[:,:,2], 'Slow fuel')
        plt.clim([0.0, 0.2])        
        plt.colorbar()
        plt.subplot(3,3,4)
        render_spatial_field(m, lon, lat, E, 'Equilibrium')
        plt.clim([0.0, 0.5])        
        plt.colorbar()
        plt.subplot(3,3,5)
#        render_spatial_field(m, lon, lat, rain[t,:,:], 'Rain')
        render_spatial_field(m, lon, lat, Kg, 'Kalman gain')       
        plt.clim([0.0, 1.0])        
        plt.colorbar()
        plt.subplot(3,3,6)
        render_spatial_field(m, lon, lat, T2[t,:,:] - 273.15, 'Temperature')
        plt.colorbar()
        plt.subplot(3,3,7)
        render_spatial_field(m, lon, lat, K, 'Kriged observations')
        plt.colorbar()
        plt.subplot(3,3,8)
        render_spatial_field(m, lon, lat, V, 'Kriging variance')
        plt.clim([0.0, np.max(V)])
        plt.colorbar()
        plt.subplot(3,3,9)
        render_spatial_field(m, lon, lat, mS, 'Model variance')
        plt.clim([0.0, np.max(mS)])
        plt.colorbar()
        
        plt.savefig('model_outputs/moisture_model_t%03d.png' % t) 
        
                
    
if __name__ == '__main__':
#    profile.run('run_module(); print', 'spatial_model.stats')
#    
#    stats = pstats.Stats('spatial_model.stats')
#    stats.strip_dirs()
#    stats.sort_stats('cumulative')
#    stats.print_stats()

    run_module()
