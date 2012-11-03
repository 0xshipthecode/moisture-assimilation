# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:36 2012

@author: martin
"""

from spatial_model_utilities import load_station_data, render_spatial_field, \
                                    equilibrium_moisture, load_stations_from_files, \
                                    match_stations_to_gridpoints, match_sample_times, \
                                    great_circle_distance

from wrf_model_data import WRFModelData

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os


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



def moisture_model_simple(T, Tk, Q, P, m_ext, r, dt):
    """
    This model captures the moisture dynamics at each grid point independently.
    
    synopsis: m_ext = moisture_model_ext(T, Tk, Q, P, m_ext, f_info, r, dt)
    
    T - temperature [Kelvins]
    Tk - time delay nominal values for fuel classes
    Q - water vapor ratio [dimensionless]
    P - pressure [Pascals]
    m_ext - current moisture values and assimilated parameters
    r - rain intensity for time unit [mm/h]
    dt - integration step [s]
    """
    k = Tk.shape[0]                 # number of fuel classes
    r0 = 0.05                       # threshold rainfall [mm/h]
    rk = 8                          # saturation rain intensity [mm/h]
    Trk = 14 * 3600                 # time constant for wetting model [s]
    S = 2.5                         # saturation intensity [dimensionless]
    
    # first, we break the state vector into components
    m = m_ext[:k]
    dlt_Tk = m_ext[k:k+k]
    dlt_E = m_ext[2*k]
    dlt_S = m_ext[2*k+1]
    dlt_Trk = m_ext[2*k+2]
    
    # compute equilibrium moisture using model
    Ed, Ew = equilibrium_moisture(P, Q, T)
    
    # add assimilated difference, which is shared across spatial locations
    Ed = Ed + dlt_E
    Ew = Ew + dlt_E
    
    # where rainfall is above threshold (spatially different), apply
    # saturation model, equi and rlag are specific to fuel type and
    # location
    equi = m.copy()         # copy over current equilibrium levels
    rlag = np.zeros((k,))
    model_ids = np.zeros((k,))
    
    # equilibrium is equal to the saturation level (assimilated)
    if r > r0:
        equi[:] = S + dlt_S
        model_ids[:] = 3
    
        # rlag is modified by the rainfall intensity (assimilated)
        rlag = 1.0 / (Trk + dlt_Trk) * (1 - np.exp(- (r - r0) / rk))
    else:

        # equilibrium is selected according to current moisture level
        model_ids[:] = 4
        model_ids[equi > Ed] = 1
        equi[equi > Ed] = Ed
        model_ids[equi < Ew] = 2
        equi[equi < Ew] = Ew

        # the inverted time lag is constant according to fuel category
        rlag = 1.0 / (Tk + dlt_Tk)
    
    # select appropriate integration method according to change for each fuel 
    # and location
    change = dt * rlag
    m_new = np.zeros_like(m_ext)
    m_new[:k] = np.where(change < 0.01,
                          m + (equi - m) * (1 - np.exp(-change)),
                          m + (equi - m) * change * (1 - 0.5 * change))
    m_new[k:] = m_ext[k:]
    return m_new
    

def construct_correlation_matrix(gridndx, mlons, mlats):
    
    N = len(gridndx)
    D = np.zeros((N,N))
    
    # compute distances in km between locations
    for (i,j), i1 in zip(gridndx, range(N)):
        lon1, lat1 = mlons[i,j], mlats[i,j]
        for (k,l), i2 in zip(gridndx, range(N)):
            lon2, lat2 = mlons[i,j], mlats[i,j]
            C[i1,i2] = great_circle_distance(lon1, lat1, lon2, lat2)
            
    # estimate correlation coeff
    C = np.max(np.zeros_like(D), 0.8565 - 0.0063 * D)
    return C


def simple_kriging_data_to_model(data, W):
    """
    Simple kriging of data points to model points.  The kriging results in
    the matrix K, which contains mean of the kriged observations and the
    matrix V, which contains the kriging standard deviations. 
    
        synopsis: K, S = simple_kriging_data_to_model(obs_data, mlons, mlats, mS, P, Q, T)
        
    """
    K = np.zeros_like(mlons)
    S = np.zeros_like(mlons)
    Nobs = len(obs_data)
    fm_obs = np.zeros((Nobs,))
    fm_stds = np.zeros((Nobs,))
    station_lonlat = []
    
    # nominal values (mean) for model
    mu_mod = 0.5 * (Ed + Ew)
    
    # accumulate the indices of the nearest grid points
    ndx = 0
    gridndx = []
    for mr in obs_data.values():
        fm_obsi, grid_pos, lonlat, fmres_std = mr['fm_obs'], mr['nearest_grid_point'], lonlat = mr['lonlat'], mr['fm_std']
        gridndx.append(grid_pos)
        fm_obs[ndx] = fm_obsi
        fm_stds[ndx] = fmres_std
        station_lonlat.append(lonlat)
        ndx += 1

    # compute nominal state for station data
    Ed, Ew = equilibrium_moisture(P[gridndx], Q[gridndx], T[gridndx])
    mu_obs = 0.5 * (Ed + Ew)
    
    # compute nominal state for grid points
    Ed, Ew = equilibrium_moisture(P, Q, T)
    mu_mod = 0.5 * (Ed + Ew)

    # compute observation residuals
    res_obs = Z_obs - mu_obs
    
    # construct the covariance matrix and invert it
    sigmas = np.diag(V[gridndx]) ** 0.5
    C = construct_correlation_matrix(gridndx, mlons, mlats)
    Sigma = np.dot(np.dot(vars, C), vars)
    SigInv = np.linalg.inv(Sigma)
    
    # run the kriging estimator for each model grid point
    K = np.zeros_like(mlat)
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            cov = np.zeros_like(mu_obs)
            for k in range(Nobs):
                lon, lat = station_lonlat[k]
                cc = max(0.8565 - 0.0063 * great_circle_distance(mlons[i,j], mlats[i,j], lon, lat), 0.0)
                cov[k] = mS[i,j] * cc * obs_vars[k]
                csi = np.dot(cov, SigInv)
                K[i,j] = np.dot(csi, res_obs) + mu_mod[i,j]
                S[i,j] = (mS[i,j]**2 - np.dot(csi, cov)) ** 0.5
    
    
    return K, S


def build_observation_data(stations, fm_ts):
    """
    Repackage the matched time series into a time-indexed structure which gives details on the observed data and active observation stations.
    
        synopsis: obs_data = build_observation_data(stations, fm_ts)
        
    """
    Ns = len(stations)
    vec_ndx = np.zeros((Ns,))
    len_ndx = np.array([len(s['t']) for s in fm_ts.values()])
    # initialize start time of observations
    tm = min([ s['t'][0] for s in fm_ts.values()])
    farfarfuture = datetime(3000, 1, 1)
    
    # iterate over time instants and accumulate them into observation packets
    obs_data = {}
    while np.any(vec_ndx < len_ndx):
        
        # find current instant
        ctm = min([ s['t'][i] if i < len(s['t']) else farfarfuture for (s,i) in zip(fm_ts.values(), vec_ndx) ])
        
        # if we are past all arrays, we are done
        if ctm == farfarfuture:
            return obs_data
        
        # gather names of all stations which correspond to this time instant
        stat_names = [ s['name'] for s in filter(lambda (s,i): s['t'][i] == ctm,  zip(fm_ts.values(), vec_ndx)) ]
        
        # construct an observation packet
        obs_i = {}
        for sn in stat_names:
                 
        
    
     
        
        
    


if __name__ == '__main__':
    
    W = WRFModelData('../real_data/witch_creek/realfire03_d02_20071021.nc')
    
    # read in vars
    lat, lon = W.get_lats(), W.get_lons()
    tm = W.get_times()
    rain = v['RAINNC']
    Q2 = v['Q2']
    T2 = v['T2']
    P = v['PSFC']
    
    # obtain sizes
    times = rain.shape[0]
    dom_shape = lat.shape
    locs = np.prod(dom_shape)
    
    # load station data
    stations = load_stations_from_files(station_data_dir, station_list, 'US/Pacific')
    match_stations_to_gridpoints(stations, lon, lat)

    # construct time series from the stations that are matched to the times in the simulation
    fm_ts = match_time_series(stations, 'fuel_moisture', field, W)
    
    # manipulate observation data into a time indexed structure
    obs_data = build_observation_data(stations, fm_ts) 
    
    # construct initial vector
    Ed, Ew = equilibrium_moisture(P[0,:,:], Q2[0,:,:], T2[0,:,:])
    
    mi = np.zeros((dom_shape[0], dom_shape[1], 9))
    mi[:, :, 0:3] = 0.5 * (Ed + Ew)
    
    # set up parameters
    Tk = np.array([1.0, 10.0, 100.0]) * 3600
    dt = 10.0 * 60
    
    # observation vector
    mi1 = np.zeros_like(mi)
    for t in range(times):
        
        # run the model update
        for i in dom_shape[0]:
            for j in dom_shape[1]:
                T_i, Q_i, P_i, r_i = T2[t, i, j], Q2[t, i, j], P[t, i, j], rain[t, i, j]
                mi1[i, j, :] = moisture_model_simple(T_i, Tk, Q_i, P_i, mi[i, j, :], r_i, dt)
                
                # if we have an observation somewhere in time, first krig it to the current state
                
        

        mi[:] = mi1[:]

    # construct a basemap representation of the area
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    # extract values for different fuel types                
    m1 = mi[:, :, 0]
    m2 = mi[:, :, 1]
    m3 = mi[:, :, 2]
    
    # also plot equilibria
    Ed, Ew = equilibrium_moisture(P[-1,:,:], Q2[-1,:,:], T2[-1,:,:])
                
    plt.figure(figsize = (14,10))
    plt.subplot(221)
    render_spatial_field(m, lon, lat, m1, '1-hr fuel')
    plt.clim([0, 0.5])
    plt.colorbar()
    plt.subplot(222)
    render_spatial_field(m, lon, lat, m1, '10-hr fuel')
    plt.clim([0, 0.5])
    plt.colorbar()
    plt.subplot(223)
    render_spatial_field(m, lon, lat, m1, '100-hr fuel')
    plt.clim([0, 0.5])
    plt.colorbar()
    plt.subplot(224)
    render_spatial_field(m, lon, lat, 0.5*(Ed + Ew), 'Equilibrium moisture')
    plt.clim([0, 0.5])
    plt.colorbar()
    plt.show()
