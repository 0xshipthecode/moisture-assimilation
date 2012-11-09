
from spatial_model_utilities import great_circle_distance
from diagnostics import diagnostics

import numpy as np


def construct_spatial_correlation_matrix(gridndx, mlons, mlats):
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



def simple_kriging_data_to_model(obs_data, obs_stds, mu_mod, wrf_data, mod_stds, t):
    """
    Simple kriging of data points to model points.  The kriging results in
    the matrix K, which contains mean of the kriged observations and the
    matrix V, which contains the kriging variance. 
    """
    mlons, mlats = wrf_data.get_lons(), wrf_data.get_lats()
    K = np.zeros_like(mlons)
    V = np.zeros_like(mlons)
    Nobs = len(obs_data)
    obs_vals = np.zeros((Nobs,))
    station_lonlat = []
    gridndx = []
    mu_obs = np.zeros((Nobs,))
    measV = np.zeros((Nobs,))
        
    # accumulate the indices of the nearest grid points
    ndx = 0
    for obs in obs_data:
        gridndx.append(obs.get_nearest_grid_point())
        obs_vals[ndx] = obs.get_value()
        station_lonlat.append(obs.get_position())
        mu_obs[ndx] = mu_mod[obs.get_nearest_grid_point()]
        measV[ndx] = obs.get_measurement_variance()
        ndx += 1

    # compute observation residuals (using model fit from examine_station_data)
    res_obs = obs_vals - mu_obs
    
    diagnostics().push("skdm_obs_res", res_obs)
    
    # construct the covariance matrix and invert it
    C = np.asmatrix(construct_spatial_correlation_matrix(gridndx, mlons, mlats))
    oS = np.asmatrix(np.diag(obs_stds))
    Sigma = oS.T * C * oS + np.diag(measV)
    SigInv = np.linalg.inv(Sigma)
    
    diagnostics().push("skdm_cov_cond", np.linalg.cond(Sigma))
    
    # run the kriging estimator for each model grid point
    K = np.zeros_like(mlats)
    cov = np.zeros_like(mu_obs)
    for p in np.ndindex(K.shape):
        # compute the covariance array anew for each grid point
        for k in range(Nobs):
            lon, lat = station_lonlat[k]
            cc = max(0.8565 - 0.0063 * great_circle_distance(mlons[p], mlats[p], lon, lat), 0.0)
            cov[k] = mod_stds[p] * cc * obs_stds[k]
        csi = np.dot(cov, SigInv)
        K[p] = mu_mod[p] + np.dot(csi, res_obs) 
        V[p] = mod_stds[p]**2 - np.dot(csi, cov)
    
    return K, V
