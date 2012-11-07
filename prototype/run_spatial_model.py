# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:36 2012

@author: martin
"""

from spatial_model_utilities import load_wrf_data, equilibrium_moisture, render_spatial_field
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np


def moisture_model_ext(T, Tk, Q, P, m_ext, f_info, r, dt):
    """
    Compute the next state of the distributed moisture model using environmental
    parameters and the current state.
    
    synopsis: m_ext = moisture_model_ext(T, Tk, Q, P, m_ext, f_info, r, dt)
    
    T - temperature [Kelvins, location specific]
    Tk - time delay nominal values for fuel classes
    Q - water vapor ratio [dimensionless, location specific]
    P - pressure [Pascals, location specific]
    m_ext - current moisture values and assimilated parameters
    f_info - information on fuel type and location for each fuel
             in m_ext
    r - rain intensity for time unit [mm/h, location specific]
    dt - integration step [s]
    """
    
    # number of fuels modelled (not necessarily at the same spatial location)
    f_type = f_info[:,0]
    f_loc = f_info[:,1]

    k = Tk.shape[0]                 # number of fuel classes
    n = f_info.shape[0] / k         # number of fuel locations 
    r0 = 0.05                       # threshold rainfall [mm/h]
    rk = 8                          # saturation rain intensity [mm/h]
    Trk = 14 * 3600                 # time constant for wetting model [s]
    S = 2.5                         # saturation intensity [dimensionless]
    
    # number of fuel states (fuel types times number of locations)
    nk = n * k
    
    # first, we break the state vector into components
    m = m_ext[:nk]
    dlt_Tk = m_ext[nk:nk+k]
    dlt_E = m_ext[nk+k]
    dlt_S = m_ext[nk+k+1]
    dlt_Trk = m_ext[nk+k+2]
    
    # compute equilibrium moisture using model
    Ed, Ew = equilibrium_moisture(P, Q, T);
    
    # add assimilated difference, which is shared across spatial locations
    Ed = Ed + dlt_E;
    Ew = Ew + dlt_E;
    
    # where rainfall is above threshold (spatially different), apply
    # saturation model, equi and rlag are specific to fuel type and
    # location
    equi = m.copy()         # copy over current equilibrium levels
    rlag = np.zeros((nk,))
    model_ids = np.zeros((nk,))
    
    # equilibrium is equal to the saturation level (assimilated)
    rain_model = r[f_loc] > r0
    equi[rain_model] = S + dlt_S
    model_ids[rain_model] = 3
    
    # rlag is modified by the rainfall intensity (assimilated)
    rlag[rain_model] = 1.0 / (Trk + dlt_Trk) * (1 - np.exp(- (r[f_loc[rain_model]] - r0) / rk))

    # equilibrium is selected according to current moisture level
    sel = np.logical_and(np.logical_not(rain_model), equi > Ed[f_loc])
    equi[sel] = Ed[f_loc][sel]
    sel = np.logical_and(np.logical_not(rain_model), equi < Ew[f_loc])
    equi[sel] = Ew[f_loc][sel]
    notrain = np.logical_not(rain_model)
    model_ids[notrain] = (m[notrain] > Ed[f_loc[notrain]]) * 1 + (m[notrain] < Ew[f_loc[notrain]]) * 2;
    model_ids[model_ids == 0] = 4;

    # the inverted time lag is constant according to fuel category
    rlag[notrain] = 1.0 / (Tk[f_type[notrain]] + dlt_Tk[f_type[notrain]])
    
    # select appropriate integration method according to change for each fuel 
    # and location
    change = dt * rlag
    m_new = np.zeros_like(m_ext)
    m_new[:nk] = np.where(change < 0.01,
                          m + (equi - m) * (1 - np.exp(-change)),
                          m + (equi - m) * change * (1 - 0.5 * change))
    m_new[nk:] = m_ext[nk:]
    return m_new
    
    
if __name__ == '__main__':
    
    v = load_wrf_data('../real_data/witch_creek/realfire03_d02_20071021.nc')

    # read in vars
    lat = v['XLAT'][0,:,:]
    lon = v['XLONG'][0,:,:]
    rain = v['RAINNC']
    Q2 = v['Q2']
    T2 = v['T2']
    P = v['PSFC']
    
    # obtain sizes
    times = rain.shape[0]
    locs = np.prod(lat.shape)
    
    # construct initial vector
    Ed,Ew = equilibrium_moisture(P[0,:,:], Q2[0,:,:], T2[0,:,:])
    m0l = 0.5 * (Ed + Ew)
    domain_shape = m0l.shape
    m0l = m0l.ravel()
    
    # inform simulation of location & fuel type structure
    f_info = np.zeros((3*locs,2), dtype = int)
    f_info[:,0] = [ i % 3 for i in range(locs * 3) ]
    f_info[:,1] = [ i // 3 for i in range(locs * 3) ]

    # initialize each fuel type at each location with same value
    m0 = m0l[f_info[:,1]]
    
    # construct initial extended state
    miext = np.r_[ m0, np.zeros((6,)) ]
    
    # set up parameters
    Tk = np.array([1.0, 10.0, 100.0]) * 3600
    dt = 10.0 * 60
    
    for i in range(times):
        T_i, Q_i, P_i, r_i = T2[i,:,:], Q2[i,:,:], P[i,:,:], rain[i,:,:]
        mi1ext = moisture_model_ext(T_i.ravel(), Tk, Q_i.ravel(), P_i.ravel(), miext, f_info, r_i.ravel(), dt)
        miext[:] = mi1ext


    # construct a basemap representation of the area
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    # extract values for different fuel types                
    m1 = mi1ext[[i*3 for i in range(locs)]]
    m1 = m1.reshape(domain_shape)
    m2 = mi1ext[[i*3+1 for i in range(locs)]]
    m2 = m2.reshape(domain_shape)
    m3 = mi1ext[[i*3+1 for i in range(locs)]]
    m3 = m3.reshape(domain_shape)
    
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
