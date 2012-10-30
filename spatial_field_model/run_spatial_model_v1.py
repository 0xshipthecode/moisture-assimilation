# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:36 2012

@author: martin
"""

from spatial_model_utilities import load_wrf_data, equilibrium_moisture, \
                                    render_spatial_field, load_station_data
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os


station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
                  "Mt Laguna_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Mositure",
                  "Palomar_Fuel_Moisture",
                  "Alpine_Moisture",
                  "Valley_Center_Moisture",
                  "Ranchita_Mositure",
                  "Camp_Elliot_Moisture",
                  "Pine Hills_Moisture" ]
                  
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
    
    
if __name__ == '__main__':
    
    v = load_wrf_data('../real_data/witch_creek/realfire03_d02_20071021.nc')
    
    stations = {}
    for s in station_list:
        print("Loading station: " + s)
        stations[s] = load_station_data(os.path.join(station_data_dir, s))

    # read in vars
    lat = v['XLAT'][0,:,:]
    lon = v['XLONG'][0,:,:]
    rain = v['RAINNC']
    Q2 = v['Q2']
    T2 = v['T2']
    P = v['PSFC']
    
    # obtain sizes
    times = rain.shape[0]
    dom_shape = lat.shape
    locs = np.prod(dom_shape)
    
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
                
        # if we have an observation somewhere, first krig it to the current state
        

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
