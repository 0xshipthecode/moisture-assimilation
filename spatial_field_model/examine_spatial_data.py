# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 17:40:24 2012

@author: martin
"""



import numpy as np
import os.path
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def equilibrium_moisture(P, Q, T):

    # saturated vapor pressure (at each location, size n x 1)
    Pws = np.exp(54.842763 - 6763.22/T - 4.210 * np.log(T) + 0.000367*T + np.tanh(0.0415*(T - 218.8)) 
        * (53.878 - 1331.22/T - 9.44523 * np.log(T) + 0.014025*T))
      
    # water vapor pressure (at each location, size n x 1)
    Pw = P * Q / (0.622 + (1 - 0.622) * Q)
    
    # relative humidity (percent, at each location, size n x 1)
    H = 100 * Pw / Pws;
    
    # drying/wetting fuel equilibrium moisture contents (location specific,
    # n x 1)
    Ed = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    Ew = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))

    # remap values
    return Ed * 0.01, Ew * 0.01
    

def render_field(m, lon, lat, field, title):
    """
    Render a geo field over the map.
    """
    m.drawparallels(np.arange(lat_rng[0],lat_rng[1], 0.2))
    m.drawmeridians(np.arange(lon_rng[0],lon_rng[1], 0.2))
    x,y = m(lon, lat)
    m.pcolormesh(x, y, field)
    plt.axis('equal')
    plt.title(title)
    plt.xlabel('Longitude [deg]');
    plt.ylabel('Lattitude [deg]');


def load_data():
    
    # list of variables to load
    data_list = { 'wrfout_d03_T2.nc' : ['T2'],
                  'wrfout_d03_Q2.nc' : ['Q2'],
                  'wrfout_d03_PSFC.nc' : ['PSFC'],
                  'wrfout_d03_latlon.nc' : [ 'XLAT', 'XLONG' ],
                  'wrfout_d03_rainnc.nc' : ['RAINNC']}
    
    data_dir = '../real_data'

    v = {}
    for fname, vlist in data_list.iteritems():
        d = netCDF4.Dataset(os.path.join(data_dir, fname))
        for vname in vlist:
            v[vname] = d.variables[vname]
    return v


if __name__ == '__main__':
    
    # we will refer to these often
    v = load_data()
    lat = v['XLAT'][0,:,:]
    lon = v['XLONG'][0,:,:]
    rain = v['RAINNC'][:,:,:]

    # construct our basemap
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')
    

    # compute the equilibrium moisture on grid points
    Ed, Ew = equilibrium_moisture(v['PSFC'][:,:,:], v['Q2'][:,:,:], v['T2'][:,:,:])

    plt.ion()
    plt.figure(figsize = (6, 12))
        
    for i in range(Ed.shape[0]):
        plt.subplot(311)
        render_field(m, lon, lat, Ed[i,:,:], 'Drying equilibrium')
        if i == 0:
            plt.colorbar()
        plt.subplot(312)
        render_field(m, lon, lat, Ew[i,:,:], 'Wetting equilibrium')
        if i == 0:
            plt.colorbar()
        plt.subplot(313)
        render_field(m, lon, lat, rain, 'Rain')
        if i == 0:
            plt.colorbar()
        plt.draw()

    plt.ioff()