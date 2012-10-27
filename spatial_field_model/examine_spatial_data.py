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
    dx = (np.max(lon) - np.min(lon)) / 5
    dy = (np.max(lat) - np.min(lat)) / 5
#    m.drawcoastlines()
    m.drawparallels(np.arange(lat_rng[0],lat_rng[1], dy))
    m.drawmeridians(np.arange(lon_rng[0],lon_rng[1], dx))
#    m.etopo(1.0)
    x, y = m(lon, lat)
    m.pcolormesh(x, y, field, alpha = 0.6)
    plt.axis('equal')
    plt.title(title)
#    plt.xlabel('Longitude [deg]');
#    plt.ylabel('Lattitude [deg]');


def load_data(data_file):
    """
    Load required variables from the file data_file.
    """
    v = {}
    d = netCDF4.Dataset(os.path.join(data_file))
    for vname in [ 'T2', 'Q2', 'PSFC', 'XLAT', 'XLONG', 'RAINNC' ]:
        v[vname] = d.variables[vname][:,...]
    d.close()
    return v


if __name__ == '__main__':
    
    # we will refer to these often
    v = load_data('../real_data/california/realfire05_d03_20120409.nc')
#    v = load_data('../real_data/witch_creek/realfire03_d02_20071021.nc')

    # read in vars
    lat = v['XLAT'][0,:,:]
    lon = v['XLONG'][0,:,:]
    rain = v['RAINNC']
    Q2 = v['Q2']

    # construct our basemap
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    # compute the equilibrium moisture on grid points
    Ed, Ew = equilibrium_moisture(v['PSFC'][:,:,:], Q2, v['T2'][:,:,:])

    plt.ion()
    plt.figure(figsize = (10, 6))
    for i in range(Ed.shape[0]):
        if i > 0 and np.all(rain[i,:,:] == 0):
            continue

        plt.subplot(221)
        render_field(m, lon, lat, Ed[i,:,:], 'Drying equilibrium [%d]' % i)
        plt.clim([np.min(Ed), np.max(Ed)])
        if i == 0:
            plt.colorbar()
        plt.subplot(222)
        render_field(m, lon, lat, Ew[i,:,:], 'Wetting equilibrium')
        plt.clim([np.min(Ew), np.max(Ew)])
        if i == 0:
            plt.colorbar()
        plt.subplot(223)
        render_field(m, lon, lat, rain[i,:,:], 'Rain')
        plt.clim([np.min(rain), np.max(rain)])
        if i == 0:
            plt.colorbar()
        plt.subplot(224)
        render_field(m, lon, lat, Q2[i,:,:], 'Water vapor ratio')
        plt.clim([np.min(Q2), np.max(Q2)])
        if i == 0:
            plt.colorbar()
        plt.draw()

    plt.ioff()
