# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 17:40:24 2012

@author: martin
"""

from spatial_model_utilities import load_wrf_data, equilibrium_moisture, render_spatial_field

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt



if __name__ == '__main__':
    
    # we will refer to these often
    v = load_wrf_data('../real_data/california/realfire05_d03_20120409.nc')
#    v = load_wrf_data('../real_data/witch_creek/realfire03_d02_20071021.nc')

    # read in vars
    lat = v['XLAT'][0,:,:]
    lon = v['XLONG'][0,:,:]
    rain = v['RAINNC']
    Q2 = v['Q2']
    T2 = v['T2']
    P = v['PSFC']

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
        render_spatial_field(m, lon, lat, Ed[i,:,:], 'Drying equilibrium [%d]' % i)
        plt.clim([np.min(Ed), np.max(Ed)])
        if i == 0:
            plt.colorbar()
        plt.subplot(222)
        render_spatial_field(m, lon, lat, Ew[i,:,:], 'Wetting equilibrium')
        plt.clim([np.min(Ew), np.max(Ew)])
        if i == 0:
            plt.colorbar()
        plt.subplot(223)
        render_spatial_field(m, lon, lat, rain[i,:,:], 'Rain')
        plt.clim([np.min(rain), np.max(rain)])
        if i == 0:
            plt.colorbar()
        plt.subplot(224)
        render_spatial_field(m, lon, lat, Q2[i,:,:], 'Water vapor ratio')
        plt.clim([np.min(Q2), np.max(Q2)])
        if i == 0:
            plt.colorbar()
        plt.draw()

    plt.ioff()
