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
    


def load_data(data_dir, data_list):
    
    v = {}
    for fname, vlist in data_list.iteritems():
        d = netCDF4.Dataset(os.path.join(data_dir, fname))
        for vname in vlist:
            v[vname] = d.variables[vname]
    return v



# list of variables to load
data_list = { 'wrfout_d03_T2.nc' : ['T2'],
              'wrfout_d03_Q2.nc' : ['Q2'],
              'wrfout_d03_PSFC.nc' : ['PSFC'],
              'wrfout_d03_latlon.nc' : [ 'XLAT', 'XLONG' ]}

v = load_data('../real-data', data_list)

# we will refer to these often
lat = v['XLAT'][:,:,0]
lon = v['XLONG'][:,:,0]

# compute the equilibrium moisture
Ed, Ew = equilibrium_moisture(v['PSFC'][:,:,:], v['Q2'][:,:,:], v['T2'][:,:,:])

# plot over basemap plots

lat_rng = (np.min(lat), np.max(lat))
lon_rng = (np.min(lon), np.max(lon))
m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
            urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
            projection='mill')
m.drawcoastlines()
m.drawparallels(np.arange(lat_rng[0],lat_rng[1], 0.1))
m.drawmeridians(np.arange(lon_rng[0],lon_rng[1], 0.1))

nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1
Ed0 = Ed[:, :, 0]
m.contourf(lat, lon, Ed0, latlon = True, levels = np.linspace(np.min(Ed0), np.max(Ed0), 15))

plt.figure()
plt.imshow(Ed0)
plt.title('Plain display')
plt.colorbar()

plt.show()