# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 17:40:24 2012

@author: martin
"""

from spatial_model_utilities import load_wrf_data, equilibrium_moisture, render_spatial_field, load_station_data

import numpy as np
from mpl_toolkits.basemap import Basemap
import os
import pylab as pb

station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
#                  "Mt Laguna_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Mositure",
#                  "Palomar_Fuel_Moisture",
                  "Alpine_Moisture",
#                  "Valley_Center_Moisture",
                  "Ranchita_Mositure",
                  "Camp_Elliot_Moisture",
#                  "Pine Hills_Moisture" 
                ]
                  
station_data_dir = "../real_data/witch_creek/"


if __name__ == '__main__':
    
    # we will refer to these often
#    v = load_wrf_data('../real_data/california/realfire05_d03_20120409.nc')
#    v = load_wrf_data('../real_data/witch_creek/realfire03_d03_20071021.nc')
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
    tm = v['Times']

    # construct our basemap
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    # compute the equilibrium moisture on grid points
    Ed, Ew = equilibrium_moisture(v['PSFC'][:,:,:], Q2, v['T2'][:,:,:])

    pb.ion()
    pb.figure(figsize = (14, 7))
    for i in range(Ed.shape[0]):

#        if i > 0 and np.all(rain[i,:,:] == 0):
#            continue
        pb.clf()
        
        pb.subplot(221)
        render_spatial_field(m, lon, lat, Ed[i,:,:], 'Drying equilibrium [%d]' % i)
        pb.clim([np.min(Ed), np.max(Ed)])
        pb.colorbar()
        pb.subplot(222)
        render_spatial_field(m, lon, lat, Ew[i,:,:], 'Wetting equilibrium')
        pb.clim([np.min(Ew), np.max(Ew)])
        pb.colorbar()
        pb.subplot(223)
        render_spatial_field(m, lon, lat, rain[i,:,:], 'Rain')
        pb.clim([np.min(rain), np.max(rain)])
        pb.colorbar()
        pb.subplot(224)
        render_spatial_field(m, lon, lat, Q2[i,:,:], 'Water vapor ratio')
        pb.clim([np.min(Q2), np.max(Q2)])
        pb.colorbar()
        pb.draw()
        
        pb.savefig('image%03d.png' % i)

    pb.ioff()