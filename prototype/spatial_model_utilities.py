# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:14:57 2012

@author: martin
"""

import numpy as np
import os.path
import matplotlib.pyplot as plt
import codecs
import re
from datetime import datetime, timedelta
import unicodedata
import pytz
import math



def great_circle_distance(lon1, lat1, lon2, lat2):
    """
    Computes the great circle distance between two points given as (lon1,lat1), (lon2,lat2)
    in kilometers.
    
        d = great_circle_distance(lon1, lat1, lon2, lat2)
    """
    rlat1, rlat2 = np.pi * lat1 / 180.0, np.pi * lat2 / 180.0
    rlon1, rlon2 = np.pi * lon1 / 180.0, np.pi * lon2 / 180.0
    
    a = math.sin(0.5*(rlat1 - rlat2))**2 + math.cos(rlat1)*math.cos(rlat2)*math.sin(0.5*(rlon1 - rlon2))**2
    c = 2 * math.atan2(a**0.5, (1-a)**0.5)
    return 6371.0 * c



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
    


def render_spatial_field(m, lon, lat, field, title):
    """
    Render a geo field over a basemap.  However the standard basemap
    does not seem to have sufficient resolution.
    TODO: can I obtain more detailed topography?
    """
    dx = (np.max(lon) - np.min(lon)) / 5
    dy = (np.max(lat) - np.min(lat)) / 5
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m.drawparallels(np.arange(lat_rng[0],lat_rng[1], dy))
    m.drawmeridians(np.arange(lon_rng[0],lon_rng[1], dx))
    m.drawcoastlines()
#    imf = m.transform_scalar(field.ravel(), lon.ravel(), lat.ravel(), lon.shape[0], lat.shape[1])
    x, y = m(lon, lat)
    m.pcolormesh(x, y, field, alpha = 0.6, edgecolor = 'none')
#    m.imshow(imf)
    plt.axis('equal')
    plt.title(title)
    


def find_closest_grid_point(slon, slat, glon, glat):
    """
    Finds the closest grid point to the given station longitude/lattitude.
    """
    closest = np.argmin((slon - glon)**2 + (slat - glat)**2)
    return np.unravel_index(closest, glon.shape)



def match_stations_to_gridpoints(sts, lon, lat):
    """
    Finds the nearest grid point for each station and stores it in the station dictionary.
    The nearest grid point is stored as the value of 'nearest_grid_point'.
    """
    for s in sts.values():
        i, j = find_closest_grid_point(s['lon'], s['lat'], lon, lat)
        d = great_circle_distance(s['lon'], s['lat'], lon[i,j], lat[i,j])
        s['nearest_grid_point'] = (i,j)
        s['dist_nearest_grid_point'] = d
        


def load_stations_from_files(station_data_dir, station_list, tz_name, W):
    """
    Loads station data from files given to me by AK and stores all station data
    in a dictionary indexed by station file names.  All the times in the files
    are considered to be relative to the time zone tz_name.
    FIXME: assuming all stations are in the same time zone is dangerous! 
    
        stations = load_stations_from_files(station_data_dir, station_list, tz_name)
    
    """
    tz = pytz.timezone(tz_name)
    return [Station(os.path.join(station_data_dir, s), tz, W) for s in station_list]


