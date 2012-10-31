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



def load_station_data(station_file, tz):
    """
    Load all available fuel moisture data from the station information.
    """
    f = codecs.open(station_file, 'r', encoding = 'utf-8')
    station_data = {}
    station_data['name'] = f.readline().strip()

    # next line is location string, not interesting    
    f.readline()
    
    # next 2 lines are lattitude & longitude
    lat_lon_re = re.compile("\D+(\d{2,3})\D+(\d{2})\D+(\d{2})")

    l = f.readline()
    mo = lat_lon_re.match(l)
    lat_info = map(lambda x: int(x), mo.groups())
    station_data['lat'] = lat_info[0] + lat_info[1] / 60.0 + lat_info[2] / 3600.0

    l = f.readline()
    mo = lat_lon_re.match(l)
    lon_info = map(lambda x: int(x), mo.groups())
    station_data['lon'] = -(lon_info[0] + lon_info[1] / 60.0 + lon_info[2] / 3600.0)
    
    # read lines 5 through 14
    for i in range(5,8):
        f.readline()

    moisture = {}
    rh = {}
    rain = {}
    temp = {}
    fuel_temp = {}
    while True:

        # read in and parse date
        l = f.readline()
        date = datetime.strptime(l.strip(), '%B %d, %Y').replace(tzinfo = tz)
        
        # read lines until a line starts with daily
        while l[0] < '0' or l[0] > '9' and len(l) > 0:
            l = f.readline()
            
        if len(l) == 0:
            break
        
        while l[0] >= '0' and l[0] <= '9' and len(l) > 0:
            fields = filter(lambda x: len(x) > 0, l.split('\t'))
            time = datetime.strptime(fields[0], "%I %p")
            timed = timedelta(0, time.hour * 3600)
            mtime = date + timed
            if len(fields) != 12:
                print fields
            temp[mtime] = float(fields[5])
            fuel_temp[mtime] = float(fields[6])
            moisture[mtime] = float(fields[7])
            rh[mtime] = float(fields[8])
            rain[mtime] = float(fields[11])
            l = f.readline()
        
        while l[:5] != 'Daily' and len(l) > 0:
            l = f.readline()
            
        if len(l) == 0:
            break
                
    station_data['fuel_moisture'] = moisture
    station_data['relative_humidity'] = rh
    station_data['rain'] = rain
    station_data['T'] = temp
    station_data['fuel_T'] = fuel_temp
    
    f.close()
    
    return station_data
    


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
        


def load_stations_from_files(station_data_dir, station_list, tz_name):
    """
    Loads station data from files given to me by AK and stores all station data
    in a dictionary indexed by station file names.  All the times in the files
    are considered to be relative to the time zone tz_name.
    FIXME: assuming all stations are in the same time zone is dangerous! 
    
        stations = load_stations_from_files(station_data_dir, station_list, tz_name)
    
    """
    stations = {}
    tz = pytz.timezone(tz_name) 
    for s in station_list:
        st = load_station_data(os.path.join(station_data_dir, s), tz)
        if s.endswith('_Moisture'):
            s = s[:-9]
        if s.endswith('_Fuel'):
            s = s[:-5]
        stations[s] = st

    return stations



def match_sample_times(tm1, tm2):
    """
    Match times assuming both times are sorted datetime arrays.  Returns
    the matching times and the indices of the matching times in the first
    and in the second array.
    
       isect, indx1, indx2 = match_sample_times(tm1, tm2) 
        
    """
    i, j = 0, 0
    isect = []
    indx1 = []
    indx2 = []
    while i < len(tm1) and j < len(tm2):
        if tm1[i] < tm2[j]:
            i += 1
        elif tm1[i] > tm2[j]:        
            j += 1
        else:
            isect.append(tm1[i])
            indx1.append(i)
            indx2.append(j)
            i += 1
            j += 1
            
    return isect, indx1, indx2
