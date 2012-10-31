
from spatial_model_utilities import load_wrf_data, load_station_data, render_spatial_field, \
                                    equilibrium_moisture, load_stations_from_files, \
                                    match_stations_to_gridpoints, match_sample_times
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
import math


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


def match_time_series(stations, st_field_name, field, mlon, mlat, mtimes):
    """
    Matches the time series of the field with the values of st_field in each station in stations.
    Returns the matched time series indexed by station name.  The field passed in has to be valid
    for the model longitude and latitude given.
    
        matched =  match_time_series(stations, st_field_name, field, mlon, mlat, mtimes)
        
    """
    matched = {}
    for s in stations.values():
        st_fm = s[st_field_name]
        st_times = sorted(st_fm.keys())
        i, j = s['nearest_grid_point']
        m_ts = field[:, i, j]
        match_times, indx1, _ = match_sample_times(mtimes, st_times)
        m_ts = m_ts[indx1]
        st_ts = [ st_fm[d] for d in match_times ]
        matched[s['name']] = (match_times, m_ts, st_ts)
        
    return matched
        

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


if __name__ == '__main__':

    # load the smallest domain
    v = load_wrf_data('../real_data/witch_creek/realfire03_d04_20071021.nc')

    # read in vars
    lat = v['XLAT'][0,:,:]
    lon = v['XLONG'][0,:,:]
    rain = v['RAINNC']
    Q2 = v['Q2']
    T2 = v['T2']
    P = v['PSFC']
    tm = v['Times']
    
    # adjust times to match california time
    tm = [t - timedelta(0, 8 * 3600) for t in tm]
    
    # load stations and match then to grid points
    stations = load_stations_from_files(station_data_dir, station_list)
    match_stations_to_gridpoints(stations, lon, lat)

    # construct basemap for rendering
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    # compute the equilibrium moisture on grid points (for all times t)
    Ed, Ew = equilibrium_moisture(P, Q2, T2)

    # show this and render stations on top
    render_spatial_field(m, lon, lat, 0.5 * (Ed[0,:,:] + Ew[0,:,:]), 'Equilibrium moisture')
    for s in station_list:
        st = stations[s]
        slon, slat = st['lon'], st['lat']
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')

    x, y = m(lon.ravel(), lat.ravel())
    plt.plot(x, y, 'k+', markersize = 4)

    # part B, compare values at the same times
    for st_field_name, field in [ ('T', T2 - 273.15), ('rain', rain), ('fuel_moisture', 50.0 * (Ed + Ew))]:
        ms_ts = match_time_series(stations, st_field_name, field, lon, lat, tm)

        # extract station time series and corresponding grid point time series
        f = plt.figure()
        f.subplots_adjust(hspace = 1.2)
        i = 1
        for station_name in ms_ts.keys():
            match_times, m_ts, st_ts = ms_ts[station_name]
            ax = plt.subplot(4, 2, i)
            ax.xaxis.set_major_formatter(DateFormatter('%H:%m')) 
            plt.plot(match_times, st_ts, 'ro-', match_times, m_ts, 'gx-', linewidth = 2)
            plt.title('%s vs. model %s' % (station_name, st_field_name))
            for l in ax.get_xticklabels():
                l.set_rotation(90)
            i += 1
        
        
    for s in stations.values():
        i, j = s['nearest_grid_point']
        mlon = lon[i, j]
        mlat = lat[i, j]
        print("Distance of %s from nearest grid point: %g km" 
                % (s['name'], great_circle_distance(mlon, mlat, st['lon'], st['lat'])))
        
        
    plt.show()
    