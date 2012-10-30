
from spatial_model_utilities import load_wrf_data, load_station_data, render_spatial_field, equilibrium_moisture
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

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


def find_closest_grid_point(slon, slat, glon, glat):
    
    closest = np.argmin((slon - glon)**2 + (slat - glat)**2)
    return np.unravel_index(closest, glon.shape)


def match_times(tm1, tm2):
    """
    Match times assuming both are sorted.
    """
    i, j = 0, 0
    isect = []
    indx = []
    while i < len(tm1) and j < len(tm2):
        while i < len(tm1) and tm1[i] < tm2[j]:
            i += 1
        while j < len(tm2) and i < len(tm1) and tm1[i] > tm2[j]:
            j += 1
        if i < len(tm1) and j < len(tm2) and tm1[i] == tm2[j]:
            isect.append(tm1[i])
            indx.append(i)
            i += 1
            j += 1
            
    return isect, indx

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
    
    stations = {}
    for s in station_list:
        print("Loading station: " + s)
        st = load_station_data(os.path.join(station_data_dir, s))
        st['closest_gridpt'] = find_closest_grid_point(st['lon'], st['lat'], lon, lat)
        stations[s] = st

    # construct our basemap
    lat_rng = (np.min(lat), np.max(lat))
    lon_rng = (np.min(lon), np.max(lon))
    m = Basemap(llcrnrlon=lon_rng[0],llcrnrlat=lat_rng[0],
                urcrnrlon=lon_rng[1],urcrnrlat=lat_rng[1],
                projection = 'mill')

    # compute the equilibrium moisture on grid points
    Ed, Ew = equilibrium_moisture(P, Q2, T2)

    # show this and render stations on top
    render_spatial_field(m, lon, lat, 0.5 * (Ed[0,:,:] + Ew[0,:,:]), 'Equilibrium moisture')
    for s in station_list:
        st = stations[s]
        slon, slat = st['lon'], st['lat']
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')
        
#    x, y = m(lon.ravel(), lat.ravel())
#    plt.plot(x, y, 'k+', markersize = 4)
    
    # part B, compare values at the same time
    plt.figure()
    for s, ndx in zip(station_list, range(len(station_list))):
        st = stations[s]
        i, j = st['closest_gridpt']
        st_fm = st['fuel_moisture']
#        st_fm = st['rain']
        isect, indx = match_times(tm, sorted(st_fm.keys()))
        s_fm = [ st_fm[d] / 100.0 for d in isect ]
        wrf_fm = [ 0.5 * (Ed[t,i,j] + Ew[t,i,j]) for t in indx ]
#        wrf_fm = [ rain[t,i,j] for t in indx ]

        plt.subplot(4,2,ndx+1)
        plt.plot(indx, s_fm, 'ro-', indx, wrf_fm, 'gx-', linewidth = 2)
        plt.title('%s vs. grid point (%d,%d) equilibrium' % (st['name'], i, j))
        
    plt.show()
    