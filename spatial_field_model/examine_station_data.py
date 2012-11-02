
from spatial_model_utilities import load_station_data, render_spatial_field, \
                                    equilibrium_moisture, load_stations_from_files, \
                                    match_stations_to_gridpoints, match_sample_times, \
                                    great_circle_distance
                                    
from wrf_model_data import WRFModelData

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter


station_list = [  "Julian_Moisture",
                  "Goose_Valley_Fuel_Moisture",
#                  "Mt Laguna_Moisture",
                  "Oak_Grove_Moisture",
                  "Decanso_Moisture",
#                  "Palomar_Fuel_Moisture",
                  "Alpine_Moisture",
#                  "Valley_Center_Moisture",
                  "Ranchita_Moisture",
                  "Camp_Elliot_Moisture",
#                  "Pine Hills_Moisture" 
                ]
                  
station_data_dir = "../real_data/witch_creek/"


def match_time_series(stations, st_field_name, field, W):
    """
    Matches the time series of the field with the values of st_field in each station in stations.
    Returns the matched time series indexed by station name.  The field must be located on the same
    grid points and times as the WRF model data in W.
    
        matched =  match_time_series(stations, st_field_name, field, W)
        
    """
    mlon = W.get_lons()
    mlat = W.get_lats()
    mtimes = W.get_times()
    
    matched = {}
    for st_name in stations:
        s = stations[st_name]
        st_data = s[st_field_name]
        st_times = sorted(st_data.keys())
        i, j = s['nearest_grid_point']
        m_ts = field[:, i, j]
        match_times, indx1, _ = match_sample_times(mtimes, st_times)
        m_ts = m_ts[indx1]
        st_ts = [ st_data[d] for d in match_times ]
        matched[st_name] = (match_times, m_ts, st_ts, great_circle_distance(mlon[i,j], mlat[i,j], s['lon'], s['lat']))
        
    return matched


def plot_stations_vs_model_ts(stations, field_name, field, W):
    
    ms_ts = match_time_series(stations, st_field_name, field, W)

    # extract station time series and corresponding grid point time series
    f = plt.figure()
    f.subplots_adjust(hspace = 1.2)
    i = 1
    for station_name in ms_ts.keys():
        match_times, m_ts, st_ts, _ = ms_ts[station_name]
        ax = plt.subplot(4, 2, i)
        ax.xaxis.set_major_formatter(DateFormatter('%H:%m')) 
        plt.plot(match_times, st_ts, 'ro-', match_times, m_ts, 'gx-', linewidth = 2)
        plt.title('%s vs. model %s' % (station_name, st_field_name))
        for l in ax.get_xticklabels():
            l.set_rotation(90)
        i += 1
    


if __name__ == '__main__':

    # load the smallest domain
    W = WRFModelData('../real_data/witch_creek/realfire03_d04_20071021.nc',
                     [ 'T2', 'PSFC', 'Q2', 'RAINNC' ], 'US/Pacific')

    # read in vars
    lat = W.get_lats()
    lon = W.get_lons()
    rain = W['RAINNC']
    Q2 = W['Q2']
    T2 = W['T2']
    P = W['PSFC']
    tm = W.get_times()
        
    # load stations and match then to grid points
    stations = load_stations_from_files(station_data_dir, station_list, 'US/Pacific')
    match_stations_to_gridpoints(stations, lon, lat)

    # construct basemap for rendering
    domain_rng = W.get_domain_extent()
    m = Basemap(llcrnrlon=domain_rng[0],llcrnrlat=domain_rng[1],
                urcrnrlon=domain_rng[2],urcrnrlat=domain_rng[3],
                projection = 'mill')

    # compute the equilibrium moisture on grid points (for all times t)
    Ed, Ew = equilibrium_moisture(P, Q2, T2)

    # show the equilibrium field and render position of stations on top
    render_spatial_field(m, lon, lat, 0.5 * (Ed[0,:,:] + Ew[0,:,:]), 'Equilibrium moisture')
    for s in stations.keys():
        st = stations[s]
        slon, slat = st['lon'], st['lat']
        x, y = m(slon, slat)
        plt.plot(x, y, 'o', markersize = 8, markerfacecolor = 'magenta')
        plt.text(x, y, s)

    x, y = m(lon.ravel(), lat.ravel())
    plt.plot(x, y, 'k+', markersize = 4)

    # part B, compare values at the same times
    for st_field_name, field in [ ('T', T2 - 273.15), ('fuel_moisture', 50.0 * (Ed + Ew))]:
        plot_stations_vs_model_ts(stations, st_field_name, field, W)
        
    # examine fuel moisture residuals
    fm_ts  = match_time_series(stations, 'fuel_moisture', field, W)
    
    # recode all of this into a matrix to fit the model
    model_data = []
    stat_data = []
    weights = []
    for s in fm_ts.keys():
        _, mts, stts, d = fm_ts[s]
        model_data.extend(mts)
        stat_data.extend(stts)
        weights.extend([1.0 / d for i in range(len(mts))])
        
    mdata = np.array(model_data)
    sdata = np.array(stat_data)
    weights = np.array(weights)
    
    # compute the best fit of the data
    beta = np.sum(weights * mdata * sdata) / np.sum(weights * mdata ** 2)
    print("Weighted lin. regression parameter is %g." % beta)
    
    # plot fitted model
    plot_stations_vs_model_ts(stations, 'fuel_moisture', beta * 50 * (Ed + Ew), W)
    
    # now get all the residuals of the fit in each station
    sres = {}
    for s in fm_ts:
        print s
        st = stations[s]
        _, mts, stts, _ = fm_ts[s]
        sres[s] = (np.array(stts) - beta * np.array(mts), st['lon'], st['lat'])

    # **********************************************************************************
    # compute COVARIANCE between station residuals and plot this vs. distance
    Ns = len(fm_ts.keys())
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    for s1, i in zip(fm_ts.keys(), range(Ns)):
        r1, lon1, lat1 = sres[s1]
        for s2, j in zip(fm_ts.keys(), range(Ns)):
            r2, lon2, lat2 = sres[s2]
            cc = np.cov(r1, r2)
            C[i,j] = cc[0,1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)

    f = plt.figure(figsize = (8,8))
    f.subplots_adjust(hspace = 0.5, wspace = 0.5)
    ax = plt.subplot(221)
    plt.imshow(C, interpolation = 'nearest')
    plt.title('Covariance [-]')
    plt.colorbar()
    ax.set_xticklabels(fm_ts.keys(), rotation = 90)
    ax.set_yticklabels(fm_ts.keys())
    ax = plt.subplot(222)
    plt.imshow(D, interpolation = 'nearest')
    plt.title('Distance [km]')
    ax.set_xticklabels(fm_ts.keys(), rotation = 90)
    ax.set_yticklabels(fm_ts.keys())
    plt.colorbar()
    plt.subplot(223)
    colors = 'rgbcmyk'
    for r in range(D.shape[0]):
        ndx = np.setdiff1d(np.arange(Ns), [r])
        plt.plot(D[r,ndx], C[r,ndx], '%so' % colors[r])
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Covariance [-]')
    plt.title('Aggregate plot of covar vs. dist')

    plt.subplot(224)
    iu1 = np.triu_indices(Ns, 1)
    Dt = D[iu1]
    Ct = C[iu1]
    plt.plot(Dt, Ct, 'ro')
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Covariance [-]')
    plt.title('Aggregate plot of covar vs. distance')

    f = open('distance_vs_covariance.txt', 'w')
    for i in range(Dt.shape[0]):
        f.write('%g,%g\n' % (Dt[i], Ct[i]))
    f.close() 


    # **********************************************************************************
    # compute CORRELATION COEFFICIENT between station residuals and plot this vs. distance
    Ns = len(fm_ts.keys())
    C = np.zeros((Ns,Ns))
    D = np.zeros((Ns,Ns))
    for s1, i in zip(fm_ts.keys(), range(Ns)):
        r1, lon1, lat1 = sres[s1]
        for s2, j in zip(fm_ts.keys(), range(Ns)):
            r2, lon2, lat2 = sres[s2]
            cc = np.corrcoef(r1, r2)
            C[i,j] = cc[0,1]
            D[i,j] = great_circle_distance(lon1, lat2, lon2, lat2)

    f = plt.figure(figsize = (8,8))
    f.subplots_adjust(hspace = 0.5, wspace = 0.5)
    ax = plt.subplot(221)
    plt.imshow(C, interpolation = 'nearest')
    plt.clim([0.0, 1.0])
    plt.title('Correlation coefficient [-]')
    plt.colorbar()
    ax.set_xticklabels(fm_ts.keys(), rotation = 90)
    ax.set_yticklabels(fm_ts.keys())
    ax = plt.subplot(222)
    plt.imshow(D, interpolation = 'nearest')
    plt.title('Distance [km]')
    ax.set_xticklabels(fm_ts.keys(), rotation = 90)
    ax.set_yticklabels(fm_ts.keys())
    plt.colorbar()
    plt.subplot(223)
    colors = 'rgbcmyk'
    for r in range(D.shape[0]):
        ndx = np.setdiff1d(np.arange(Ns), [r])
        plt.plot(D[r,ndx], C[r,ndx], '%so' % colors[r])
    plt.ylim([0.0, 1.0])
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation coefficient [-]')
    plt.title('Aggregate plot of cc vs. dist')

    plt.subplot(224)
    iu1 = np.triu_indices(Ns, 1)
    Dt = D[iu1]
    Ct = C[iu1]
    plt.plot(Dt, Ct, 'ro')
    plt.plot(Dt, 0.8166 - 0.0052 * Dt, 'r-')
    plt.plot(Dt, 0.8565 - 0.0063 * Dt, 'g-')
    plt.plot(Dt, 0.8729 - 0.0068 * Dt, 'b-')
    plt.ylim([0.0, 1.0])
    plt.legend(['Data', 'OLS', 'Huber', 'Bisquare'], loc = 'lower left', prop = { 'size' : 10 } )
    plt.grid()
    plt.xlabel('Distance [km]')
    plt.ylabel('Correlation coefficient [-]')
    plt.title('Aggregate plot of cc vs. distance')
    
    f = open('distance_vs_covariance.txt', 'w')
    for i in range(Dt.shape[0]):
        f.write('%g,%g\n' % (Dt[i], Ct[i]))
    f.close() 
    
    plt.show()
    