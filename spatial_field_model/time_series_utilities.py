



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



def match_time_series(stations, st_field_name, field, W):
    """
    Matches the time series of the field with the values of st_field in each station in stations.
    Returns the matched time series indexed by station name.  The field must be located on the same
    grid points and times as the WRF model data in W.
    
        matched =  match_time_series(stations, st_field_name, field, W)
        
    matched is a dictionary indexed by station name with tuple values
    (match_times, model_ts, station_ts, dist_to_nearest_grid_point, index_of_nearest_grid_point). 
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
        matched[st_name] = { 't' : match_times, 'model_ts' : m_ts, 'station_ts' : st_ts, 'name' : st_name }
        
    return matched


def build_observation_data(stations, obs_type, wrf_data):
    """
    Repackage the matched time series into a time-indexed structure which gives details on the observed data and active observation stations.
    
        synopsis: obs_data = build_observation_data(stations, wrf_data)
        
    """
    Ns = len(stations)
    
    # accumulate all observations from stations
    observations = []
    for s in stations:
        observations.extend(s.get_observations_for_times(obs_type, wrf_data.get_times())[1])

    # repackage all the observations into a time-indexed structure which groups
    # observations at the same time together
    obs_data = {}
    for obs in observations:
        t = obs.get_time()
        o = obs_data[t] if obs_data.has_key(t) else []
        o.append(obs)
        obs_data[t] = o
        
    return obs_data
