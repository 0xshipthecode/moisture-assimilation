#!/usr/bin/env julia

#
#  This is the main file for data assimilation runs.
#
#
# 1. loads up a configuration file,
# 2. obtains data from a WRF model,
# 3. construct covariate vectors
# 3. reads in observations and metadata for a list of stations,
# 4. runs the moisture model and the assimilation mechanism.
#
#

using Calendar
import Calendar.CalendarTime

using Storage
import Storage.setup_tag, Storage.spush, Storage.next_frame, Storage.flush_frame

using Stations
import Stations.Station, Stations.Observation, Stations.load_station_info,
       Stations.load_station_data, Stations.build_observation_data, Stations.register_to_grid, Stations.nearest_grid_point, Stations.obs_value, Stations.obs_station_id

using Kriging
import Kriging.trend_surface_model_kriging

using WRF

using FM
import FM.FMModel, FM.advance_model, FM.kalman_update


function main(args)

    # the arguments passed to the julia program do not include program name
    if length(args) != 1
        println("Usage: julia run_data_assimilation.jl cfg_file")
        exit(1)
    end
    
    ### Read configuration file and setup the system

    cfg = evalfile(args[1])

    # configure Storage mechanism
    Storage.sopen(cfg["output_dir"], "moisture_model_v2_diagnostics.txt", "frame")

    # setup Storage & output policies for interesting quantities
    setup_tag("mt", true, true, true)

    setup_tag("fm10_model_state", false, false, false)
    setup_tag("fm10_model_na_state", false, false, false)
    setup_tag("fm10_model_var", false, false, false)

    setup_tag("kriging_beta", true, true, true)
    setup_tag("kriging_xtx_cond", true, true, true)
    setup_tag("kriging_field", false, false, false)
    setup_tag("kriging_variance", false, false, false)

    setup_tag("model_raws_mae", true, true, true)
    setup_tag("model_na_raws_mae", true, true, true)

    # co-located model/model_na/kriging field/observation
    setup_tag("kriging_obs", false, false, false)
    setup_tag("kriging_obs_station_ids", false, false, false)
    setup_tag("kriging_obs_ngp", false, false, false)

    ### Load WRF output data

    # read in data from the WRF output file pointed to by cfg
    w = WRF.load_wrf_data(cfg["wrf_output"], ["HGT"])

    # the terrain height need not be stored for all time points
    WRF.slice_field(w, "HGT")

    # extract WRF fields
    lat, lon = WRF.lat(w), WRF.lon(w)
    wtm = WRF.times(w)
    dsize = size(lat)

    # retrieve equilibria and rain
    Ed, Ew = WRF.field(w, "Ed"), WRF.field(w, "Ew")
    rain = WRF.field(w, "RAIN")
    hgt = WRF.field(w, "HGT")

    ### Load observation data from stations
    io = open(join([cfg["station_data_dir"], cfg["station_list_file"]], "/"), "r")
    station_ids = filter(x -> x[1] != '#', map(x -> strip(x), readlines(io)))
    close(io)

    # load each station from its info and observation files
    stations = Station[]
    for sid in station_ids
        s = load_station_info(join([cfg["station_data_dir"], string(sid, ".info")], "/"))
        load_station_data(s, join([cfg["station_data_dir"], string(sid, ".obs")], "/"))
	register_to_grid(s, lat, lon)
        push!(stations, s)
    end

    # build the observation data from stations
    obs_fm10 = build_observation_data(stations, "FM")
    obs_times = keys(obs_fm10)

    ### Initialize model

    # maximum equilibrium moisture for visualization
    maxE = cfg["maxE"]

    # construct initial conditions (FIXME: can we do better here?)
    E = squeeze(0.5 * (Ed[2,:,:] + Ew[2,:,:]), 1)

    # set up parameters
    Q = eye(9) * cfg["Q"]
    P0 = eye(9) * cfg["P0"]
    dt = (wtm[2] - wtm[1]).millis / 1000
    mV = zeros(Float64, dsize)
    pred = zeros(Float64, dsize)
    mresV = zeros(Float64, dsize)
    mid = zeros(Int32, dsize)
    Kg = zeros(Float64, (dsize[1], dsize[2], 9))

    println("INFO: time step from WRF is $dt s.")

    # build static part of covariates
    covar_ids = cfg["static_covariates"]
    covar_map = [ :lon => lon, :lat => lat, :elevation => hgt, :constant => ones(Float64, dsize) ]
    Xd3 = length(covar_ids) + 1                         
    X = zeros(Float64, (dsize[1], dsize[2], Xd3))
    for i in 2:Xd3
        v = covar_map[covar_ids[i-1]]
        if covar_ids[i-1] != :constant
            X[:,:,i] = v - mean(v)  # ensure zero mean for each covariate (except for the constant)
        end
    end
    println("INFO: there are $Xd3 covariates (including model state).")

    # construct model grid from fuel parameters
    Tk = [ 1.0, 10.0, 100.0 ]
    models = Array(FMModel, size(E))
    models_na = Array(FMModel, size(E))
    for i in 1:dsize[1]
    	for j in 1:dsize[2]
            geo_loc = (lat[i,j], lon[i,j])
	    models[i,j] = FMModel(geo_loc, 3, E[i,j], P0, Tk)
	    models_na[i,j] = FMModel(geo_loc, 3, E[i,j], P0, Tk)
	end
    end

    ###  Run the model
    for t in 2:length(wtm)
    	mt = wtm[t]
        spush("mt", mt)

        # run the model update
	for i in 1:dsize[1]
	    for j in 1:dsize[2]
	        advance_model(models[i,j], Ed[t-1, i, j], Ew[t-1, i, j], rain[t-1, i, j], dt, Q)
		advance_model(models_na[i,j], Ed[t-1, i, j], Ew[t-1, i, j], rain[t-1, i, j], dt, Q)
	    end
	end

        # store the model state in an array (and store in output frame)
        fm10_model_state = [ models[i,j].m_ext[2] for i=1:dsize[1], j=1:dsize[2] ]
        fm10_model_na_state = [ models_na[i,j].m_ext[2] for i=1:dsize[1], j=1:dsize[2] ]
        fm10_model_var = [ models[i,j].P[2,2] for i=1:dsize[1], j=1:dsize[2] ]

        spush("fm10_model_state", fm10_model_state)
        spush("fm10_model_na_state", fm10_model_na_state)
        spush("fm10_model_var", fm10_model_var)

        # if observation data for this timepoint is available
        obs_i = Observation[]
        tm_valid_now = filter(x -> abs((mt - x).millis) / 1000.0 < 30*60, obs_times)
        if length(tm_valid_now) > 0

            # gather all observations
            for t in tm_valid_now append!(obs_i, obs_fm10[t]) end

            # set the current fm10 model state as the covariate
            X[:,:,1] = fm10_model_state

            # scale each covariate to have approximately the same norm as fm10
            # to improve condition number of X'*X
            s = sum(X[:,:,1].^2)^0.5
            for i in 2:Xd3
                X[:,:,i] = X[:,:,i] / sum(X[:,:,i].^2)^0.5 * s
            end

            # store diagnostic information
            ngp_list = map(x -> nearest_grid_point(x), obs_i)
            stat_ids = map(x -> obs_station_id(x), obs_i)
            m_at_obs = Float64[X[p[1], p[2], 1] for p in  ngp_list]
            m_na_at_obs = Float64[models_na[p[1], p[2]].m_ext[2] for p in ngp_list]
            raws = Float64[obs_value(o) for o in obs_i]

            spush("model_raws_mae", mean(abs(m_at_obs - raws)))
            spush("model_na_raws_mae", mean(abs(m_na_at_obs - raws)))

            spush("kriging_obs", raws)
            spush("kriging_obs_station_ids", stat_ids)
            spush("kriging_obs_ngp", ngp_list)

            # compute the kriging estimate
            K, V, y = trend_surface_model_kriging(obs_i, X)

            # push diagnostic outputs
            spush("kriging_field", K)
            spush("kriging_variance", V)

            # execute the krigin update at each model
            Kp = zeros(1)
            Vp = zeros(1,1)
            fuel_types = [2]
            for i in 1:dsize[1]
                for j in 1:dsize[2]
                    Kp[1] = K[i,j]
                    Vp[1,1] = V[i,j]
                    kalman_update(models[i,j], Kp, Vp, fuel_types)
                end
            end

            # move to the next storage frame
            next_frame()
        end
    end

    # Close down the storage system
    Storage.sclose()

end


main(ARGS)
