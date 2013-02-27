#!/usr/bin/env julia

#
#  This is the main file for data assimilation runs.
#
#
# 1. loads up a configuration file,
# 2. obtains data from a WRF model,
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
       Stations.load_station_data, Stations.build_observation_data

using Kriging

using WRF

using FM
import FM.FMModel



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
    setup_tag("assim_K0", false, true, true)
    setup_tag("assim_K1", true, true, true)
    setup_tag("assim_data", false, false, true)

    setup_tag("obs_residual_var", true, true, true)

    setup_tag("fm10_model_residual_var", true, true, true)
    setup_tag("fm10_model_var", false, true, true)
    setup_tag("fm10_kriging_var", false, true, true)

    ### Load WRF output data

    # read in data from the WRF output file pointed to by cfg
    w = WRF.load_wrf_data(cfg["wrf_output"])

    # extract WRF fields
    lat, lon = WRF.lat(w), WRF.lon(w)
    wtm = WRF.times(w)
    dom_shape = size(lat)

    # retrieve equilibria and rain
    Ed, Ew = WRF.field(w, "Ed"), WRF.field(w, "Ew")
    rain = WRF.field(w, "RAIN")

    ### Load observation data from stations
    io = open(join([cfg["station_data_dir"], cfg["station_list_file"]], "/"), "r")
    station_ids = filter(x -> x[1] != '#', map(x -> strip(x), readlines(io)))
    close(io)

    # load each station from its info and observation files
    stations = Station[]
    for sid in station_ids
        s = load_station_info(join([cfg["station_data_dir"], string(sid, ".info")], "/"))
        load_station_data(s, join([cfg["station_data_dir"], string(sid, ".obs")], "/"))
        push!(stations, s)
    end

    # build the observation data from stations
    obs_fm10 = build_observation_data(stations, "FM")

    ### Initialize model

    # maximum equilibrium moisture for visualization
    maxE = cfg["maxE"]

    # construct initial conditions (FIXME: can we do better here?)
    E = squeeze(0.5 * (Ed[1,:,:] + Ew[1,:,:]), 1)

    # set up parameters
    Q = eye(9) * cfg["Q"]
    P0 = eye(9) * cfg["P0"]
    dt = (wtm[2] - wtm[1]).millis / 1000
    println("INFO: time step from WRF is $dt s.")
    K = zeros(Float64, dom_shape)
    V = zeros(Float64, dom_shape)
    mV = zeros(Float64, dom_shape)
    pred = zeros(Float64, dom_shape)
    mresV = zeros(Float64, dom_shape)
    Kfi = zeros(Float64, dom_shape)
    Vfi = zeros(Float64, dom_shape)
    mid = zeros(Int32, dom_shape)
    Kg = zeros(Float64, (dom_shape[1], dom_shape[2], 9))

    # construct model grid from fuel parameters
    Tk = [ 1.0, 10.0, 100.0 ]
    models = [ FMModel((lat[x,y], lon[x,y]), 3, E[x,y], P0, Tk) for x=1:dom_shape[1], y=1:dom_shape[2] ]
    models_na = [ FMModel((lat[x,y], lon[x,y]), 3, E[x,y], P0, Tk) for x=1:dom_shape[1], y=1:dom_shape[2] ]

    println("Done")
    
end


main(ARGS)
