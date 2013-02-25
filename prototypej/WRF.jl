module WRF

#
# Module for reading in WRF outputs and re-packaging selected variables
# for use with the moisture simulation/data assimilation code.
#
#

using Calendar
import Calendar.CalendarTime

using netcdf


type WRFData

    # source file path
    file_path :: String

    # the variable fields
    fields :: Dict{String,Array{Float64}}

    # the times of recording history
    tm :: Array{CalendarTime}

    WRFData(file_path) = new(file_path, Dict{String,Array{Float64}}(), Array(CalendarTime,0))

end



function load_wrf_data(file_path::String)

    w = WRFData(file_path)
    w.file_path = file_path

    nc_file = netcdf.open(file_path)

    for v in ["T2", "Q2", "PSFC", "RAINC", "RAINNC"]
        w.fields[v] = netcdf.ncread(file_path, v)
    end

    # time is stored as a sequence of 19-byte datetime markers,
    # however the orientation of the read-in variable is incorrect
    # size is T x 19 but should be 19 x T, thus the reshape.
    tm = netcdf.ncread(file_path, "Times")
    s = size(tm)
    tm = reshape(tm, (s[2], s[1]))
    for i in 1:size(tm,2)
        push!(w.tm, Calendar.parse("yyyy-MM-dd_HH:mm:ss", ascii(squeeze(tm[:,i], 2)), "GMT"))
    end

    # read in grid dimensions
    w.fields["lon"] = netcdf.ncread(file_path, "XLONG")
    w.fields["lat"] = netcdf.ncread(file_path, "XLAT")

    # compute derived variables
    compute_rainfall_per_hour(w)
    compute_equilibrium_moisture(w)

    netcdf.close(nc_file)

    return w
end


function lat(w::WRFData)
    return w.fields["lat"]
end

function lon(w::WRFData)
    return w.fields["lon"]
end

function times(w::WRFData)
    return w.tm
end

function field(w::WRFData, fname)
    return w.fields[fname]
end


function compute_rainfall_per_hour(w::WRFData)

    rainnc = field(w, "RAINNC")
    rainc = field(w, "RAINC")
    rain = zeros(Float64, size(rainnc))
    rain_old = zeros(Float64, size(rain[1,:,:]))

    for i in 2:length(w.tm)
        dt = (w.tm[i] - w.tm[i-1]).millis / 1000
        rain[i,:,:] = ((rainc[i,:,:] + rainnc[i,:,:]) - rain_old) * 3600.0 / dt
        rain_old[:,:,:] = rainc[i,:,:]
        rain_old += rainnc[i,:,:]
     end

     w.fields["RAIN"] = rain
     delete!(w.fields, "RAINC")
     delete!(w.fields, "RAINNC")
end


function domain_extent(w::WRFData)
    lon = field(w, "LON")
    lat = field(w, "LAT")

    return (min(lon), min(lat), max(lon), max(lat))
end


function compute_equilibrium_moisture(w::WRFData)
    P, Q, T = field(w, "PSFC"), field(w, "Q2"), field(w, "T2")

    N = size(P,1)

    # Xi stands for X interpolated
    Pi = 0.5 * (P[1:N-1,:,:] + P[2:N,:,:])
    Qi = 0.5 * (Q[1:N-1,:,:] + Q[2:N,:,:])
    Ti = 0.5 * (T[1:N-1,:,:] + T[2:N,:,:])

    # compute saturated vapor pressure
    Pws = exp(54.842763 - 6763.22/Ti - 4.210 * log(Ti) + 0.000367*Ti
              + tanh(0.0415*(Ti - 218.8)) .*
              (53.878 - 1331.22/Ti - 9.44523 * log(Ti) + 0.014025*Ti))

    # compute current water vapor pressure
    Pw = Pi .* Qi ./ (0.622 + (1 - 0.622) * Qi)

    # compute relative humidity in percent
    RH = 100 * Pw ./ Pws

    # drying/wetting fuel equilibrium moisture contents
    Ed = zeros(Float64, size(P))
    Ew = zeros(Float64, size(P))

    Ed[2:,:,:] = 0.924*RH.^0.679 + 0.000499*exp(0.1*RH) + 0.18*(21.1 + 273.15 - Ti).*(1 - exp(-0.115*RH))
    Ew[2:,:,:] = 0.618*RH.^0.753 + 0.000454*exp(0.1*RH) + 0.18*(21.1 + 273.15 - Ti).*(1 - exp(-0.115*RH))

    Ed *= 0.01
    Ew *= 0.01

    w.fields["Ed"] = Ed
    w.fields["Ew"] = Ew

end


end