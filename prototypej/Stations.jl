module stations

using Calendar



type Station
    # station id
    name::String

    # station location (lat/lon decimal)
    ll::Location

    # station elevation (meters)
    elevation::Float64

    # list of observed variable names
    observed_vars::Array{String}

    # variance of observed variables (measurements)
    m_variance::Array{Float64}
end


type Observation
    # station of origin
    station::Station

    # observation time
    tm::Calendar.CalendarTime

    # observed value
    value::Float64

    # name of the variable
    obs_type::String

    # variance of the observation
    var :: Float64

end

function Observation(s::Station, tm::)
        # construct initial state
        m_ext = zeros(2*k+3)
        m_ext[1:k] = m0

        # initialize with defaults
        new(k, m_ext, P0, Tk * 3600, 14.0 * 3600, 5.0, 8.0, 2.5, latlon)
    end

end


end