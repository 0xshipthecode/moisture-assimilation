import Base.isequal, Base.isless

type DateTime
    year::Int16
    mon::Int8
    day::Int8
    hour::Int8
    min::Int8
    sec::Int8

    function DateTime(y,m,d,h,mn,s)
        @assert 0 <= y <= 9999
        @assert 1 <= m <= 12
        @assert 1 <= d <= 31
        @assert 0 <= h <= 23
        @assert 0 <= mn <= 59
        @assert 0 <= s <= 59
        new(y,m,d,h,mn,s)
    end

    function DateTime(y,m,d,h,mn,s)
        @assert 0 <= y <= 9999
        @assert 1 <= m <= 12
        @assert 1 <= d <= 31
        @assert 0 <= h <= 23
        @assert 0 <= mn <= 59
        @assert 0 <= s <= 59
        new(y,m,d,h,mn,s)
    end

end

function Date(y,m,d)
    DateTime(y,m,d,0,0,0)
end

function Time(h,mn,s)
    DateTime(1,1,1,h,mn,s)
end



function parse_wrf_date(str::String)
    tok = map(x -> int(x), split(str, ['-', '_', ':']))
    DateTime(tok[1], tok[2], tok[3], tok[4], tok[5], tok[6])
end


isequal(x::DateTime, y::DateTime) =
    x.year == y.year && x.mon == y.mon && x.day == y.day && x.hour == y.hour && x.min == y.min && x.sec == y.sec


function isless(x::DateTime, y::DateTime)
    """
    Standard comparison for sort() and general '<' use.
    """
    (x.year < y.year) && return true
    (x.year > y.year) && return false
    (x.mon < y.mon) && return true
    (x.mon > y.mon) && return false
    (x.day < y.day) && return true
    (x.day > y.day) && return false
    (x.hour < y.hour) && return true
    (x.hour > y.hour) && return false
    (x.min < y.min) && return true
    (x.min > y.min) && return false
    (x.sec < y.sec) && return true
    return false
end

