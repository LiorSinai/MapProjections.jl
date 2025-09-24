"""
    AbstractProjection

An abstract projection. Subtypes should implement the following methods:

```
project(proj, coordinate::Tuple; options...) → Tuple
inv(proj) → AbstractProjection
```
The `coordinate` is a tuple of `(longitude, latitude)` and the outputs are `(x, y)` pairs.
"""
abstract type AbstractProjection end

function (proj::AbstractProjection)(xy::Tuple{T, T}; options...) where {T <: AbstractFloat}
    project(proj, xy; options...)
end

(proj::AbstractProjection)(xy::Tuple{T1, T2}) where {T1, T2} = proj(promote(xy...))

function (proj::AbstractProjection)(xy::SVector{2, T}; options...) where {T <: AbstractFloat}
    project(proj, (xy[1], xy[2]); options...)
end

"""
    project(projection, coordinate)

Project the coordinate onto a new set of coordinates.
"""
function project end

function bound_longitude(x::AbstractFloat)
    if (x < -180)
        return x + 360
    elseif  x > 180
        return x - 360
    end
    x
end

"""
    earth_radius_at_latitude(latitude=35.196f0)

Earth radius at latitude (in degrees) based on the WGS-84 ellipsoid.
The default is set so that it is close to the mean radius: 6,371,008.5 m.

The radius ``r`` and latitude ``ϕ`` satisfy the following equations:
```math
 (x/a)^2 + (y/b)^2 = 1 
    x = r*cos(ϕ)
    y = r*sin(ϕ)
```
"""
function earth_radius_at_latitude(latitude::AbstractFloat=35.196f0)
    a = SEMI_MAJOR_AXIS_WGS_84
    b = SEMI_MINOR_AXIS_WGS_84
    latitude = latitude * π / 180
    radius = a * b / sqrt(b^2 * cos(latitude)^2 + a^2 * sin(latitude)^2)
    radius
end
