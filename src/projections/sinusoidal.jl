"""
    Sinusoidal(radius, long0)
    Sinusoidal(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

Projection from a sphere which preserves area.

It is used for calculating areas.

The equations satisfy:
```math
    x = k*r*(λ-λ_0)cos(ϕ)
    y = k*r*ϕ
```

Inputs should be in degrees:
```
proj = Sinusoidal()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (1.93f6, 3.34f6)
```
"""
struct Sinusoidal{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
end

inv(proj::Sinusoidal) = InverseSinusoidal(proj.radius, proj.long0)

function show(io::IO, mime::MIME"text/plain", proj::Sinusoidal)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0)
    print(io, ")")
end

function Sinusoidal(;long0::AbstractFloat=0.0f0, k::AbstractFloat=1.0f0, radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    Sinusoidal(promote(radius * k, long0)...)
end

function project(proj::Sinusoidal{T1}, coordinate::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude -= proj.long0
    longitude = degree_to_radian(T1, longitude)
    latitude  = degree_to_radian(T1, latitude)
    x = proj.radius * longitude * cos(latitude)
    y = proj.radius * latitude
    (x, y)
end

"""
    InverseSinusoidal(radius, long0)

Convert `(x, y)` co-ordinates in a `Sinusoidal` projection to `(longitude, latitude)`. 
"""
struct InverseSinusoidal{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
end

inv(proj::InverseSinusoidal) = Sinusoidal(proj.radius, proj.long0)

function show(io::IO, mime::MIME"text/plain", proj::InverseSinusoidal)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0)
    print(io, ")")
end

function project(proj::InverseSinusoidal{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    x, y = xy
    latitude = y / proj.radius
    longitude = x / (proj.radius * cos(latitude))
    longitude = radian_to_degree(T1, longitude)
    latitude = radian_to_degree(T1, latitude)
    longitude += proj.long0
    (longitude, latitude)
end