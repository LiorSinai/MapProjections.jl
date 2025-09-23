function angular_distance(
    lat0::AbstractFloat, 
    longitude::AbstractFloat, 
    latitude::AbstractFloat;
    )
    acos(sin(lat0) * sin(latitude) + cos(lat0) * cos(latitude) * cos(longitude))
end

function in_angular_range(
    lat0::AbstractFloat, 
    longitude::AbstractFloat, 
    latitude::AbstractFloat;
    atol::AbstractFloat = 1e-4
    )
    c = acos(sin(lat0) * sin(latitude) + cos(lat0) * cos(latitude) * cos(longitude))
    abs(c) < π/2 - atol
end

function _project_azimuthal(
    radius::T, lat0::T, longitude::T, latitude::T;
    ) where {T <: AbstractFloat}
    x = radius * cos(latitude) * sin(longitude)
    y = radius * (cos(lat0) * sin(latitude) - sin(lat0) * cos(latitude) * cos(longitude))
    (x, y)
end

function _inv_azimuthal(
    radius::T, lat0::T, x::T, y::T, max_radius::Number, angular_distance_func
    ) where {T <: AbstractFloat}
    x = x / (radius)
    y = y / (radius)
    ρ = sqrt(x * x + y * y)
    if ρ > max_radius
        nan = convert(T, NaN)
        return (nan, nan)
    end
    c = angular_distance_func(ρ)
    longitude = atan(x * sin(c), (ρ * cos(c) * cos(lat0) - y * sin(c) * sin(lat0)))
    latitude = asin(cos(c) * sin(lat0) + y * sin(c) * cos(lat0) / ρ )
    (longitude, latitude)
end

"""
    Orthographic(radius, long0, lat0)
    Orthographic(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

The "view from space". An orthographic perspective of the earth, resulting in a disk with the illusion of a 3D sphere.

The equations satisfy:
```math
x = k*r*cos(ϕ) * sin(λ-λ0)
y = k*r[cos(ϕ0) * sin(ϕ) - sin(ϕ0) * cos(ϕ) * cos(λ-λ0)]
```

All points on the front disk will have an angular distance between -π/2 and +π/2:
```math
cos(c) = sin(ϕ0) * sin(ϕ) + cos(ϕ0) * cos(ϕ) * cos(λ-λ0)
```

Inputs should be in degrees:
```
proj = Orthographic()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (1.89e6, 3.19e6)
```

References
- https://en.wikipedia.org/wiki/Orthographic_map_projection
- https://mathworld.wolfram.com/OrthographicProjection.html 
"""
struct Orthographic{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    lat0::T
end

function Orthographic(;
    long0::AbstractFloat=0.0f0, lat0::AbstractFloat=0.0f0,
    k::AbstractFloat=1.0f0,
    radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    Orthographic(promote(k * radius, long0, lat0)...)
end

inv(proj::Orthographic) = InverseOrthographic(proj.radius, proj.long0, proj.lat0)

"""
    project(proj::Orthographic, coordinate; wrap=false, clip=true, atol=1e-4)

Project coordinate onto an Orthographic projection.

Options:
- `wrap`: if true, return all coordinates including those out of view. Takes precedence.
- `clip`: if true return `(NaN, NaN)` for all coordinates out of view.
- `atol`: tolerance for `in_angular_range`.

If `wrap` is false and so is `clip`, project coordinates to the edge for a smooth appearance. 
"""
function project(
        proj::Orthographic{T1}, coordinate::Tuple{T2, T2};
        wrap::Bool=false, clip::Bool=true, atol=1e-4
        ) where {T1,T2 <: AbstractFloat}
        longitude, latitude = coordinate
    longitude, latitude = coordinate
    longitude = convert(T1, deg2rad(longitude - proj.long0))
    latitude  = convert(T1, deg2rad(latitude))
    lat0 = convert(T1, deg2rad(proj.lat0))
    x, y = _project_azimuthal(proj.radius, lat0, longitude, latitude)
    if !wrap && !(in_angular_range(lat0, longitude, latitude; atol=atol))
        if clip
            nan = convert(T1, NaN)
            return (nan, nan)
        else
            # project to edge for smooth curved boundary.
            # warning: objects around the polar opposite (c=π) will be projected to the entire circle.
            # it is suggested that these are objects are removed entirely.
            angle = atan(y, x)
            x = proj.radius * cos(angle)
            y = proj.radius * sin(angle)
        end
    end
    (x, y)
end

function show(io::IO, mime::MIME"text/plain", proj::Orthographic)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "lat0=", proj.lat0)
    print(io, ")")
end

"""
    InverseOrthographic(radius, long0, lat0)

Convert `(x, y)` co-ordinates in a `Orthographic` projection back to `(longitude, latitude)`. 
"""
struct InverseOrthographic{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    lat0::T
end

inv(proj::InverseOrthographic) = Orthographic(proj.radius, proj.long0, proj.lat0)

function project(proj::InverseOrthographic{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    lat0 = deg2rad(proj.lat0)
    longitude, latitude = _inv_azimuthal(proj.radius, lat0, xy[1], xy[2], one(T1), asin)
    longitude = rad2deg(longitude)
    latitude = rad2deg(latitude)
    longitude += proj.long0
    if longitude < -180
        longitude += 360
    elseif  longitude > 180
        longitude -= 360
    end
    (longitude, latitude)
end

"""
    AzimuthalEquidistant(radius, long0, lat0)
    AzimuthalEquidistant(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

An azimuthal projection.

The equations satisfy:
```math
x = k*c/sin(c)*r*cos(ϕ) * sin(λ-λ0)
y = k*c/sin(c)*r*[cos(ϕ0) * sin(ϕ) - sin(ϕ0) * cos(ϕ) * cos(λ-λ0)]
cos(c) = sin(ϕ0) * sin(ϕ) + cos(ϕ0) * cos(ϕ) * cos(λ-λ0)
```

Inputs should be in degrees:
```
proj = AzimuthalEquidistant()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (2.01f6, 3.40f6)
```

References
- https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
- https://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
"""
struct AzimuthalEquidistant{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    lat0::T
end

function AzimuthalEquidistant(;
    long0::AbstractFloat=0.0f0, lat0::AbstractFloat=0.0f0,
    k::AbstractFloat=1.0f0,
    radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    AzimuthalEquidistant(promote(k * radius, long0, lat0)...)
end

inv(proj::AzimuthalEquidistant) = InverseAzimuthalEquidistant(proj.radius, proj.long0, proj.lat0)

function project(
    proj::AzimuthalEquidistant{T1}, coordinate::Tuple{T2, T2};
    ) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude, latitude = coordinate
    longitude = deg2rad(longitude - proj.long0)
    latitude  = deg2rad(latitude)
    lat0 = deg2rad(proj.lat0)
    if lat0 ≈ π/2
        # Antarctica is not rendered properly for 84°-90°.
        # This provides a fix for lat0=90° and is also more computationally efficient.
        radius = proj.radius * (π/2 - latitude)
        x = radius * sin(longitude) 
        y = -radius * cos(longitude)
        return (x, y)
    end
    c = angular_distance(lat0, longitude, latitude)
    radius = proj.radius * c / sin(c)
    _project_azimuthal(radius, lat0, longitude, latitude)
end

function show(io::IO, mime::MIME"text/plain", proj::AzimuthalEquidistant)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "lat0=", proj.lat0)
    print(io, ")")
end

"""
    InverseOrthographic(radius, long0, lat0)

Convert `(x, y)` co-ordinates in a `Orthographic` projection back to `(longitude, latitude)`. 
"""
struct InverseAzimuthalEquidistant{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
    lat0::T
end

inv(proj::InverseAzimuthalEquidistant) = AzimuthalEquidistant(proj.radius, proj.long0, proj.lat0)

function project(proj::InverseAzimuthalEquidistant{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    lat0 = deg2rad(proj.lat0)
    longitude, latitude = _inv_azimuthal(proj.radius, lat0, xy[1], xy[2], π, identity)
    longitude = rad2deg(longitude)
    latitude = rad2deg(latitude)
    longitude += proj.long0
    if longitude < -180
        longitude += 360
    elseif  longitude > 180
        longitude -= 360
    end
    (longitude, latitude)
end
