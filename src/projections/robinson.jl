"""
    Robinson(radius, long0, interpolatorX, interpolatorY)
    Robinson(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84, interpolator=CubicSpline())

Artistic projection of the entire sphere. 
It is used for visually appealing maps of the whole earth.

The interpolator should either be `CubicSpline()` or `LinearSpline()`. 
`CubicSpline()` gives slightly smoother results near the poles but is about `2x` slower in the forward direction
and about `5x` slower in the inverse direction.
The relative results differ by less than 0.08% below 80° and up to a maximum of 0.16% near 82.2°.

The equations satisfy:
```math
    x = k*r*0.8487*(λ-λ_0)*X(|ϕ|)
    y = sign(ϕ)*k*r*1.3523*Y(|ϕ|)
```
where ``X(ϕ)`` and ``Y(ϕ)`` are interpolatd from a table.

Inputs should be in degrees:
```
proj = Robinson()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (1.81e6, 3.20e6)
```
"""
struct Robinson{T<:AbstractFloat, IX, IY} <: AbstractProjection
    radius::T
    long0::T
    interpolatorX::IX
    interpolatorY::IY
    function Robinson(radius::T, long0::T, interpolator) where {T}
        latitudes = 0.0:5.0:90.0
        X = [
            1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216, 0.8962,
            0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322
        ]
        Y = [
            0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 
            0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000
        ]
        interpolatorX = similar(interpolator, latitudes, X)
        interpolatorY = similar(interpolator, latitudes, Y)
        IX = typeof(interpolatorX)
        IY = typeof(interpolatorY)
        new{T, IX, IY}(radius, long0, interpolatorX, interpolatorY)
    end
end

function inv(proj::Robinson) 
    interpolatorX = proj.interpolatorX
    invInterpolatorY = inv(proj.interpolatorY)
    InverseRobinson(proj.radius, proj.long0, interpolatorX, invInterpolatorY)
end

function Robinson(
    ;long0::AbstractFloat=0.0f0, k::AbstractFloat=1.0f0, radius::AbstractFloat=MEAN_RADIUS_WGS_84,
    interpolator=CubicSpline()
    ) 
    Robinson(promote(radius * k, long0)..., interpolator)
end

function show(io::IO, mime::MIME"text/plain", proj::Robinson)
    print(io, "Robinson(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "interpolatorX=") 
    show(io, proj.interpolatorX)
    print(io, ", interpolatorY=") 
    show(io, proj.interpolatorY)
    print(io, ")")
end

function project(proj::Robinson{T1}, coordinate::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude = bound_longitude(longitude - proj.long0)
    x_interp = interpolate(proj.interpolatorX, abs(latitude))
    y_interp = interpolate(proj.interpolatorY, abs(latitude))
    longitude = degree_to_radian(T1, longitude)
    x = proj.radius * 0.8487 * longitude * x_interp
    y = sign(latitude) * proj.radius * 1.3523 * y_interp
    (x, y)
end

"""
    InverseRobinson(radius, long0, interpolatorX, invInterpolatorY)

Convert `(x, y)` co-ordinates in a `Robinson` projection back to `(longitude, latitude)`. 
"""
struct InverseRobinson{T<:AbstractFloat, IX, IY} <: AbstractProjection
    radius::T
    long0::T
    interpolatorX::IX
    invInterpolatorY::IY
end

function project(proj::InverseRobinson{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    # Not recommended for long0 ≠ 0
    x, y = xy
    x = x / (proj.radius * 0.8487)
    y = y / (proj.radius * 1.3523)
    latitude = sign(y) * proj.invInterpolatorY(abs(y))
    x_interp = interpolate(proj.interpolatorX, abs(latitude))
    longitude = x / x_interp
    longitude = radian_to_degree(T1, longitude)
    longitude += proj.long0
    promote(longitude, latitude)
end