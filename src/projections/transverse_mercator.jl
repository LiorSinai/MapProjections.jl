"""
    TransverseMercator(radius, long0)
    TransverseMercator(;long0=0.0f0, k=1.0f0, radius=MEAN_RADIUS_WGS_84)

Projection from a sphere  which corrects distortion along the longitude `long0` at the expense of large distortion away from it.

It is used to correct local distortion in small regions within 15° of `long0`.

The equations satisfy:
```math
    t = sin(λ-λ_0)cos(ϕ)
    x = k*r*0.5*log((1+t) / (1-t))
    y = k*r*atan(tan(ϕ)*sec(λ-λ_0))
```

Inputs should be in degrees:
```
proj = TransverseMercator()
coords_src = (20.0f0, 30.0f0) # (longitude, latitude) degrees
coords_dest = proj(coords_src) # (1.95e6, 3.51e6) #m
```
"""
struct TransverseMercator{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
end

function show(io::IO, mime::MIME"text/plain", proj::TransverseMercator)
    print(io, typeof(proj), "(")
    print(io, "radius=", proj.radius, ", ")
    print(io, "long0=", proj.long0)
    print(io, ")")
end

function TransverseMercator(;long0::AbstractFloat=0.0f0, k::AbstractFloat=1.0f0, radius::AbstractFloat=MEAN_RADIUS_WGS_84) 
    TransverseMercator(promote(radius * k, long0)...)
end

inv(proj::TransverseMercator) = InverseTransverseMercator(proj.radius, proj.long0)

function project(proj::TransverseMercator{T1}, coordinate::Tuple{T2, T2}; extend::Bool=false) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    longitude = bound_longitude(longitude - proj.long0)
    longitude = degree_to_radian(T1, longitude)
    latitude  = degree_to_radian(T1, latitude)
    t = sin(longitude) * cos(latitude)
    x = proj.radius * convert(T1, 0.5) * log((one(T1) + t) / (one(T1) - t))
    y = proj.radius * atan(tan(latitude) * sec(longitude))
    if extend && (abs(longitude) > π/2) # rear of sphere
        y += proj.radius * π
    end
    (x, y)
end

"""
    InverseTransverseMercator(radius, long0)

Convert `(x, y)` co-ordinates in a `SphericalTransverseMercator` projection back to `(longitude, latitude)`. 
"""
struct InverseTransverseMercator{T<:AbstractFloat} <: AbstractProjection
    radius::T
    long0::T
end

inv(proj::InverseTransverseMercator) = TransverseMercator(proj.radius, proj.long0)

function project(proj::InverseTransverseMercator{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    x, y = xy
    x = x / proj.radius
    y = y / proj.radius
    longitude = atan(sinh(x) * sec(y))
    latitude = asin(sech(x) * sin(y))
    longitude = radian_to_degree(T1, longitude)
    latitude = radian_to_degree(T1, latitude)
    longitude += proj.long0
    if abs(y) > π/2
        if longitude < 0
            longitude += 180
        else
            longitude -= 180
        end
    end
    (longitude, latitude)
end

"""
    EllipsoidalTransverseMercator(;
        semi_major_axis=SEMI_MAJOR_AXIS_WGS_84,
        flattening=1/INV_FLATTENING_WGS_84
        long0=0.0f0, 
        k=1.0f0, 
        )

A version of the transverse mercator which maps from an ellipsoid to a cylinder.
It is more accurate than `TransverseMercator` for removing distortion in small local regions.

This implementation uses the Gauss-Krüger equations which approximate the exact soluton.
Because of numerical instability it is not recommended for more than 45° on either side of `long0`.
While the exact solution has a finite map, these approximations will tend to infinity at the points at `long0±90°`.

Source:
- https://proj.org/en/9.3/operations/projections/tmerc.html
- https://pro.arcgis.com/en/pro-app/latest/help/mapping/properties/gauss-kruger.htm
"""
struct EllipsoidalTransverseMercator{T<:AbstractFloat, TM<:TransverseMercator} <: AbstractProjection
    semi_major_axis::T
    flattening::T
    long0::T
    spherical_crs::TM
    ## constants
    n::T # third flattening
    e::T # eccentricity
    A::T # 2πA is the circumfernce of a meridian
    C::Matrix{T}
    N::Matrix{T}
    function EllipsoidalTransverseMercator(a::T, f::T , long0::T) where T<:AbstractFloat
        n = f / (2 - f)
        e = sqrt(2f - f*f)
        n2 = n * n
        n4 = n2 * n2
        n6 = n2 * n4
        N = reshape([n; n2; n2*n; n4; n4*n; n6], :, 1)
        C = [
            1/2 -2/3   5/16   41/180       -127/288      7891/37800;
            0   13/48 -3/5    557/1400      281/630     -1983433/1935360;
            0   0     61/240 -103/140       15061/26880  167603/181440;
            0   0     0       49561/161280 -179/168      6601661/7257600;
            0   0     0       0             34729/80640 -3418889/1995840;
            0   0     0       0             0            212378941/319334400 ; 
        ]
        spherical_crs = TransverseMercator(;radius=one(T), long0=long0)
        A = a / (1 + n) * (1 + n2/4 + n4/64 + n6/256)
        new{T, TransverseMercator{T}}(a, f, long0, spherical_crs, n, e, A, C, N)
    end
end

inv(proj::EllipsoidalTransverseMercator) = 
    InverseEllipsoidalTransverseMercator(proj.semi_major_axis, proj.flattening, proj.long0)

function show(io::IO, mime::MIME"text/plain", proj::EllipsoidalTransverseMercator)
    print(io, typeof(proj), "(")
    print(io, "semi_major_axis=", proj.semi_major_axis, ", ")
    print(io, "flattening=", proj.flattening, ", ")
    print(io, "long0=", proj.long0, ", ")
    print(io, "spherical_crs=", proj.spherical_crs, ", ")
    print(io, "n=", proj.n, ", ")
    print(io, "e=", proj.e, ", ")
    print(io, "A=", proj.A, ", ")
    print(io, "C=", proj.C, ", ")
    print(io, "N=", proj.N)
    print(io, ")")
end

function EllipsoidalTransverseMercator(
    ;long0::AbstractFloat=0.0f0, 
    k::AbstractFloat=1.0f0,
    semi_major_axis::AbstractFloat=SEMI_MAJOR_AXIS_WGS_84,
    flattening::AbstractFloat=1/INV_FLATTENING_WGS_84
    )
    EllipsoidalTransverseMercator(promote(semi_major_axis * k, flattening, long0)...)
end

function project(proj::EllipsoidalTransverseMercator{T1}, coordinate::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    longitude, latitude = coordinate
    latitude  = degree_to_radian(T1, latitude)
    e = sqrt(4*proj.n)/(1+proj.n)
    latitude_conformal = atan(sinh(asinh(tan(latitude)) - e*atanh(e*sin(latitude))))
    latitude_conformal = radian_to_degree(latitude_conformal)
    x, y = proj.spherical_crs((longitude, latitude_conformal))
    z_sphere = Complex(y, x)
    S = Sfunc(z_sphere) 
    z = proj.A * (z_sphere + first(S * proj.C * proj.N))
    (imag(z), real(z))
end

Sfunc(z::Number) = [
    sin(2 * z) sin(4 * z) sin(6 * z) sin(8 * z) sin(10 * z) sin(12 * z);
]

function _trig_complex_addition(s1::T, c1::T, sh1::T, ch1::T, s2::T, c2::T, sh2::T, ch2::T) where T<:Real
    s3 = s1*c2 + c1*s2
    c3 = c1*c2 - s1*s2
    sh3 = sh1*ch2 + ch1*sh2
    ch3 = ch1*ch2 + sh1*sh2
    (s3, c3, sh3, ch3)
end

function Sfunc_identities(z::Complex)
    s = sin(real(z))
    c = cos(real(z))
    sh = sinh(imag(z))
    ch = cosh(imag(z))
    # sin(2z)
    s2, c2, sh2, ch2 = _trig_complex_addition(s, c, sh, ch, s, c, sh, ch)
    sin2z = s2 * ch2 + c2 * sh2 * im
    # sin(4z)
    s4, c4, sh4, ch4 = _trig_complex_addition(s2, c2, sh2, ch2, s2, c2, sh2, ch2)
    sin4z = s4 * ch4 + c4 * sh4 * im
    # sin(6z)
    s6, c6, sh6, ch6 = _trig_complex_addition(s4, c4, sh4, ch4, s2, c2, sh2, ch2)
    sin6z = s6 * ch6 + c6 * sh6 * im
    # sin(8z)
    s8, c8, sh8, ch8 = _trig_complex_addition(s4, c4, sh4, ch4, s4, c4, sh4, ch4)
    sin8z = s8 * ch8 + c8 * sh8 * im
    # sin(10z)
    s10, c10, sh10, ch10 = _trig_complex_addition(s8, c8, sh8, ch8, s2, c2, sh2, ch2)
    sin10z = s10 * ch10 + c10 * sh10 * im
    # sin(12z)
    s12, c12, sh12, ch12 = _trig_complex_addition(s6, c6, sh6, ch6, s6, c6, sh6, ch6)
    sin12z = s12 * ch12 + c12 * sh12 * im
    [
        sin2z sin4z sin6z sin8z sin10z sin12z
    ]
end

function Sfunc_identities(t::Real)
    s = sin(t)
    c = cos(t)
    sin2t = 2 * s * c
    cos2t = c * c - s * s
    sin4t = 2 * sin2t * cos2t
    cos4t = cos2t * cos2t - sin2t * sin2t
    sin6t = sin4t * cos2t + cos4t * sin2t
    cos6t = cos4t * cos2t - sin4t * sin2t
    sin8t = 2 * sin4t * cos4t
    sin10t = sin6t * cos4t + cos6t * sin4t
    sin12t = 2 *sin6t * cos6t
    [
        sin2t sin4t sin6t sin8t sin10t sin12t
    ]
end

"""
    InverseEllipsoidalTransverseMercator(radius, long0)

Convert `(x, y)` co-ordinates in a `TransverseMercator` projection back to `(longitude, latitude)`. 
"""
struct InverseEllipsoidalTransverseMercator{T<:AbstractFloat, ITM<:InverseTransverseMercator} <: AbstractProjection
    semi_major_axis::T
    flattening::T
    long0::T
    spherical_crs::ITM
    ## constants
    n::T # third flattening
    e::T # eccentricity
    A::T # 2πA is the circumfernce of a meridian
    Cϕχ::Matrix{T}
    Cχμ::Matrix{T}
    N::Matrix{T}
    function InverseEllipsoidalTransverseMercator(a::T, f::T , long0::T) where T<:AbstractFloat
        n = f / (2 - f)
        e = sqrt(2f - f*f)
        n2 = n * n
        n4 = n2 * n2
        n6 = n2 * n4
        N = reshape([n; n2; n2*n; n4; n4*n; n6], :, 1)
        Cϕχ = [
            2 -2/3 -2      116/45    26/45    -2854/675 ;
            0  7/3 -8/5   -227/45    2704/315  2323/945 ;
            0  0    56/15 -136/35    1262/105  73814/2835;
            0  0    0      4279/630 -332/35   -399572/14175 ;
            0  0    0      0         4174/315  144838/6237 ;
            0  0    0      0         0         601676/22275;
        ]
        Cχμ = [
            -1/2  2/3  -37/96   1/360        81/512      -96199/604800 ;
            0    -1/48 -1/15    437/1440    -46/105       1118711/3870720;
            0     0    -17/480  37/840       209/4480    -5569/90720;
            0     0     0      -4397/161280  11/504       830251/7257600;
            0     0     0       0           -4583/161280  108847/3991680;
            0     0     0       0            0           -20648693/638668800 ; 
        ]
        spherical_crs = InverseTransverseMercator(1.0, long0)
        A = a / (1 + n) * (1 + n2/4 + n4/64 + n6/256)
        new{T, InverseTransverseMercator{T}}(a, f, long0, spherical_crs, n, e, A, Cϕχ, Cχμ, N)
    end
end

function project(proj::InverseEllipsoidalTransverseMercator{T1}, xy::Tuple{T2, T2}) where {T1,T2 <: AbstractFloat}
    z = xy[2] / proj.A  + xy[1] / proj.A * im
    z_sphere = z + first(Sfunc(z) * proj.Cχμ * proj.N)
    coord_sphere = (imag(z_sphere), real(z_sphere))
    longitude, latitude_conformal = proj.spherical_crs(coord_sphere)
    latitude_conformal = degree_to_radian(T1, latitude_conformal)
    latitude = latitude_conformal + first(Sfunc(latitude_conformal) * proj.Cϕχ * proj.N)
    latitude = radian_to_degree(T1, latitude)
    #longitude += proj.long0
    (longitude, latitude)
end
