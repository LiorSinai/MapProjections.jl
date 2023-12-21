module MapProjections

using StaticArrays
using ImageTransformations
import Base: inv, show, similar

export AbstractProjection, earth_radius_at_latitude
export AzimuthalEquidistant, CylindricalEqualArea, Mercator, 
        Orthographic, Robinson, TransverseMercator, EllipsoidalTransverseMercator,
        WorldGeodeticSystem84
export project
export reproject, reproject_warp, calculate_suggested_transform

include("Interpolation.jl")
using .Interpolation


include("affine.jl")
include("bounds.jl")
include("projection.jl")
include("recentre.jl")
include("reproject.jl")

end 
