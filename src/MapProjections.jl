module MapProjections

using StaticArrays
using ImageTransformations
import Base: inv, show, similar

export AbstractProjection, earth_radius_at_latitude
export AzimuthalEquidistant, CylindricalEqualArea, Mercator, 
        Orthographic, Robinson, TransverseMercator, EllipsoidalTransverseMercator,
        Sinusoidal, WorldGeodeticSystem84
export project
export reproject, reproject_warp, calculate_suggested_transform

include("Interpolation.jl")
using .Interpolation

include("affine.jl")
include("bounds.jl")
include("projection.jl")
include("recentre.jl")
include("reproject.jl")

include("projections\\WGS_84.jl")
include("projections\\azimuthal.jl")
include("projections\\cylindrical_equal_area.jl")
include("projections\\mercator.jl")
include("projections\\sinusoidal.jl")
include("projections\\robinson.jl")
include("projections\\transverse_mercator.jl")

end 
