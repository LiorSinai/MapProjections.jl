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

### constants
SEMI_MAJOR_AXIS_WGS_84 = 6_378_137.0f0
INV_FLATTENING_WGS_84 = 298.257223563f0
SEMI_MINOR_AXIS_WGS_84 = SEMI_MAJOR_AXIS_WGS_84 * (1 - 1/INV_FLATTENING_WGS_84)
## The oblate ellipsoid has 2 of the same major semi-axis (x and y radii) and a slightly flatter minor axis (z radius).
MEAN_RADIUS_WGS_84 = (2 * SEMI_MAJOR_AXIS_WGS_84 + SEMI_MINOR_AXIS_WGS_84) / 3
AUTHALIC_RADIUS_WGS_84 = 6_371_007.2f0 # radius of sphere with same surface area as ellipsoid

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
