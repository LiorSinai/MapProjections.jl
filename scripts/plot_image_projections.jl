using Plots
using Images
using ImageTransformations
using ImageTransformations: BSpline, Linear
using MapProjections
using MapProjections: affine_from_bounds

## Config
bilinear = ImageTransformations.BSpline(ImageTransformations.Linear())
output_dir = "../gallery/blue_marble/"
img_path = "../data/Blue_Marble_2002.png"
long_min, long_max = -180.0, 180.0;
lat_min, lat_max = -90.0, 90.0;
alpha = 0.5 # transparency for grid lines
max_figure_size = 500

## Functions
include("grid.jl")

## data
print("Loading image at $img_path ... ")
img = load(img_path)
println("done")
height, width = size(img)
println("size: $height × $width")
println("")

## Source
src_proj = WorldGeodeticSystem84()
height, width = size(img)
src_affine = affine_from_bounds(long_min, lat_min, long_max, lat_max, width, height)

## Projections
include("image_projections/azimuthal_equal_distant.jl")
include("image_projections/cylindrical_equal_area.jl")
include("image_projections/equirectangular.jl")
include("image_projections/mercator.jl")
include("image_projections/orthographic.jl")
include("image_projections/robinson.jl")
include("image_projections/sinusoidal.jl")
include("image_projections/transverse_mercator.jl")
include("image_projections/transverse_mercator_ellipsoidal.jl")
