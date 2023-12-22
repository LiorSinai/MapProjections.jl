using Plots
using Images
using ImageTransformations
using LinearAlgebra: norm

using MapProjections
using MapProjections: AffineTransform, affine_from_bounds, bounds

function break_discontinuities!(vertices; threshold=0.5)
    distance_max = largest_diagonal(vertices) * threshold
    v1s = vertices[1:(end-1)]
    v2s = vertices[2:end]
    offset = 1
    for (idx, pair) in enumerate(zip(v1s, v2s))
        d = norm(pair[2] .- pair[1])
        if d > distance_max
            #println("breaking at $idx")
            insert!(vertices, idx + offset, (NaN, NaN))
            offset += 1
        end
    end        
end

function largest_diagonal(vertices::AbstractVector)
    xmin, ymin, xmax, ymax = bounds(vertices)
    norm((xmax - xmin), (ymax - ymin))
end

# Config 
bilinear = ImageTransformations.BSpline(ImageTransformations.Linear())
img_path = "../data/Blue_Marble_2002.png"
out_dir = "../gallery"
long0, lat0 = (29.0, -26.2) # Johannesburg
max_size = 500

## Data
println("Loading image at $img_path")
img = load(img_path)
height, width = size(img)
println("Pixels: $height × $width")
println("")

## WGS-84 
src_proj = WorldGeodeticSystem84()
long_min, long_max = -180.0, 180.0;
lat_min, lat_max = -90.0, 90.0;
height, width = size(img)
src_affine = affine_from_bounds(long_min, lat_min, long_max, lat_max, width, height)

## Azimuthal
println("Azimuthal projection")
azimuthal_proj = AzimuthalEquidistant(;lat0=lat0, long0=long0, k=0.001) # covnert to km

height, width = size(img)
width_az, height_az = 500, 500
r = π * azimuthal_proj.radius
azimuthal_affine = affine_from_bounds(-r, -r, r, r, width_az, height_az);
println("Source:    ($height, $width)")
println("Azimuthal: ($height_az, $width_az)")

azimuthal_img = reproject_warp(
    img, (height_az, width_az), src_proj, azimuthal_proj, src_affine, azimuthal_affine; method=bilinear)
azimthal_img = map(clamp01nan, azimuthal_img)

lines_az = Vector{Vector{Tuple{Float64, Float64}}}()
t= 0.0:0.01:2π
for r in 2500:2500:20_000
    x_az = r .* cos.(t)
    y_az = r .* sin.(t)
    line = collect(zip(x_az, y_az))
    push!(lines_az, line)
end

canvas = plot(
    azimuthal_img, 
    xlims=(1, width_az),
    ylims=(1, height_az),
    size=(500, 500),
    ticks =:none,
)
inv_affine = inv(azimuthal_affine)
for line in lines_az
    coords_pix = map(inv_affine, line)
    x_pix = [xy[1] for xy in coords_pix]
    y_pix = [xy[2] for xy in coords_pix]
    plot!(canvas, x_pix, y_pix, label="", color=:white)
end
coords_centre = (azimuthal_proj.long0, azimuthal_proj.lat0)
coords_az = inv_affine(coords_centre)
scatter!(canvas, [coords_az[1]], [coords_az[2]], color=:white, markershape=:x, label="")

outpath = joinpath(out_dir, "contours_azimuthal.png")
result = savefig(canvas, outpath)
println("Saved Azimuthal plot to $result")
println("")

## Destination
# dest_proj = Orthographic(;k=0.001, long0=long0)
# output_name = "contours_orthographic.png"
dest_proj = Robinson(;k=1.0)
output_name = "contours_robinson.png"

projection_type = typeof(dest_proj).name.wrapper
println("projection_type projection")

height, width = size(img)
## Orthographic
# width_dest, height_dest = 500, 500
# r = dest_proj.radius
# dest_affine = affine_from_bounds(-r, -r, r, r, width_dest, height_dest);
## Robinson
width_dest, height_dest, dest_affine = calculate_suggested_transform(
    src_proj, dest_proj, width, height, src_affine)

println("Source:   ($height, $width)")
println("$projection_type: ($height_dest, $width_dest)")

dest_img = reproject_warp(
    img, (height_dest, width_dest), src_proj, dest_proj, src_affine, dest_affine; method=bilinear)
dest_img = map(clamp01nan, dest_img)

ratio = height_dest / width_dest
figsize = ratio > 1 ? (max_size, round(Int, max_size/ratio)) : (max_size, round(Int, max_size*ratio))
println("Figure size: $figsize")

canvas = plot(
    dest_img,
    xlims=(1, width_dest),
    ylims=(1, height_dest),
    size=figsize,
    ticks =:none,
    )
inv_affine = inv(dest_affine)
for line in lines_az
    coords_dest = reproject(line, azimuthal_proj, dest_proj);
    break_discontinuities!(coords_dest)
    coords_pix = map(inv_affine, coords_dest)
    x_pix = [xy[1] for xy in coords_pix]
    y_pix = [xy[2] for xy in coords_pix]
    plot!(canvas, x_pix, y_pix, label="", color=:white)
end
coords_centre_dest = reproject(coords_centre, azimuthal_proj, dest_proj) |> inv_affine
scatter!(
    canvas, 
    [coords_centre_dest[1]], [coords_centre_dest[2]], 
    color=:white, markershape=:x, label=""
)

outpath = joinpath(out_dir, output_name)
result = savefig(canvas, outpath)
println("Saved $projection_type plot to $result")