using MapProjections: bound_longitude

println("EllipsoidalTransverseMercator")
dest_proj = EllipsoidalTransverseMercator(;long0=0.0, semi_major_axis=1.0)

height, width = size(img)
west_boundary = -45.0 + dest_proj.long0
east_boundary = 45.0 + dest_proj.long0
j, i = inv(src_affine)((west_boundary, 0.0))
j_west = floor(Int, j)
j, i = inv(src_affine)((east_boundary, 0.0))
j_east = floor(Int, j)
if j_west < 0 
    js = vcat((width + j_west):width, 1:j_east)
    j_west += width
elseif j_east > width
    js = vcat(j_west:width, 1:(j_east - width))
    j_east -= width
else
    js = j_west:j_east
end
println("cropped j: ", j_west, "-", j_east)

img_cropped = img[1:end, js];

height_crop, width_crop = size(img_cropped)
crop_affine = affine_from_bounds(west_boundary, lat_min, east_boundary, lat_max, width_crop, height_crop)
width_dest, height_dest, dest_affine = 
    calculate_suggested_transform(src_proj, dest_proj, width_crop, height_crop, crop_affine)
println("Source:      ($height, $width)")
println("Destination: ($height_dest, $width_dest)")

out_img = reproject_warp(
    img_cropped, (height_dest, width_dest), src_proj, dest_proj, crop_affine, dest_affine; method=bilinear)

ratio = width_dest / height_dest
figure_size = (ceil(Int, max_figure_size * ratio), ceil(Int, max_figure_size * ratio))
println("figure_size: $figure_size")

canvas = plot(
    out_img, 
    xlims=(1, width_dest),
    ylims=(1, height_dest),
    size=figure_size, 
    ticks=:none,
    )

gridlines = make_grid_lines(
    dest_affine, dest_proj, west_boundary:15.0:east_boundary, -90:30.0:90; num_points=1000
    );
for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:white, label="", alpha=alpha)
end

output_path = joinpath(output_dir, "transverse_mercator_ellipsoidal.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")