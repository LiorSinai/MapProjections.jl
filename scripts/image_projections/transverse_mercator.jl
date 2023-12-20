println("TransverseMercator")

dest_proj = TransverseMercator(;k=1.0, long0=0.0)

height, width = size(img)
west_boundary = -90.0 + dest_proj.long0
east_boundary = 90.0 + dest_proj.long0
j, i = inv(src_affine)((west_boundary, 0.0))
j_west = floor(Int, j)
j, i = inv(src_affine)((east_boundary, 0.0))
j_east = floor(Int, j)
if j_west < 0 
    js = vcat((width + j_west):width, 1:(j_east))
    js_back = (j_east - 1):(width + j_west + 1) # slight overlap
    j_west += width
elseif j_east > width
    js = vcat(j_west:width, 1:(j_east - width))
    js_back = (j_east - width - 1):(j_west + 1) # slight overlap
    j_east -= width
else
    js = j_west:j_east
    js_back = vcat((j_east-1):width, 1:(j_west+1)) # slight overlap
end
println("cropped j: ", j_west, "-", j_east)
println(length(js), " ", length(js_back))

img_cropped = img[1:end, js];
img_back = img[1:end, js_back]
img_front = img_cropped; # front of sphere

height_crop, width_crop = size(img_cropped)
crop_affine = affine_from_bounds(west_boundary, lat_min, east_boundary, lat_max, width_crop, height_crop)
width_dest, height_dest, dest_affine = 
    calculate_suggested_transform(src_proj, dest_proj, width_crop, height_crop, crop_affine)
println("Source:      ($height, $width)")
println("Destination: ($height_dest, $width_dest)")

out_img_front = reproject_warp(
    img_front, (height_dest, width_dest), src_proj, dest_proj, crop_affine, dest_affine; method=bilinear)
out_img_back = reproject_warp(
    img_back, (height_dest, width_dest), src_proj, dest_proj, crop_affine, dest_affine; method=bilinear)
out_img = vcat(out_img_back[end:-1:1, end:-1:1], out_img_front)
out_img = map(clamp01nan, out_img)

ratio = (width_dest-400) / (2 * height_dest)
figure_size = ratio > 1 ? 
    (max_figure_size, ceil(Int, max_figure_size / ratio)) : 
    (ceil(Int, max_figure_size * ratio), max_figure_size)

canvas = plot(
    out_img, 
    xlims=(200, width_dest - 200),
    ylims=(1, 2 * height_dest),
    size=figure_size, 
    ticks=:none,
    )

gridlines = make_grid_lines(
    dest_affine, dest_proj, -180.0:30.0:180.0, -90:30.0:90; num_points=1000
    );
for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:white, label="", alpha=alpha)
    plot!(canvas, xs, ys .+ height_dest, color=:white, label="", alpha=alpha)
end

output_path = joinpath(output_dir, "transverse_mercator.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")