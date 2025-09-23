let src_proj=src_proj, 
    src_affine=src_affine,
    img=img,
    max_figure_size=max_figure_size,
    output_dir=output_dir
    
println("Mercator")

lat_limit = 85.0
j, i = inv(src_affine)((0.0, lat_limit))
i_cut = floor(Int, i)
height, width = size(img)
img_cropped = img[i_cut:(height - i_cut), 1:end];
height, width = size(img_cropped)
crop_affine = affine_from_bounds(long_min, -lat_limit, long_max, +lat_limit, width, height)

dest_proj = Mercator(;k=1.0)
height, width = size(img_cropped)
width_dest, height_dest, dest_affine = 
    calculate_suggested_transform(src_proj, dest_proj, width, height, crop_affine)
println("Source:      ($height, $width)")
println("Destination: ($height_dest, $width_dest)")

out_img = reproject_warp(
    img_cropped, (height_dest, width_dest), src_proj, dest_proj, crop_affine, dest_affine
    ; method=bilinear
)
out_img = map(clamp01nan, out_img)

figure_size = (max_figure_size, max_figure_size)

canvas = plot(
    out_img, 
    xlims=(1, width_dest),
    ylims=(1, height_dest),
    size=figure_size, 
    ticks=:none,
    )

gridlines = make_grid_lines(
    dest_affine, dest_proj, -180.0:30.0:180.0, -90:30.0:90; num_points=1000
    );

for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:white, label="", alpha=α)
end

output_path = joinpath(output_dir, "mercator.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end