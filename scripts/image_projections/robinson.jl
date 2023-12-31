let src_proj=src_proj, 
    src_affine=src_affine,
    img=img,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("Robinson")

dest_proj = Robinson(;k=1.0, long0=0.0, interpolator=MapProjections.LinearSpline())
height, width = size(img)
width_dest, height_dest, dest_affine = 
    calculate_suggested_transform(src_proj, dest_proj, width, height, src_affine)
println("Source:      ($height, $width)")
println("Destination: ($height_dest, $width_dest)")

out_img = reproject_warp(
    img, (height_dest, width_dest), src_proj, dest_proj, src_affine, dest_affine
    ; method=bilinear
)
out_img = map(clamp01nan, out_img)

ratio = width_dest / height_dest
figure_size = ratio > 1 ? 
    (max_figure_size, ceil(Int, max_figure_size / ratio)) : 
    (ceil(Int, max_figure_size * ratio), max_figure_size)

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
    plot!(canvas, xs, ys, color=:white, label="", alpha=alpha)
end

output_path = joinpath(output_dir, "robinson.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end