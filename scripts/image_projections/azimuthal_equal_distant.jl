let src_proj=src_proj, 
    src_affine=src_affine,
    img=img,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("AzimuthalEquidistant")

dest_proj = AzimuthalEquidistant(;lat0=90.0, long0=0.0)
width_dest, height_dest = 500, 500

height, width = size(img)
r = π * dest_proj.radius
dest_affine = affine_from_bounds(-r, -r, r, r, width_dest, height_dest);

out_img = reproject_warp(
    img, (height_dest, width_dest), src_proj, dest_proj, src_affine, dest_affine
    ; method=bilinear
)
out_img = map(clamp01nan, out_img)

canvas = plot(
    out_img, 
    xlims=(1, width_dest),
    ylims=(1, height_dest),
    size=(max_figure_size, max_figure_size), 
    ticks=:none,
    )

gridlines = make_grid_lines(
    dest_affine, dest_proj, -180.0:30.0:180.0, -60:30.0:90; num_points=1000
    );

for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:white, label="", alpha=α)
end

output_path = joinpath(output_dir, "azimuthal_equidistant.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end