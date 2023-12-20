let src_proj=src_proj, 
    src_affine=src_affine,
    img=img,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("Equirectangular")

figure_size = (max_figure_size, ceil(Int, max_figure_size/2))
height, width = size(img)

canvas = plot(
    img, 
    xlims=(1, width),
    ylims=(1, height),
    size=figure_size, 
    ticks=:none,
    )

gridlines = make_grid_lines(
    src_affine, src_proj, -180.0:30.0:180.0, -90:30.0:90;num_points=1000
    );

for line in gridlines
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:white, label="", alpha=0.5)
end

output_path = joinpath(output_dir, "equirectangular.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end