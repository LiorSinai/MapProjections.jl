let src_proj=src_proj, 
    features=features,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("CylindricalEqualArea")

dest_proj = CylindricalEqualArea(;k=1.0, radius=1.0)

xmin, ymax = dest_proj((-180.0, 90.0));
xmax, ymin = dest_proj((180.0, -90.0));

ratio = 0.95 *2π / 2
figure_size = (max_figure_size, ceil(Int, max_figure_size / ratio))

canvas = plot(
    aspect_ratio=:equal, 
    xlims=(xmin, xmax), ylims=(ymin, ymax), 
    title="Cylindrical Equal Area",
    size=figure_size
)
for (idx, shape) in enumerate(features)
    print("$idx, ")
    projected = reproject(shape, src_proj, dest_proj)
    plot_geometry!(canvas, projected.geometry; label="", color=:black, fillalpha=0.3)
end
println("")

output_path = joinpath(output_dir, "cylindrical_equal_area.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end