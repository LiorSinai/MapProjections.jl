let src_proj=src_proj, 
    features=features,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("Robinson")

dest_proj = Robinson(;radius=1.0, interpolator=MapProjections.CubicSpline())

ratio = 0.95 * 2
figure_size = (max_figure_size, ceil(Int, max_figure_size / ratio))

canvas = plot(
    aspect_ratio=:equal, 
    title="Robinson",
    size=figure_size,
)

for (idx, shape) in enumerate(features)
    print("$idx, ")
    geometry = shape.geometry
    projected = reproject(geometry, src_proj, dest_proj)
    plot_geometry!(canvas, projected; label="", color=:black, fillalpha=0.3)
end
println("")

## boundary
num_points = 1000
boundary = [
    [(-180.0 + dest_proj.long0, lat) for lat in range(-90.0, 90.0; length=num_points)]...,
    [(180.0 + dest_proj.long0, lat) for lat in range(90.0, -90.0; length=num_points)]...,
]
projected_boundary = reproject(boundary, src_proj, dest_proj)
plot!(canvas, Shape(projected_boundary), fillalpha=0.0, label="");

output_path = joinpath(output_dir, "robinson.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end