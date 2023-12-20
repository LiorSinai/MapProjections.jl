let src_proj=src_proj, 
    features=features,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("AzimuthalEquidistant")

dest_proj = AzimuthalEquidistant(;long0=0.0, lat0=90.0, radius=1.0)

figure_size = (max_figure_size, max_figure_size)

canvas = plot(
    aspect_ratio=:equal, 
    title="Azimuthal Equidistant",
    size=figure_size
)
for (idx, shape) in enumerate(features)
    print("$idx, ")
    projected = reproject(shape, src_proj, dest_proj)
    plot_geometry!(canvas, projected.geometry; label="", color=:black, fillalpha=0.3)
end
println("")

## boundary
num_points = 1000
t = range(0, 2π, length=num_points)
xs = dest_proj.radius * π * cos.(t)
ys = dest_proj.radius * π * sin.(t)
plot!(canvas, xs, ys, color=:black, label="");

output_path = joinpath(output_dir, "azimuthal_equidistant.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

x = 5

end