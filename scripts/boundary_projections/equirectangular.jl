println("Equirectangular")

figure_size = (max_figure_size, ceil(Int, 1.12*max_figure_size/2))

canvas = plot(
    aspect_ratio=:equal,
    size=figure_size,
    xlims=(-180, 180), ylims=(-90, 90),
    xlabel="longitude (°)",
    ylabel="latitude (°)",
    title="Equirectangular",
    margin=5Plots.mm,
    )

for (idx, shape) in enumerate(features)
    print("$idx, ")
    plot_geometry!(canvas, shape.geometry; label="", color=:black, fillalpha=0.3)
end
println("")

output_path = joinpath(output_dir, "equirectangular.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")