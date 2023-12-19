function filter_latitudes(
    geometry::GeoJSON.Polygon{D, T}, latitude_min::Float64, latitude_max::Float64
    ) where {D, T}
    coords = Vector{Vector{NTuple{D, T}}}()
    for region in geometry.coordinates
        coords_filtered = filter(x -> latitude_min <= x[2] <= latitude_max , region)
        push!(coords, coords_filtered)
    end
    GeoJSON.Polygon{D, T}(nothing, coords)
end

function filter_latitudes(
    geometry::GeoJSON.MultiPolygon{D, T}, latitude_min::Float64, latitude_max::Float64
    ) where {D, T}
    coords = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in geometry.coordinates
        regions = Vector{Vector{NTuple{D, T}}}()
        for region in polygon
            coords_filtered = filter(x -> latitude_min <= x[2] <= latitude_max, region)
            push!(regions, coords_filtered)
        end
        push!(coords, regions)
    end
    GeoJSON.MultiPolygon{D, T}(nothing, coords)
end

println("Mercator")

dest_proj = Mercator(;radius=1.0)

xmin, ymin = dest_proj((-180.0, -85.0));
xmax, ymax = dest_proj((180.0, +85.0));
width = xmax - xmin
height = ymax - ymin

ratio = 0.95 * width / height
figure_size = (max_figure_size, ceil(Int, max_figure_size / ratio))

canvas = plot(
    aspect_ratio=:equal, 
    xlims=(xmin, xmax), ylims=(ymin, ymax), 
    title="Mercator",
    size=figure_size,
)

for (idx, shape) in enumerate(features)
    print("$idx, ")
    geometry = shape.geometry
    geometry = filter_latitudes(geometry, -86.0, 86.0)
    projected = reproject(geometry, src_proj, dest_proj)
    plot_geometry!(canvas, projected; label="", color=:black, fillalpha=0.3)
end
println("")

output_path = joinpath(output_dir, "mercator.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")