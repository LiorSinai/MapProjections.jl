using MapProjections: bound_longitude

function filter_longitudes(
    geometry::GeoJSON.Polygon{D, T},
    longitude_min::Float64,
    longitude_max::Float64,
    long0::AbstractFloat
    ) where {D, T}
    coords = Vector{Vector{NTuple{D, T}}}()
    for region in geometry.coordinates
        coords_filtered = filter(x -> longitude_min <= bound_longitude(x[1] - long0) <= longitude_max , region)
        push!(coords, coords_filtered)
    end
    bbox = geometry.bbox
    GeoJSON.Polygon{D, T}(bbox, coords)
end

function filter_longitudes(
    geometry::GeoJSON.MultiPolygon{D, T},
    longitude_min::Float64,
    longitude_max::Float64,
    long0::AbstractFloat
    ) where {D, T}
    coords = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in geometry.coordinates
        regions = Vector{Vector{NTuple{D, T}}}()
        for region in polygon
            coords_filtered = filter(x -> longitude_min <= bound_longitude(x[1] - long0) <= longitude_max , region)
            push!(regions, coords_filtered)
        end
        push!(coords, regions)
    end
    bbox = geometry.bbox
    GeoJSON.MultiPolygon{D, T}(bbox, coords)
end

println("EllipsoidalTransverseMercator")

# Antartica not rendered properly for abs(long)>=135.0
dest_proj = EllipsoidalTransverseMercator(;long0=0.0, semi_major_axis=1.0)

canvas = plot(
    aspect_ratio=:equal, 
    title="Ellipsoidal Transverse\nMercator $(dest_proj.long0)°",
)
west_long = -45.0 + dest_proj.long0
east_long = dest_proj.long0 + 45.0
for (idx, shape) in enumerate(features)
    print("$idx, ")
    geometry = shape.geometry
    geometry = filter_longitudes(geometry, -45.0, 45.0, dest_proj.long0)
    projected = reproject(geometry, src_proj, dest_proj)
    plot_geometry!(canvas, projected, fillalpha=0.3, color=:black, label="")
end
println("")

boundaries = [
    [(west_long, lat) for lat in -90:0.1:90 ],
    [(east_long, lat) for lat in -90:0.1:90 ]
]
boundaries_projected = reproject(boundaries, src_proj, dest_proj)
for line in boundaries_projected
    xs = [xy[1] for xy in line]
    ys = [xy[2] for xy in line]
    plot!(canvas, xs, ys, color=:black, label="")
end

xlimits = xlims(canvas)
ylimits = ylims(canvas)
width = xlimits[2] - xlimits[1]
height = ylimits[2] - ylimits[1]
ratio = width / height
figure_size = ratio > 1 ? 
    (max_figure_size, ceil(Int, max_figure_size / ratio)) : 
    (ceil(Int, max_figure_size * ratio), max_figure_size)
plot!(canvas, size=figure_size)

output_path = joinpath(output_dir, "transverse_mercator_ellipsoidal.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")