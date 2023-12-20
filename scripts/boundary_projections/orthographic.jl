using MapProjections: degree_to_radian, angular_distance

function outside_angular_range(long0, lat0, coordinate; min_radius=0.5)
    longitude, latitude = coordinate
    longitude = degree_to_radian(longitude - long0)
    latitude = degree_to_radian(latitude)
    lat0 = degree_to_radian(lat0)
    c = angular_distance(lat0, longitude, latitude)
    #print(c, " ")
    (abs(c) > π/2) && abs(sin(c)) < min_radius  
end

function filter_angular_range(polygon::GeoJSON.Polygon{D, T}, long0, lat0; options...) where {D, T}
    regions = Vector{Vector{NTuple{D, T}}}()
    for region in polygon.coordinates
        outside = any(xy -> outside_angular_range(long0, lat0, xy; options...), region)
        if !outside
            push!(regions, region)
        end
    end
    GeoJSON.Polygon{D, T}(nothing, regions)
end

function filter_angular_range(multipolygon::GeoJSON.MultiPolygon{D, T}, long0, lat0; options...) where {D, T}
    collection = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in multipolygon.coordinates
        regions = Vector{Vector{NTuple{D, T}}}()
        for region in polygon
            outside = any(xy -> outside_angular_range(long0, lat0, xy; options...), region)
            if !outside
                push!(regions, region)
            end
        end
        if length(regions) > 0
            push!(collection, regions)
        end
    end
    GeoJSON.MultiPolygon{D, T}(nothing, collection)
end

function filter_angular_range(feature::GeoJSON.Feature, long0, lat0; options...)
    filter_angular_range(feature.geometry, long0, lat0; options...)
end

let src_proj=src_proj, 
    features=features,
    max_figure_size=max_figure_size,
    output_dir=output_dir

println("Orthographic")

dest_proj = Orthographic(;long0=0.0, lat0=15.0, radius=1.0)
back_long0 = dest_proj.long0 <= 0 ? dest_proj.long0 + 180 : dest_proj.long0 - 180
back_lat0 = -dest_proj.lat0
dest_proj_back = Orthographic(;long0=back_long0, lat0=back_lat0, radius=dest_proj.radius)

canvas_front = plot(aspectratio=:equal)
canvas_back = plot(aspectratio=:equal)
for (idx, shape) in enumerate(features)
    print("$idx, ")
    front = filter_angular_range(shape,  dest_proj.long0, dest_proj.lat0)
    projected_front = reproject(front, src_proj, dest_proj; clip=false)
    plot_geometry!(canvas_front, projected_front; label="", color=:black, fillalpha=0.3)
    back = filter_angular_range(shape, dest_proj_back.long0, dest_proj_back.lat0)
    projected_back = reproject(back, src_proj, dest_proj_back; clip=false)
    plot_geometry!(canvas_back, projected_back; label="", color=:black, fillalpha=0.3)
end
println("")

## boundary
num_points = 1000
t = range(0, 2π, length=num_points)
xs = dest_proj.radius * cos.(t)
ys = dest_proj.radius * sin.(t)
plot!(canvas_front, xs, ys, color=:black, label="");
plot!(canvas_back, xs, ys, color=:black, label="");

figure_size = (max_figure_size, ceil(Int, max_figure_size/2))
canvas = plot(
    canvas_front,
    canvas_back,
    plot_title="Orthographic",
    size=figure_size
    );

output_path = joinpath(output_dir, "orthographic.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")

end