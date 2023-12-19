function split_discontinuities(vertices::AbstractVector{NTuple{D, T}}, distance_max::AbstractFloat) where {D, T}
    idx_prev = 1
    idxs = Vector{Vector{Int}}()
    for (idx, pair) in enumerate(zip(vertices[1:(end-1)], vertices[2:end]))
        d = norm(pair[2] .- pair[1])
        if d > distance_max
            #print("breaking at $idx ")
            push!(idxs, idx_prev:idx)
            idx_prev = idx + 1
        end
    end
    push!(idxs, idx_prev:length(vertices))
    idxs = merge_splits(vertices, idxs, distance_max)
    [[vertices[idx]] for idx in idxs] # no support for holes
end

function merge_splits(
        vertices::AbstractVector{NTuple{D, T}}, idxs::Vector{Vector{Int}}, distance_max::AbstractFloat
    ) where {D, T}
    groups = Dict{Int, Tuple{NTuple{D, T}, Vector{Int}}}()
    for i in 1:length(idxs)
        merged = false
        avg_point = reduce(.+, vertices[idxs[i]]) ./ length(idxs[i])
        for j in keys(groups)
            d = norm(avg_point .- groups[j][1])
            if d < distance_max
                #print("merging $i and $j ")
                push!(groups[j][2], idxs[i]...)
                merged = true
                break
            end
        end
        if !merged
            groups[i] = (avg_point, idxs[i])
        end
    end
    map(avg_idx->avg_idx[2], values(groups))
end

function split_discontinuities(polygon::GeoJSON.Polygon{D, T}, distance_max::AbstractFloat) where {D, T}
    splits = split_discontinuities(polygon.coordinates[1], distance_max) # no support for holes
    GeoJSON.MultiPolygon{D, T}(nothing, splits)
end

function split_discontinuities(multipolygon::GeoJSON.MultiPolygon{D, T}, distance_max::AbstractFloat) where {D, T}
    coords = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in multipolygon.coordinates
        splits = split_discontinuities(polygon[1], distance_max) # no support for holes
        push!(coords, splits...)
    end
    GeoJSON.MultiPolygon{D, T}(nothing, coords)
end

function add_infinity!(coordinates::Vector{<:Tuple}; xmax=5.0)
    print("adding infinity ... ")
    idx = 1
    diff_max = 0.0
    for i in 2:length(coordinates)
        y_diff = abs(coordinates[i][2] - coordinates[i - 1][2])
        if y_diff > diff_max
            diff_max = y_diff
            idx = i
        end
    end
    is_increasing = coordinates[idx][2] > coordinates[idx-1][2]
    ymax = max(coordinates[idx][2], coordinates[idx - 1][2])
    ymin = min(coordinates[idx][2], coordinates[idx - 1][2])
    y1 = is_increasing ? ymin : ymax
    y2 = is_increasing ? ymax : ymin
    insert!(coordinates, idx, (xmax, y1))
    insert!(coordinates, idx + 1, (xmax, y2))
end

function maybe_add_infinity!(geometry::GeoJSON.Polygon,  projected::GeoJSON.Polygon, pole; options...)
    added_infinity = false
    for (region, region_proj) in zip(geometry.coordinates, projected.coordinates)
        if contains(region, pole)
            add_infinity!(region_proj; options...)
            added_infinity = true
        end
    end
    added_infinity
end

function maybe_add_infinity!(geometry::GeoJSON.MultiPolygon, projected::GeoJSON.MultiPolygon, pole; options...)
    added_infinity = false
    for (polygon, polygon_proj) in zip(geometry.coordinates, projected.coordinates)
        for (region, region_proj) in zip(polygon, polygon_proj)
            if contains(region, pole)
                added_infinity = add_infinity!(region_proj; options...)
                added_infinity = true
            end
        end
    end
    added_infinity
end

println("TransverseMercator")

# Doesn't work well for all angles. (-170, -140) is problematic
dest_proj = TransverseMercator(;long0=0.0, radius=1.0)
extend = true

canvas = plot(
    aspect_ratio=:equal, 
    ylims=(-π/2, 3π/2), 
    title="Transverse Mercator $(dest_proj.long0)°",
)
west_pole = (-90.0 + dest_proj.long0, 0.0)
east_pole = (90.0 + dest_proj.long0, 0.0)
xmax=5.0
for (idx, shape) in enumerate(features)
    print("$idx, ")
    geometry = shape.geometry
    projected = reproject(geometry, src_proj, dest_proj; extend=extend)
    added_infinity1 = maybe_add_infinity!(geometry, projected, west_pole; xmax=-xmax)
    added_infinity2 = maybe_add_infinity!(geometry, projected, east_pole; xmax=+xmax)
    if !(added_infinity1 || added_infinity2)
        projected = split_discontinuities(projected, 2 * dest_proj.radius)
    end
    plot_geometry!(canvas, projected, fillalpha=0.3, color=:black, label="")
end
println("")
plot!([0, 0], [-π/2, π/2+π], c=:black, linestyle=:dash, label="");
xlimits = xlims(canvas)
xlimits = (max(xlimits[1], -xmax), min(xlimits[2], xmax))
plot!(collect(xlimits), [π/2, π/2], c=:black, linestyle=:dash, label="", xlims=xlimits);

width = xlimits[2] - xlimits[1]
height = 2π
ratio = width / height
figure_size = ratio > 1 ? 
    (max_figure_size, ceil(Int, max_figure_size / ratio)) : 
    (ceil(Int, max_figure_size * ratio), max_figure_size)
plot!(canvas, size=figure_size)

output_path = joinpath(output_dir, "transverse_mercator.png")
result = savefig(canvas, output_path)
println("saved image to $result")
println("")