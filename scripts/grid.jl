using MapProjections: AffineTransform

function make_grid_lines(
    affine::AffineTransform, crs, range_long::StepRangeLen, range_lat::StepRangeLen
    ; num_points=1000
    )
    gridlines = Vector{Vector{Tuple{Float64, Float64}}}()
    lat_min, lat_max = extrema(range_lat)
    long_min, long_max = extrema(range_long)
    inv_affine = inv(affine)
    for long in range_long
        coords = map(lat -> crs((long, lat)), range(lat_min, lat_max, length=num_points))
        line = map(x -> inv_affine(x), coords)
        push!(gridlines, line)
    end
    for lat in range_lat
        coords = map(long -> crs((long, lat)), range(long_min, long_max, length=2*num_points))
        line = map(x -> inv_affine(x), coords)
        push!(gridlines, line)
    end
    gridlines
end