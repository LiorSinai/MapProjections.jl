function bounds(coordinates::AbstractVector{Tuple{T, T}}) where T
    start = coordinates[1]
    xmin = start[1]
    xmax = start[1]
    ymin = start[2]
    ymax = start[2]
    for (x, y) in coordinates
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    (xmin, ymin, xmax, ymax)
end