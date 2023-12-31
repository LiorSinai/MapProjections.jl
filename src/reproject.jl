"""
    reproject(coordinate, src_crs, dest_crs; options)
    reproject(coordinates, src_crs, dest_crs; options)

Reproject coordinates in the `src_crs` to the `dest_crs`.

There must be a definition for `inv(src_crs)`.
"""
function reproject end

function reproject(coordinate::Tuple, src_crs, dest_crs; options...)
    inv_src_crs = inv(src_crs)
    long_lat = inv_src_crs(coordinate)
    dest_crs(long_lat; options...)
end

function reproject(coordinates::Vector{<:Tuple}, src_crs, dest_crs; options...)
    inv_src_crs = inv(src_crs)
    long_lats = map(inv_src_crs, coordinates)
    map(coord -> dest_crs(coord; options...), long_lats)
end

function reproject(polygon::Vector{Vector{NTuple{D, T}}}, src_crs, dest_crs; options...) where {D, T}
    coordinates = Vector{Vector{NTuple{D, T}}}()
    for region in polygon
        projected = reproject(region, src_crs, dest_crs; options...)
        push!(coordinates, projected)
    end
    coordinates
end

function reproject(multipolygon::Vector{Vector{Vector{NTuple{D, T}}}}, src_crs, dest_crs; options...) where {D, T}
    coordinates = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in multipolygon.coordinates
        projected_polygon = reproject(polygon, src_crs, dest_crs; options...)
        push!(coordinates, projected_polygon)
    end
    coordinates
end

"""
    calculate_suggested_transform(src_crs, dest_crs, width, height, src_transform; num_steps=20)

Inspired by GDAL's `GDALSuggestedWarpOutput()`. Outputs an image shape and a transform for warping.

The scale is computed so that the distance from the top left corner to the bottom right corner is the same number of pixels as the original image.
This is intended to approximately preserve the resolution of the input data in the output file.

Samples are taking on the borders with approximately `num_steps` per border.
The border is assumed to defined the extent of the projection.
Note this assumption does not hold for all projections.
"""
function calculate_suggested_transform(
    src_crs, dest_crs, width::Int, height::Int, src_affine::AffineTransform; num_steps::Int=20
    )
    step = max((width - 1)/ num_steps, (height - 1)/num_steps)
    inv_src_crs = inv(src_crs)
    transform = xy -> xy |> src_affine |> inv_src_crs |> dest_crs
    width = convert(Float64, width)
    height = convert(Float64, height)
    ## sample points on the boundaries.
    corners = map(transform, [(1.0, 1.0), (width, 1.0), (width, height), (1.0, height)])
    xs = (1 + step):step:(width - step)
    ys = (1 + step):step:(height - step)
    top = [transform((j, 1.0)) for j in xs]
    bottom = [transform((j, height)) for j in xs]
    left = [transform((1.0, i)) for i in ys]
    right = [transform((width, i)) for i in ys]
    points = vcat(corners, top, bottom, left, right)
    filter!(p -> all(isfinite.(p)), points)
    ## get bounds
    xmin, ymin, xmax, ymax = _bounds(points)
    width_dest = xmax - xmin
    height_dest = ymax - ymin
    ratio = width_dest / height_dest
    ## Same distance for diagonal:
    ## hs*hs + ws*ws = hd*hd + wd*wd = hd*hd + (hd*r)*(hd*r)
    diagonal_src = width * width + height * height
    height_dest = floor(Int, sqrt( diagonal_src / (ratio * ratio + 1.0)))
    width_dest = floor(Int, ratio * height_dest)
    dest_affine = affine_from_bounds(xmin, ymin, xmax, ymax, width_dest, height_dest)
    width_dest, height_dest, dest_affine
end

function _bounds(points::Vector{<:Tuple{T, T}}) where T
    start = points[1]
    xmin = start[1]
    xmax = start[1]
    ymin = start[2]
    ymax = start[2]
    for (x, y) in points
        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
    end
    (xmin, ymin, xmax, ymax)
end

"""
    reproject_warp(
        img, inds, src_crs, dest_crs, src_affine, out_affine, src_transform
        ; options...
    )

Wrapper around `ImageTransformations.warp`.

Transform the coordinates of `img`, returning a new `imgw` satisfying `imgw[I] = img[tform(I)]`.

There must be a definition for `inv(dest_crs)`.
"""
function reproject_warp(
    img::AbstractArray, 
    inds::Tuple, 
    src_crs,
    dest_crs,
    src_affine::AffineTransform,
    dest_affine::AffineTransform,
    ; options...
    )
    # `imgw[I] = img[tform(I)]`
    # `I` are indices in `imgw`. The steps to convert `I` back to indices in `img` are:
    # 1. reverse - Reverse indices: (i, j) to (j, i) because require (x, y) tuples.
    # 2. dest_affine - Convert to destination CRS co-ordinates: (j, i) to (x, y).
    # 3. inv_dest_crs - Convert to spherical co-ordinates: (x, y) to (long, lat).
    # 4. src_crs - convert to source CRS co-ordinates: (long, lat) to (x, y)
    # 5. invA - Convert to pixels: (x, y) to (j, i).
    # 6. reverse - Reverse: (j, i) to (i, j).
    invA = inv(src_affine)
    inv_dest_crs = inv(dest_crs)
    tform = idx -> idx |> reverse |> dest_affine |> inv_dest_crs |> src_crs |> invA |> reverse
    out_img = ImageTransformations.warp(img, tform, inds; options...);
    out_img
end
