using GeoJSON
using LinearAlgebra: norm
using PolygonAlgorithms: area_polygon
import PolygonAlgorithms: bounds
import Base: similar
import MapProjections: reproject

function filter_features(shape_data::GeoJSON.FeatureCollection, property_name::Symbol, property_set::Set)
    features = GeoJSON.Feature[]
    for feature in shape_data.features
        if getproperty(feature, property_name) in property_set
            push!(features, feature)
        end
    end
    features
end

function filter_largest_features(shape_data::GeoJSON.FeatureCollection; keep_top::Int=1)
    features = GeoJSON.Feature[]
    for feature in shape_data.features
        if feature.geometry isa GeoJSON.MultiPolygon
            areas = map(coords -> area_polygon(coords[1]), feature.geometry)
            idxs = sortperm(areas; rev=true)
            idx_max = min(keep_top, length(idxs))
            new_coords = feature.geometry.coordinates[idxs[1:idx_max]]
            new_geom = GeoJSON.MultiPolygon(; bbox=nothing, coordinates=new_coords)
            new_feature = similar(feature, new_geom)
            push!(features, new_feature)
        else
            push!(features, feature)
        end
    end
    features
end

function filter_features_by_size(
    shape_data::GeoJSON.FeatureCollection, min_area::AbstractFloat)
    features = GeoJSON.Feature[]
    for feature in shape_data.features
        if feature.geometry isa GeoJSON.MultiPolygon
            areas = Vector{Float64}(undef, length(feature.geometry))
            for (idx, polygon) in enumerate(feature.geometry.coordinates)
                # only work with the first region
                xy = map(project_to_cylinder, polygon[1])
                areas[idx] = area_polygon(xy)
            end
            idxs = areas .>= min_area
            new_coords = feature.geometry.coordinates[idxs]
            new_geom = GeoJSON.MultiPolygon(; bbox=nothing, coordinates=new_coords)
            new_feature = similar(feature, new_geom)
            push!(features, new_feature)
        else
            polygon = feature.geometry.coordinates[1] # ignore holes
            area = area_polygon(map(project_to_cylinder, polygon))
            if area >= min_area
                push!(features, feature)
            end
        end
    end
    features
end

function project_to_cylinder(long_lat::Tuple{T, T}) where T <: AbstractFloat
    (long_lat[1] * π / 180, sin(long_lat[2] * π / 180))
end

function similar(feature::GeoJSON.Feature{D, T}, geometry::GeoJSON.AbstractGeometry) where {D, T}
    id = getfield(feature, :id)
    bbox = getfield(feature, :bbox)
    GeoJSON.Feature{D, T}(id, bbox, geometry, feature.properties)
end

"""
    translate(feature, offset)
    translate(multipolygon, offset)
    translate(polygon, offset)

Return a similar GeoJSON object with co-ordinates translated by `offset`.
"""
function translate end

function translate(polygon::GeoJSON.Polygon{D, T}, offset::Tuple) where {D, T}
    new_regions = Vector{NTuple{D, T}}[]
    for region in polygon.coordinates
        coords = [(xy .+ offset) for xy in region]
        push!(new_regions, coords)
    end
    GeoJSON.Polygon{D, T}(nothing, new_regions)
end

function translate(multipolygon::GeoJSON.MultiPolygon{D, T}, offset::Tuple) where {D, T}
    collection = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in multipolygon.coordinates
        new_regions = Vector{Vector{NTuple{D, T}}}()
        for region in polygon
            coords = [(xy .+ offset) for xy in region]
            push!(new_regions, coords)
        end
        push!(collection, new_regions)
    end
    GeoJSON.MultiPolygon{D, T}(nothing, collection)
end

function translate(feature::GeoJSON.Feature, offset::NTuple)
    new_geom = translate(feature.geometry, offset)
    similar(feature, new_geom)
end

"""
    bounds(feature_collection)
    bounds(multipolygon)
    bounds(polygon)

Return rectangular bounds for a GeoJSON object as `(xmin, ymin, xmax, ymax)`.
"""
function bounds(shape_data::GeoJSON.FeatureCollection)
    xmin, xmax, ymin, ymax = bounds(shape_data.features[1].geometry)
    for feature in shape_data.features
        x0, y0, x1, y1 = bounds(feature.geometry)
        xmin = min(xmin, x0, x1)
        xmax = max(xmax, x0, x1)
        ymin = min(ymin, y0, y1)
        ymax = max(ymax, y0, y1)
    end
    (xmin, ymin, xmax, ymax)
end

function bounds(polygon::GeoJSON.Polygon)
    start = polygon[1][1]
    xmin = start[1]
    xmax = start[1]
    ymin = start[2]
    ymax = start[2]
    for region in polygon
        x0, y0, x1, y1 = bounds(region)
        xmin = min(xmin, x0, x1)
        xmax = max(xmax, x0, x1)
        ymin = min(ymin, y0, y1)
        ymax = max(ymax, y0, y1)
    end
    (xmin, ymin, xmax, ymax)
end

function bounds(multipolygon::GeoJSON.MultiPolygon)
    start = multipolygon[1][1][1]
    xmin = start[1]
    xmax = start[1]
    ymin = start[2]
    ymax = start[2]
    for polygon in multipolygon
        for region in polygon
            x0, y0, x1, y1 = bounds(region)
            xmin = min(xmin, x0, x1)
            xmax = max(xmax, x0, x1)
            ymin = min(ymin, y0, y1)
            ymax = max(ymax, y0, y1)
        end
    end
    (xmin, ymin, xmax, ymax)
end

"""
    reproject(feature, src_crs, dest_crs)
    reproject(multipolygon, src_crs, dest_crs)
    reproject(polygon, src_crs, dest_crs)

Reproject coordinates in the `src_crs` to the `dest_crs`.

There must be a definition for `inv(src_crs)`.
"""
function reproject(polygon::GeoJSON.Polygon{D, T}, src_crs, dest_crs; options...) where {D, T}
    coords = [reproject(region, src_crs, dest_crs; options...) for region in polygon.coordinates]
    GeoJSON.Polygon{D, T}(nothing, coords)
end

function reproject(multipolygon::GeoJSON.MultiPolygon{D, T}, src_crs, dest_crs; options...) where {D, T}
    collection = Vector{Vector{Vector{NTuple{D, T}}}}()
    for polygon in multipolygon.coordinates
        projected_polygon = Vector{Vector{NTuple{D, T}}}()
        for region in polygon
            coords = reproject(region, src_crs, dest_crs; options...)
            push!(projected_polygon, coords)
        end
        push!(collection, projected_polygon)
    end
    GeoJSON.MultiPolygon{D, T}(nothing, collection)
end
   
function reproject(feature::GeoJSON.Feature{D, T}, src_crs, dest_crs; options...) where {D, T}
    geometry = reproject(feature.geometry, src_crs, dest_crs; options...)
    id = getfield(feature, :id)
    bbox = getfield(feature, :bbox)
    GeoJSON.Feature{D, T}(id, bbox, geometry, feature.properties)
end

