using Shapefile
using Shapefile.Tables
using GeoJSON
using PolygonAlgorithms: is_counter_clockwise

# Shapefiles are more space efficient but GeoJSON is easier to work with
function shapefile_to_geojson(shape_data::Shapefile.Table)
    features = GeoJSON.Feature{2, Float64}[]
    for (idx, row) in enumerate(shape_data)
        properties = Dict{Symbol, Any}()
        for property_name in propertynames(shape_data)
            if property_name == :geometry
                continue
            end
            properties[property_name] = _format_property(getproperty(row, property_name))
        end
        geometry = _format_geometry(Shapefile.shape(row))
        feature = GeoJSON.Feature(;id=idx, geometry=geometry, properties=properties, bbox=geometry.bbox)
        push!(features, feature)
    end
    GeoJSON.FeatureCollection(;features=features)
end

function _format_geometry(geometry::Shapefile.Polygon)
    parts = vcat(geometry.parts .+ 1, length(geometry.points))
    mbr = geometry.MBR
    bounding_box = [mbr.left, mbr.bottom, mbr.right, mbr.top]
    points = map(p -> (p.x, p.y), geometry.points)
    ## The Shapefile specification does not have a MultiPolygon equivalent.
    ## Simple checks are done here to diffentiate them from Polygons.
    if length(parts) == 2
        return GeoJSON.Polygon(; bbox=bounding_box, coordinates=[points])
    end
    coordinates = Vector{Vector{Vector{NTuple{2, Float64}}}}()
    for (idx_start, idx_end) in zip(parts[1:(end-1)], parts[2:end])
        region = points[idx_start:(idx_end - 1)]
        if is_hole(region)
            # holes must always come after the main polygon
            push!(coordinates[end], region)
        else
            push!(coordinates, [region])
        end
    end
    if length(coordinates) == 1
        return GeoJSON.Polygon(; bbox=bounding_box, coordinates=coordinates[1])
    end
    GeoJSON.MultiPolygon(; bbox=bounding_box, coordinates=coordinates)
end

is_hole(vertices::Vector{<:Tuple}) = is_counter_clockwise(vertices)

function _format_property(property::String)
    property = replace(property, "\0"=>"")
    isempty(property) ? nothing : property
end
    
_format_property(property) = property
