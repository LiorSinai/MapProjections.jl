using Shapefile
using Shapefile.Tables

function shapefile_to_dict(shape_data::Shapefile.Table)
    features = Dict{Symbol, Any}[]
    for (idx, row) in enumerate(shape_data)
        properties = Dict{Symbol, Any}()
        for property_name in propertynames(shape_data)
            if property_name == :geometry
                continue
            end
            properties[property_name] = _format_property(getproperty(row, property_name))
        end
        geometry = _get_geometry(Shapefile.shape(row))
        feature = Dict{Symbol, Any}(
            :type => "Feature",
            :id => idx,
            :properties => properties,
            :geometry => geometry
        )
        push!(features, feature)
    end
    Dict(
        :type => "FeatureCollection",
        :features => [features[1]]
    )
end

function _get_geometry(geometry::Shapefile.Polygon)
    coordinates = Vector{Vector{NTuple{2, Float32}}}()
    parts = vcat(geometry.parts .+ 1, length(geometry.points))
    for (idx_start, idx_end) in zip(parts[1:(end-1)], parts[2:end])
        region = map(p->(p.x, p.y), geometry.points[idx_start:idx_end])
        push!(coordinates, region)
    end
    Dict{Symbol, Union{String, Vector{Vector{NTuple{2, Float32}}}}}(
        :type => "Polygon",
        :coordinates => coordinates,
    )
end

function _format_property(property::String)
    replace(property, "\0"=>"")
end
    
_format_property(property) = property
