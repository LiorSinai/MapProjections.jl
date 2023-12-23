using GeoJSON
using Shapefile
using Plots
using PolygonAlgorithms

using MapProjections
include("geoJSON.jl")
include("shapefile_to_geojson.jl")

function plot_geometry!(canvas, geometry::GeoJSON.Polygon; options...)
    if isempty(geometry)
        return canvas
    end
    polygon = geometry.coordinates[1]
    plot!(canvas, [Shape(polygon)]; options...)
    for hole in geometry.coordinates[2:end]
        xs = map(xy->xy[1], hole)
        ys = map(xy->xy[2], hole)
        plot!(canvas, xs, ys; options...)
    end
    canvas
end

function plot_geometry!(canvas, geometry::GeoJSON.MultiPolygon; options...)
    if isempty(geometry)
        return canvas
    end
    for polygon in geometry.coordinates
        plot!(canvas, [Shape(polygon[1])]; options...)
        for hole in polygon[2:end]
            xs = map(xy->xy[1], hole)
            ys = map(xy->xy[2], hole)
            plot!(canvas, xs, ys; options...)
        end
    end
    canvas
end

## Config
data_dir = "../data/ne_10m_admin_0_countries"
filename = "ne_10m_admin_0_countries.shp"
output_dir = "../gallery/boundaries"
data_path = joinpath(data_dir, filename);
max_figure_size = 500
# Percentiles: 80th → 8e-6 ; 90th → 3.2e-5, 95th → 3.2e-4
min_earth_area = 3.2e-5 # On a cylinder with radius 1 => total area=4π
max_figure_size = 500

## Data
print("loading data at $data_path ...")
shape_data = Shapefile.Table(data_path) |> shapefile_to_geojson
println("done")
print("filtering data ...")
features = filter_features_by_size(shape_data, min_earth_area); 
println("done")
println("")

## Source
src_proj = WorldGeodeticSystem84()

## Projections
include("boundary_projections/azimuthal.jl")
include("boundary_projections/cylindrical_equal_area.jl")
include("boundary_projections/equirectangular.jl")
include("boundary_projections/mercator.jl")
include("boundary_projections/orthographic.jl")
include("boundary_projections/robinson.jl")
include("boundary_projections/sinusoidal.jl")
include("boundary_projections/transverse_mercator.jl")
include("boundary_projections/transverse_mercator_ellipsoidal.jl")