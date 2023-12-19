using MapProjections
using Test

@testset verbose = true "MapProjections" begin
    include("interpolation.jl")
    include("projections.jl")
end