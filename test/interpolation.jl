using MapProjections.Interpolation

using Test

@testset verbose = true "Interpolation" begin

@testset "linear interpolation" begin
    x1, y1, x2, y2 = 0.6, 10.0, 0.8, 15.5

    # middle
    y_interp = linear_interpolation(x1, y1, x2, y2, 0.7)
    @test y_interp == 12.75

    # boundary
    y_interp = linear_interpolation(x1, y1, x2, y2, x1)
    @test y_interp == y1
    y_interp = linear_interpolation(x1, y1, x2, y2, x2)
    @test y_interp == y2

    # extrapolation
    y_interp = linear_interpolation(x1, y1, x2, y2, 1.5)
    @test y_interp ≈ 34.75
end

@testset "LinearSpline" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = [5.0, 7.0, 11.0, 34.0]

    spline = LinearSpline(nodes, values)
    ys = map(spline, nodes)
    @test ys == values

    y_interp = spline(0.0)
    @test y_interp == (6 + 1/3)
end

@testset "LinearSpline - reverse" begin
    nodes = 10:-1.0:0.0
    values = nodes.^2 .+ nodes .+ 3

    spline = LinearSpline(nodes, values)
    ys = map(spline, nodes)
    @test ys == values

    y_interp = spline(6.3)
    actual = 48.99
    @test isapprox(y_interp, 49.2)
end

@testset "CubicSpline - coefficients" begin
    nodes = [0.0, 1.0, 2.0, 2.5]
    values = [0.0, 1.0, 8.0, 9.0]

    spline = CubicSpline(nodes, values)
    expected = [
        0    11   88
        -12  57   48
        0    69   -78
        23   -49  52
    ]
    @test spline.coefficients * 11 ≈ expected 

    ys = map(spline, nodes)
    @test ys ≈ values
end

@testset "CubicSpline" begin
    nodes = [-2.0, 1.0, 3.0, 7.0]
    values = [5.0, 7.0, 11.0, 34.0]
    spline = CubicSpline(nodes, values)

    ys = map(spline, nodes)
    @test ys == values

    y_interp = spline(0.0)
    @test y_interp ≈ 6.0890804597701145

    inv_spline = SplineRoots(spline.coefficients, spline.intervals, map(spline, spline.intervals))
    xs = map(inv_spline, values)
    @test xs == nodes

    root = inv_spline(y_interp)
    @test isapprox(root, 0.0; atol=1e-9)

    y_interp = spline(5.32)
    @test y_interp ≈ 22.678032
    root = inv_spline(y_interp)
    @test isapprox(root, 5.32; atol=1e-6)
end

end