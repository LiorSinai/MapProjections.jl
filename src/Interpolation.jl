module Interpolation

export interpolate, linear_interpolation
export LinearSpline, CubicSpline, SplineRoots
export polynomial, polynomial_grad, polynomial_root
import Base: inv, similar, show

function linear_interpolation(x1::Real, y1::Real, x2::Real, y2::Real, x::Real)
    grad = (y2 - y1) / (x2 - x1)
    inter = (y1 * x2 - y2 * x1) / (x2 - x1)
    grad * x + inter
end

function _get_interval_idx(
    nodes::AbstractVector, x::Real
    )
    idx = 2
    # assume nodes are ordered.
    # could do this with binary search in log(n) time.
    while (idx) < length(nodes) && (x > nodes[idx])
        idx += 1
    end
    idx # idx of first value greater than x
end

abstract type AbstractSpline end

"""
    interpolate(spline, x)

Evaluate the spline's polynomial `j` at `x`, where:

- `spline.intervals[j] <= x <= spline.intervals[j + 1]`.
- the coefficients are given by `spline.coefficients[:, j]`.

```
y = c0 + c1 * x + c2 * x^2 + ... + cn * x^n
```
"""
function interpolate(spline::AbstractSpline, x::Real)
    j = _get_interval_idx(spline.intervals, x) - 1
    coeffs = spline.coefficients[:, j]
    xj = spline.intervals[j]
    polynomial(coeffs, x - xj)
end


"""
    LinearSpline(intervals, values)

Linear interpolation between intervals.

```
    y = (y2-y1)/(x2-x1)*(x-x1) + y1 where xs[idx] < x < xs[idx+1]
```
"""
struct LinearSpline{V<:AbstractVector, M<:AbstractMatrix} <: AbstractSpline
    coefficients::M
    intervals::V
end

function LinearSpline(intervals::AbstractVector{T}, values::AbstractVector{T}) where T
    idxs = sortperm(intervals)
    values = values[idxs]
    intervals = intervals[idxs]
    n = length(intervals)
    coefficients = zeros(T, 2, n - 1)    
    for (i, j) in zip(1:(n-1), 2:n)
        grad =  (values[j] - values[i]) / (intervals[j] - intervals[i])
        inter = values[i]
        coefficients[:, i] = [inter, grad]
    end
    LinearSpline(coefficients, intervals)
end

LinearSpline() = LinearSpline(Matrix{Float64}(undef, 0, 0), Float64[])
similar(::LinearSpline, xs, ys) = LinearSpline(xs, ys)

inv(spline::LinearSpline) = LinearSpline(map(spline, spline.intervals), spline.intervals)

(spline::LinearSpline)(x::Real) = interpolate(spline, x)

function show(io::IO, mime::MIME"text/plain", spline::LinearSpline)
    print(io, "LinearSpline(")
    print(io, spline.coefficients, ", ")
    print(io, spline.intervals)
    print(io, ")")
end

"""
    CubicSpline(xs, ys)
    CubicSpline(coefficients, intervals)

The coefficients for cubic spline interpolation. Extra boundary condition that the second derivatives are zero.
    
The coefficients are stored as a matrix. Each column represents a polynomial ``j`` where the coefficient of ``x^i`` is in row ``i``.

That is:
```
    y = sum(A[i, j] * (x - x[j])^i) where xs[j] < x < xs[j+1]
```

Assumes `xs` are ordered.

E.g:
```
nodes = [-2.0, 1.0, 3.0, 7.0]
values = [5.0, 7.0, 11.0, 34.0]
spline = CubicSpline(nodes, values)
spline(5.0) # 20.698
```
"""
struct CubicSpline{M<:AbstractMatrix, V<:AbstractVector} <: AbstractSpline
    coefficients::M
    intervals::V
end

function CubicSpline(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    @assert length(xs) == length(ys)
    n = length(xs)
    col = 1:(n-1)
    col_next = 2:n 
    hs = xs[col_next] - xs[col]
    A = zeros(T, n, n)
    b = zeros(T, n)
    A[1, 1] = 1
    A[n, n] = 1
    for i in 1:(n-2)
        A[i + 1, i] = hs[i]
        A[i + 1, i + 1] = 2 * (hs[i] + hs[i + 1])
        A[i + 1, i + 2] = hs[i + 1]
        b[i + 1] = 3 * (ys[i + 2] - ys[i + 1]) / hs[i + 1] - 3 * (ys[i + 1] - ys[i])  / hs[i]
    end
    c = A \ b
    coefficients = zeros(T, 4, n - 1)
    coefficients[1, :] = ys[col]
    coefficients[2, :] = (ys[col_next] - ys[col]) ./ hs[col] - hs[col] / 3 .* (2 * c[col] + c[col_next])
    coefficients[3, :] = c[col]
    coefficients[4, :] = (c[col_next] - c[col]) ./ (3 * hs[col])
    CubicSpline{Matrix{T}, typeof(xs)}(coefficients, xs)
end

CubicSpline() = CubicSpline(zeros(0, 0), zeros(0))
similar(::CubicSpline, xs, ys) = CubicSpline(xs, ys)

inv(spline::CubicSpline) = 
    SplineRoots(spline.coefficients, spline.intervals, map(spline, spline.intervals))

(spline::CubicSpline)(x::Real) = interpolate(spline, x)
 
function show(io::IO, mime::MIME"text/plain", spline::CubicSpline)
    print(io, "CubicSpline")
    print(io, "\n  intervals = ", spline.intervals)
    print(io, "\n  coefficients = ")
    show(io, mime, spline.coefficients)
end

"""
    SplineRoots(coefficients, node_intervals, value_intervals)

The inverse of an `AbstractSpline`.
Instead of evaluating the polynomial `j` at `x`, this finds the value `x` to produce `y`,
where `x` is the root of polynomial `j`.

```
nodes = [-2.0, 1.0, 3.0, 7.0]
values = [5.0, 7.0, 11.0, 34.0]
spline = CubicSpline(nodes, values)
inv_spline = inv(spline)
inv_spline(20.698) # 4.99995
```
"""
struct SplineRoots{V1<:AbstractVector, V2<:AbstractVector}
    coefficients::Matrix
    node_intervals::V1
    value_intervals::V2
end


"""
    interpolate(inv_spline, y)

Given an `inv_spline`, find the root `x` of the relevant polynomial of the spline:

```
y = c0 + c1 * x + c2 * x^2 + ... + cn * x^n
```
"""
function interpolate(inv_spline::SplineRoots, y::Real; num_guesses::Int=5, newton_options...)
    j = _get_interval_idx(inv_spline.value_intervals, y) - 1
    coeffs = inv_spline.coefficients[:, j]
    coeffs[1] -= y
    xj = inv_spline.node_intervals[j]
    xk = inv_spline.node_intervals[j + 1]
    r = range(xj, xk, num_guesses)
    root = polynomial_root(coeffs, r, xj; newton_options...)
    if !(root >= xj) || !(root <= xk)
        throw(DomainError(root, "root not in interval range ($xj, $xk)"))
    end
    root
end

(inv_spline::SplineRoots)(x::Real) = interpolate(inv_spline, x)

## Polynomials

function polynomial(coefficients::AbstractVector, x::Real)
    result = 0.0
    xpow = 1.0
    for coeff in coefficients
        result += coeff * xpow
        xpow *= x
    end
    result
end

function polynomial_grad(coefficients::AbstractVector, x::Real)
    result = 0.0
    xpow = 1.0
    for (exponent_minus_1, coeff) in enumerate(coefficients[2:end])
        result += exponent_minus_1 * coeff * xpow
        xpow *= x
    end
    result
end

"""
    polynomial_root(coefficients, xmin, xmax, stepsize=)
"""
function polynomial_root(coefficients::AbstractVector, range_::AbstractRange, offset::Real=0.0; options...)
    x0 = argmin(x -> abs(polynomial(coefficients, x - offset)), range_)
    newtons_method_polynomial(coefficients, x0, offset; options...)
end

function newtons_method_polynomial(
    coefficients::AbstractVector, x0::Real, offset::Real=0.0
    ; atol::AbstractFloat=1e-6, max_iters::Int=10
    )
    xi = x0
    yi = polynomial(coefficients, xi - offset)
    iters = 0
    while iters <= max_iters && abs(yi) > atol
        iters += 1
        ygrad = polynomial_grad(coefficients, xi - offset)
        if abs(ygrad) < atol 
            return NaN # failure
        end
        xi = xi - yi / ygrad
        yi = polynomial(coefficients, xi - offset)
    end
    xi
end

end