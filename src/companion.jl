export Companion
using LinearAlgebra: dot

"""
[`Companion` matrix](http://en.wikipedia.org/wiki/Companion_matrix)

```jldoctest
julia> A = Companion([3,2,1])
3×3 Companion{Int64}:
 0  0  -3
 1  0  -2
 0  1  -1
```
Also, directly from a polynomial:

```julia
julia> using Polynomials

julia> P = Polynomial(2:5)
Polynomial(2 + 3x + 4x^2 + 5x^3)

julia> C = Companion(P)
3×3 Companion{Float64}:
 0.0  0.0  -0.4
 1.0  0.0  -0.6
 0.0  1.0  -0.8
```
"""
struct Companion{T} <: AbstractMatrix{T}
    c :: Vector{T}
end

Companion(v::AbstractVector) = Companion(collect(v))

# Construct companion matrix from a polynomial

function Companion(P::Polynomial)
    c = P.coeffs[begin:end-1] ./ P.coeffs[end]
    return Companion(c)
end

# Basic properties

function size(C::Companion)
    n = length(C.c)
    return n, n
end

@inline Base.@propagate_inbounds function getindex(
    C::Companion{T},
    i::Int,
    j::Int,
) where T
    @boundscheck checkbounds(C, i, j)
    return (j == length(C.c)) ? (@inbounds -C.c[i]) :
        i == j+1 ? one(T) : zero(T)
end


# Linear algebra

# 3-argument mul! mutates first argument: y <= C * x
function mul!(y::Vector, C::Companion, x::AbstractVector)
    @boundscheck length(y) == length(x) == size(C, 1) ||
        throw(DimensionMismatch("mul! arguments incompatible sizes"))
    z = x[end]
    y[1] = -C.c[1] * z
    y[2:end] = x[begin:end-1] - C.c[2:end] * z
    return y
end

# A <= B * C
function mul!(A::Matrix, B::AbstractMatrix, C::Companion)
    @boundscheck (size(A) == (size(B,1), size(C,2)) && size(B, 2) == size(C,1)) ||
        throw(DimensionMismatch("mul! arguments incompatible sizes"))
    Base.require_one_based_indexing(B)
    @views for j in 1:size(A,2)-1
        @inbounds A[:,j] = B[:,j+1]
    end
    @inbounds mul!((@view A[:,end]), B, C.c, -1, 0) # A[:,end] <= - B * c
    return A
end


function inv(C::Companion{T}) where T
    M = zeros(T, size(C)...)
    for i in 1:size(C,1)-1
        M[i, i+1] = one(T)
    end
    d = M[end, 1] = -one(T) / C.c[1]
    M[1:end-1, 1] = d * C.c[2:end]
    return M
end
