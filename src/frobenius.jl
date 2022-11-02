export Frobenius

"""
    [`Frobenius` matrix](http://en.wikipedia.org/wiki/Frobenius_matrix)

Frobenius matrices or Gaussian elimination matrices of the form

```
[ 1 0 ...     0 ]
[ 0 1 ...     0 ]
[ .........     ]
[ ... 1 ...     ]
[ ... c1 1 ...  ]
[ ... c2 0 1 ...]
[ ............. ]
[ ... ck ...   1]
```

i.e., an identity matrix with possibly nonzero subdiagonal elements
along a single column.

In this implementation,
the subdiagonal of the nonzero column is stored as a dense vector,
so that the size can be inferred automatically as j+k,
where j is the column index and k is the number of subdiagonal elements.

```jldoctest
julia> using SpecialMatrices

julia> F = Frobenius(3, 4:6) # Specify subdiagonals of column 3
6×6 Frobenius{Int64}:
 1  0  0  0  0  0
 0  1  0  0  0  0
 0  0  1  0  0  0
 0  0  4  1  0  0
 0  0  5  0  1  0
 0  0  6  0  0  1

julia> inv(F) # Special form of inverse
6×6 Frobenius{Int64}:
 1  0   0  0  0  0
 0  1   0  0  0  0
 0  0   1  0  0  0
 0  0  -4  1  0  0
 0  0  -5  0  1  0
 0  0  -6  0  0  1

julia> F*F # Special form preserved if the same column has the subdiagonals
6×6 Frobenius{Int64}:
 1  0   0  0  0  0
 0  1   0  0  0  0
 0  0   1  0  0  0
 0  0   8  1  0  0
 0  0  10  0  1  0
 0  0  12  0  0  1

julia> F*Frobenius(2, 2:5) # Promotes to Matrix
6×6 Matrix{Int64}:
 1   0  0  0  0  0
 0   1  0  0  0  0
 0   2  1  0  0  0
 0  11  4  1  0  0
 0  14  5  0  1  0
 0  17  6  0  0  1

julia> F * [10.0,20,30,40,50,60.0]
6-element Vector{Float64}:
  10.0
  20.0
  30.0
 160.0
 200.0
 240.0
```
"""
struct Frobenius{T} <: AbstractMatrix{T}
    colidx :: Int
    c :: Vector{T}
    function Frobenius(colidx::Int, c::AbstractVector{T}) where {T <: Number}
        1 ≤ colidx || throw("Bad column index $colidx")
        return new{T}(colidx, collect(c))
    end
end

# size
function size(F::Frobenius)
    n = F.colidx + length(F.c)
    return n, n
end

# getindex
@inline Base.@propagate_inbounds function getindex(F::Frobenius{T}, i::Integer, j::Integer) where {T}
    @boundscheck checkbounds(F, i, j)
    return (i == j) ? one(T) :
        (i > j == F.colidx) ?
        (@inbounds F.c[i-F.colidx]) :
        zero(T)
end


# Linear algebra

# 3-argument mul! mutates first argument: y <= F * x
# *(F, x) = F * x derives from this automatically in Base
function mul!(y::Vector, F::Frobenius, x::AbstractVector)
    Base.require_one_based_indexing(x)
    @boundscheck (n = size(F, 2)) == length(x) == length(y) ||
        throw(DimensionMismatch("$n $(length(y)) $(length(x))"))
    j = F.colidx
    copyto!(y, x)
    @inbounds @. (@view y[(j+1):end]) += F.c * x[j]
    return y
end

# This product is not type stable (result could be Frobenius or Matrix).
#=
It could be stable if we put the column index "j" as a type parameter.
Here is the basic idea:
function *(F::Frobenius{TF,j}, G::Frobenius{TG,j}) where {TF, TG, j}
    size(F) == size(G) || throw(DimensionMismatch())
    return Frobenius(j, F.c + G.c)
end
=#

function *(F::Frobenius, G::Frobenius)
    (n = size(F, 2)) == size(G, 2) || throw(DimensionMismatch(""))
    if F.colidx == G.colidx #Answer is still expressible as a Frobenius
        return Frobenius(F.colidx, F.c + G.c)
    else
        T = promote_type(eltype(F), eltype(G))
        M = Matrix{T}(F)
        M[G.colidx+1:end, G.colidx] = F.colidx < G.colidx ? G.c :
            Frobenius(F.colidx-G.colidx, F.c) * G.c
        return M
    end
end

inv(F::Frobenius) = Frobenius(F.colidx, -F.c)
