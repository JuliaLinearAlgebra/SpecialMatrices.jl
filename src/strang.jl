export Strang

"""
    Strang([T=Int,] n::Int)

Construct `Strang` matrix with elements of type `T`.
This special matrix named after Gilbert Strang
is symmetric, tridiagonal, and Toeplitz.
See Fig. 1 of [Julia paper](http://doi.org/10.1137/141000671).

```jldoctest
julia> Strang(4)
4×4 Strang{Int64}:
  2  -1   0   0
 -1   2  -1   0
  0  -1   2  -1
  0   0  -1   2

julia> Strang(Int16, 3)
3×3 Strang{Int16}:
  2  -1   0
 -1   2  -1
  0  -1   2
```
"""
struct Strang{T <: Number} <: AbstractMatrix{T}
    n::Int

    function Strang(::Type{T}, n::Int) where {T <: Number}
        n > 0 || throw(ArgumentError("$n ≤ 0"))
        return new{T}(n)
     end
end

Strang(n::Int) = Strang(Int, n)

# Properties
size(s::Strang) = (s.n, s.n)
LinearAlgebra.adjoint(s::Strang) = s
LinearAlgebra.transpose(s::Strang) = s
LinearAlgebra.issymmetric(s::Strang) = true
LinearAlgebra.ishermitian(s::Strang) = true
LinearAlgebra.isposdef(s::Strang) = true

@inline Base.@propagate_inbounds function getindex(A::Strang{T}, i::Integer, j::Integer) where T
    @boundscheck checkbounds(A, i, j)
    if i == j
        return T(2)
    elseif i == j + 1 || i == j - 1
        return -one(T)
    else
        return zero(T)
    end
end

@inline Base.@propagate_inbounds function *(A::Strang{T}, x::AbstractVector) where T
    @boundscheck length(x) == A.n || throw(ArgumentError("length(x) not n=$n"))
    Base.require_one_based_indexing(x)
    y = collect(T(2) * x) # from diagonal
    y[1:end-1] .-= @view x[2:end]
    y[2:end] .-= @view x[1:end-1]
    return y
end
