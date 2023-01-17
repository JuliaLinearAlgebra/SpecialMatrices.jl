export Strang

using LinearAlgebra: LDLt, SymTridiagonal

"""
    Strang([T=Int,] n::Int)

Construct `Strang` matrix with elements of type `T`.
This special matrix named after Gilbert Strang
is symmetric, tridiagonal, and Toeplitz.
See Fig. 1 of [Julia paper](https://doi.org/10.1137/141000671).

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
    elseif abs(i-j) == 1
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


# LU factors for future reference:
#=
SparseArrays: spdiagm
   L = spdiagm(0 => ones(n), -1 => -(1:n-1) ./ (2:n))
   U = spdiagm(0 => (2:n+1) ./ (1:n), 1 => -ones(n-1))
=#

if VERSION >= v"1.8"
#=
Strang is a special case of SymTridiagonal,
and the default factorization of that class is LDLt
https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.factorize
so we use that here,
even though it does not fully exploit the sparse structure
of the LL' and LDL' factorizations of a Strang matrix.
=#
function LinearAlgebra.factorize(A::Strang)
    n = A.n
    ev = -(1:n-1) ./ (2:n) # L = spdiagm(0 => ones(n), -1 => ev)
    dv = (2:n+1) ./ (1:n) # D = Diagonal(dv)
    S = SymTridiagonal(dv, ev)
    return LDLt(S)
end
end
