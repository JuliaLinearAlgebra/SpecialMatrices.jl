#The m-by-n Hilbert matrix has matrix elements
# H_{ij} = 1/(i+j-1)
#
export Hilbert, InverseHilbert

"""
   Hilbert(m [,n])

Construct `m × m` or `m × n`
[`Hilbert` matrix](https://en.wikipedia.org/wiki/Hilbert_matrix)
from its specified dimensions,
where element `i,j` equal to `1 / (i+j-1)`.

```jldoctest hilbert1
julia> A = Hilbert(5)
5×5 Hilbert{Rational{Int64}}:
  1    1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9

julia> Matrix(A)
5×5 Matrix{Rational{Int64}}:
  1    1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9
```

Inverses are also integer matrices:

```jldoctest hilbert1
julia> inv(A)
5×5 InverseHilbert{Rational{Int64}}:
    25    -300     1050    -1400     630
  -300    4800   -18900    26880  -12600
  1050  -18900    79380  -117600   56700
 -1400   26880  -117600   179200  -88200
   630  -12600    56700   -88200   44100
```
"""
struct Hilbert{T} <: AbstractMatrix{T}
    m :: Int
    n :: Int
end

Hilbert(::Type{T}, m::Integer, n::Integer) where {T <: Number} = Hilbert{T}(m, n)
Hilbert(::Type{T}, m::Integer, n::Integer) where {T <: Integer} = Hilbert(Rational{T}, m, n)
Hilbert(m::Integer, n::Integer) = Hilbert(Int, m, n)
Hilbert(n::Integer) = Hilbert(n, n)
Hilbert(::Type{T}, n::Integer) where {T <: Number} = Hilbert(T, n, n)

# Define its size
size(H::Hilbert, dim::Integer) = dim==1 ? H.m : dim==2 ? H.n : 1
size(H::Hilbert) = size(H,1), size(H,2)

# Index into a Hilbert matrix
function getindex(H::Hilbert{T}, i::Integer, j::Integer) where {T}
    return one(T)/(i+j-1)
end

# Dense version (provided in Core)
#if VERSION < v"1.6"
#    Matrix(H::Hilbert) = [H[i,j] for i in 1:size(H,1), j in 1:size(H,2)]
#end

# Some properties
ishermitian(H::Hilbert) = H.m==H.n
isposdef(H::Hilbert) = H.m==H.n
det(H::Hilbert) = inv(det(inv(H)))


# Inverse of a Hilbert matrix
struct InverseHilbert{T} <: AbstractMatrix{T}
    n :: Int
end

InverseHilbert(::Type{T}, n::Integer) where {T <: Number} = InverseHilbert{T}(n)
InverseHilbert(n::Integer) = InverseHilbert(Int, n)

# Define its size
size(A::InverseHilbert, dim::Integer) = dim==1 || dim==2 ? A.n : 1
size(A::InverseHilbert) = size(A,1), size(A,2)

# Properties
ishermitian(::InverseHilbert) = true
isposdef(::InverseHilbert) = true

# Index into a inverse Hilbert matrix
function getindex(A::InverseHilbert{T}, i::Integer, j::Integer) where {T}
    out = (-1)^(i+j) * (i+j-1) * binomial(A.n+i-1, A.n-j) *
        binomial(A.n+j-1, A.n-i) * binomial(i+j-2, i-1)^2
    return T(out)
end

# Use bigger numbers in case of bigger types
function getindex(A::InverseHilbert{T}, i::Integer, j::Integer) where {T<:Union{BigInt,BigFloat}}
    N = big(A.n)
    out = (-1)^(i+j) * (i+j-1) * binomial(N+i-1, N-j) *
        binomial(N+j-1, N-i) * binomial(big(i+j-2), i-1)^2
    return T(out)
end

"""
    det(A::InverseHilbert)
Explicit formula for the determinant.
Caution: this function overflows easily.
"""
det(A::InverseHilbert{T}) where {T} = prod(T, (2k+1) * binomial(2k,k)^2 for k in 1:A.n-1)

# Dense version (in Core)
#if VERSION < v"1.6"
#    Matrix(A::InverseHilbert) = [A[i,j] for i in 1:size(A,1), j in 1:size(A,2)]
#end

# Define Base.inv
function inv(H::Hilbert{T}) where {T}
    H.m == H.n || throw(ArgumentError("Works only for square Hilbert matrices."))
    return InverseHilbert(T, H.n)
end

function inv(A::InverseHilbert{T}) where {T}
    HT = promote_type(T,typeof(Rational(one(T))))
    return Hilbert(HT, A.n)
end
