#The m-by-n Hilbert matrix has matrix elements
# H_{ij} = 1/(i+j-1)
export Hilbert, InverseHilbert
"""
[`Hilbert` matrix](http://en.wikipedia.org/wiki/Hilbert_matrix)

```julia
julia> A=Hilbert(5)
Hilbert{Rational{Int64}}(5,5)

julia> Matrix(A)
5x5 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9

julia> Matrix(Hilbert(5))
5x5 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9
```
Inverses are also integer matrices:

```julia
julia> inv(A)
5x5 Array{Rational{Int64},2}:
    25//1    -300//1     1050//1    -1400//1     630//1
  -300//1    4800//1   -18900//1    26880//1  -12600//1
  1050//1  -18900//1    79380//1  -117600//1   56700//1
 -1400//1   26880//1  -117600//1   179200//1  -88200//1
   630//1  -12600//1    56700//1   -88200//1   44100//1
```
"""
struct Hilbert{T} <: AbstractMatrix{T}
    m :: Int
    n :: Int
end
Hilbert(T::Type, m::Integer, n::Integer) = Hilbert{T}(m, n)
Hilbert(T::Type{<:Integer}, m::Integer, n::Integer) = Hilbert(Rational{T}, m, n)
Hilbert(m::Integer, n::Integer) = Hilbert(Int, m, n)
Hilbert(n::Integer) = Hilbert(n, n)

# Define its size
size(H::Hilbert, dim::Integer) = dim==1 ? H.m : dim==2 ? H.n : 1
size(H::Hilbert)= size(H,1), size(H,2)

# Index into a Hilbert matrix
function getindex(H::Hilbert{T}, i::Integer, j::Integer) where {T}
    return one(T)/(i+j-1)
end

# Dense version
Matrix(H::Hilbert) = [H[i,j] for i=1:size(H,1), j=1:size(H,2)]

# Some properties
ishermitian(H::Hilbert) = H.m==H.n
isposdef(H::Hilbert) = H.m==H.n
det(H::Hilbert) = inv(det(inv(H)))


# Inverse of a Hilbert matrix
struct InverseHilbert{T} <: AbstractMatrix{T}
    n :: Int
end
InverseHilbert(T::Type, n::Integer) = InverseHilbert{T}(n)
InverseHilbert(n::Integer) = InverseHilbert(Int, n)

# Define its size
size(A::InverseHilbert, dim::Integer) = dim==1 || dim==2 ? A.n : 1
size(A::InverseHilbert) = size(A,1), size(A,2)

# Index into a inverse Hilbert matrix
function getindex(A::InverseHilbert{T}, i::Integer, j::Integer) where {T}
    return T((-1)^(i+j)*(i+j-1)*binomial(A.n+i-1,A.n-j)*
              binomial(A.n+j-1,A.n-i)*binomial(i+j-2,i-1)^2)
end
# Use bigger numbers in case of bigger types
function getindex(A::InverseHilbert{T}, i::Integer, j::Integer) where {T<:Union{BigInt,BigFloat}}
    N = big(A.n)
    return T((-1)^(i+j)*(i+j-1)*binomial(N+i-1,N-j)*
              binomial(N+j-1,N-i)*binomial(big(i+j-2),i-1)^2)
end

# Explicit formula for the determinant
det(A::InverseHilbert{T}) where {T} = prod(T,(2k+1)*binomial(2k,k)^2 for k=1:A.n-1)

# Dense version
Matrix(A::InverseHilbert) = [A[i,j] for i=1:size(A,1), j=1:size(A,2)]

# Define Base.inv
function inv(H::Hilbert{T}) where {T}
    H.m == H.n || throw(ArgumentError("Works only for square Hilbert matrices."))
    return InverseHilbert(T,H.n)
end
function inv(A::InverseHilbert{T}) where {T}
    HT = promote_type(T,typeof(Rational(one(T))))
    return Hilbert(HT,A.n)
end
