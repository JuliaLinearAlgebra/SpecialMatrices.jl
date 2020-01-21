#The m-by-n Hilbert matrix has matrix elements
# H_{ij} = 1/(i+j-1)
export Hilbert, InverseHilbert

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
