export Toeplitz

immutable Toeplitz{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

getindex(T::Toeplitz, i::Int, j::Int) = T.c[i-j+div(length(T.c)+1,2)]
isassigned(T::Toeplitz, i::Int, j::Int) = isassigned(T.c,i-j+div(length(T.c)+1,2))

size(T::Toeplitz, r::Int) = (r==1 || r==2) ? 1 + div(length(T.c),2) :
    throw(ArgumentError("Invalid dimension $r"))
size(T::Toeplitz) = size(T,1), size(T,2)

# Fast matrix x vector multiplication via embedding Toeplitz() into Circulant()
function *{T}(A::Toeplitz{T},x::Vector{T})
    n=length(A.c)
    k=Int(round((n+1)/2))
    C=Circulant([A.c[k:n];A.c[1:k-1]])
    (C*[x;zeros(T,k-1)])[1:k]
end

function A_mul_B!{T}(y::StridedVector{T},A::Toeplitz{T},x::StridedVector{T})
    n=length(A.c)
    k=Int(round((n+1)/2))
    C=Circulant([A.c[k:n];A.c[1:k-1]])
    xx=[x;zeros(T,k-1)]
    yy=A_mul_B!(similar(xx),C,xx)
    copy!(y, 1, yy, 1, length(y))
    return y
end

function full{T}(To::Toeplitz{T})
	n=size(To, 1)
	M=Array(T, n, n)
	for i=1:n
		M[i:n,i] = To.c[n:2n-i]
		M[1:i-1,i] = To.c[n-i+1:n-1]
	end
	M
end
