export Circulant, Toeplitz

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



immutable Circulant{T} <: AbstractMatrix{T}
	c::Vector{T}
end

typealias Circ Circulant

getindex(C::Circ, i::Int, j::Int) = C.c[mod(i-j,length(C.c))+1]
isassigned(C::Circ, i::Int, j::Int) = isassigned(C.c,mod(i-j,length(C.c))+1)
size(C::Circ, r::Int) = (r==1 || r==2) ? length(C.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(C::Circ) = size(C,1), size(C,2)

# Generic code for multiplication by a Circulant matrix. Dispatch is used later
# on to ensure that the types returned are the correct ones. The idea is that
# if both inputs are Integer, then the output should also be integer.
# Similarly, if the inputs are both real, then the output should also be real.
function circ_dot(C::Circ, X::VecOrMat)
    if size(X, 1) != size(C.c, 1)
        throw(ArgumentError("C and X are not conformal."))
    end
    U = DFT(size(X, 1))
    return U * Diagonal(U * C.c) * U'X
end
function circ_dot(X::VecOrMat, C::Circ)
    if size(X, 2) != size(C.c, 1)
        throw(ArgumentError("C and X are not conformal."))
    end
    U = DFT(size(X, 1))
    return X * U * Diagonal(U * C.c) * U'
end

# In general, the return type of the multiplication operation should be a
# complex float.
*(C::Circ, X::VecOrMat) = circ_dot(C, X)
*(X::VecOrMat, C::Circ) = circ_dot(X, C)

# If both arguments contain reals, then the result should really contain reals.
*{T<:Real, V<:Real}(C::Circ{T}, X::VecOrMat{V}) = real(circ_dot(C, X))
*{T<:Real, V<:Real}(X::VecOrMat{V}, C::Circ{T}) = real(circ_dot(X, C))

# If both arguments contain (real) integers, then the result should be both
# Integer-valued and real. Promote to common type for the purposes of
# converting back to integer later.
*{T<:Integer, V<:Integer}(C::Circ{T}, X::VecOrMat{V}) = *(promote(C, X)...)
*{T<:Integer, V<:Integer}(X::VecOrMat{V}, C::Circ{T}) = *(promote(C, X)...)
*{T<:Integer}(C::Circ{T}, X::VecOrMat{T}) =
    convert(T, round.(real.(circ_dot(C, X))))
*{T<:Integer}(X::VecOrMat{T}, C::Circ{T}) =
    convert(T, round.(real.(circ_dot(C, X))))

# TODO: Create functionality to allow complex integers to be appropriately handled.




function A_mul_B!{T}(y::StridedVector{T},C::Circulant{T},x::StridedVector{T})
    xt=fft(x)
    vt=fft(C.c)
    yt=ifft(vt.*xt)
    if T<: Int
        map!(round,y,yt) 
    elseif T<: Real
        map!(real,y,yt)
    else
        copy!(y,yt)
    end
    return y
end

function full{T}(C::Circulant{T})
	n=size(C, 1)
	M=Array(T, n, n)
	for i=1:n
		M[i:n,i] = C.c[1:n-i+1]
		M[1:i-1,i] = C.c[n-i+2:n]
	end
	M
end
