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

function full{T}(To::Toeplitz{T})
	n=size(To, 1)
	M=Array(T, n, n)
	for i=1:n
		M[i:n,i] = To.c[n:2n-i]
		M[1:i-1,i] = To.c[n-i+1:n-1]
	end
	M
end



immutable Circulant{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

getindex(C::Circulant, i::Int, j::Int) = C.c[mod(i-j,length(C.c))+1]
isassigned(C::Circulant, i::Int, j::Int) = isassigned(C.c,mod(i-j,length(C.c))+1)
size(C::Circulant, r::Int) = (r==1 || r==2) ? length(C.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(C::Circulant) = size(C,1), size(C,2)

# Fast matrix x vector via fft()
# see Golub, van Loan, Matrix Computations, John Hopkins, Baltimore, 1996, p. 202 
function *{T}(C::Circulant{T},x::Vector{T})
    xt=fft(x)
    vt=fft(C.c)
    yt=vt.*xt
    typeof(x[1])==Int ? map(Int,round(real(ifft(yt)))): ( (T <: Real) ? map(T,real(ifft(yt))) : ifft(yt))
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
