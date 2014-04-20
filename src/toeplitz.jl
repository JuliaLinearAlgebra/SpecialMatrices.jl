export Circulant, Toeplitz

immutable Toeplitz{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

#XXX Inefficient but works
getindex(T::Toeplitz, i, j) = getindex(full(T), i, j)
isassigned(T::Toeplitz, i, j) = isassigned(full(T), i, j)

size(T::Toeplitz, r::Int) = (r==1 || r==2) ? 1 + div(length(T.c),2) :
    throw(ArgumentError("Invalid dimension $r"))
size(T::Toeplitz) = size(T,1), size(T,2)

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

#XXX Inefficient but works
getindex(C::Circulant, i, j) = getindex(full(C), i, j)
isassigned(C::Circulant, i, j) = isassigned(full(C), i, j)

size(C::Circulant, r::Int) = (r==1 || r==2) ? length(C.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(C::Circulant) = size(C,1), size(C,2)

function full{T}(C::Circulant{T})
	n=size(C, 1)
	M=Array(T, n, n)
	for i=1:n
		M[i:n,i] = C.c[1:n-i+1]
		M[1:i-1,i] = C.c[n-i+2:n]
	end
	M
end
