export Toeplitz

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
