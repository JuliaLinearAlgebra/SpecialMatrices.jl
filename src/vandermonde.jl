export Vandermonde

immutable Vandermonde{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

getindex(V::Vandermonde, i::Int, j::Int) = V.c[i]^(j-1)
isassigned(V::Vandermonde, i::Int, j::Int) = isassigned(V.c, i)

size(V::Vandermonde, r::Int) = (r==1 || r==2) ? length(V.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(V::Vandermonde) = length(V.c), length(V.c)

function full{T}(V::Vandermonde{T})
	n=size(V, 1)
	M=Array(T, n, n)
	for i=1:n
		M[:,i] = V.c.^(i-1)
	end
	M
end
