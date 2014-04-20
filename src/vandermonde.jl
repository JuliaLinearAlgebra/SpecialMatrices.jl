export Vandermonde

immutable Vandermonde{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

#XXX Inefficient but works
getindex(V::Vandermonde, i, j) = getindex(full(V), i, j)
isassigned(V::Vandermonde, i, j) = isassigned(full(V), i, j)

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
