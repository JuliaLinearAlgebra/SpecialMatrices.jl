export Toeplitz, embed
"""
[`Toeplitz` matrix](http://en.wikipedia.org/wiki/Toeplitz_matrix)

Input is a vector of odd length.

```julia
julia> Toeplitz(collect(-4:4))
5x5 Toeplitz{Int64}:
 0  -1  -2  -3  -4
 1   0  -1  -2  -3
 2   1   0  -1  -2
 3   2   1   0  -1
 4   3   2   1   0
```
"""
struct Toeplitz{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

getindex(T::Toeplitz, i::Int, j::Int) = T.c[i-j+div(length(T.c)+1,2)]
isassigned(T::Toeplitz, i::Int, j::Int) = isassigned(T.c,i-j+div(length(T.c)+1,2))

size(T::Toeplitz, r::Int) = (r==1 || r==2) ? 1 + div(length(T.c),2) :
    throw(ArgumentError("Invalid dimension $r"))
size(T::Toeplitz) = size(T,1), size(T,2)

# Fast matrix x vector multiplication via embedding Toeplitz() into Circulant()
function *(A::Toeplitz{T},x::Vector{T}) where T
    n=length(A.c)
    k=div(n+1,2)
    C=Circulant([A.c[k:n];A.c[1:k-1]])
    (C*[x;zeros(T,k-1)])[1:k]
end

function mul!(y::StridedVector{T},A::Toeplitz{T},x::StridedVector{T}) where T
    n=length(A.c)
    k=div(n+1,2)
    C=Circulant([A.c[k:n];A.c[1:k-1]])
    xx=[x;zeros(T,k-1)]
    yy=mul!(similar(xx),C,xx)
    copyto!(y, 1, yy, 1, length(y))
    return y
end

function Matrix(To::Toeplitz{T}) where T
	n=size(To, 1)
	M=Array{T}(undef,n,n)
	for i=1:n
		M[i:n,i] = To.c[n:2n-i]
		M[1:i-1,i] = To.c[n-i+1:n-1]
	end
	M
end

function embed(To::Toeplitz{T}) where T
    return Circulant([To.c[div(end+1,2):end];To.c[1:div(end-1,2)]])
end
