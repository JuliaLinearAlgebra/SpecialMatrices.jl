export Circulant
"""
[`Circulant` matrix](https://en.wikipedia.org/wiki/Circulant_matrix)

```julia
julia> Circulant([1,2,3,4])
4x4 Circulant{Int64}:
 1  4  3  2
 2  1  4  3
 3  2  1  4
 4  3  2  1
```
"""
struct Circulant{T} <: AbstractArray{T, 2}
	c :: Vector{T}
end

getindex(C::Circulant, i::Int, j::Int) = C.c[mod(i-j,length(C.c))+1]
isassigned(C::Circulant, i::Int, j::Int) = isassigned(C.c,mod(i-j,length(C.c))+1)
size(C::Circulant, r::Int) = (r==1 || r==2) ? length(C.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(C::Circulant) = size(C,1), size(C,2)

# Fast matrix x vector via fft()
# see Golub, van Loan, Matrix Computations, John Hopkins, Baltimore, 1996, p. 202
function *(C::Circulant{T},x::Vector{T}) where T
    xt=fft(x)
    vt=fft(C.c)
    yt=vt.*xt
    typeof(x[1])==Int ? map(Int,round.(real(ifft(yt)))) : ( (T <: Real) ? map(T,real(ifft(yt))) : ifft(yt))
end

function mul!(y::StridedVector{T},C::Circulant{T},x::StridedVector{T}) where T
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

function Matrix(C::Circulant{T}) where T
	n=size(C, 1)
	M=Array{T}(undef,n,n)
	for i=1:n
		M[i:n,i] = C.c[1:n-i+1]
		M[1:i-1,i] = C.c[n-i+2:n]
	end
	M
end
