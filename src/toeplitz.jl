export Circulant, Toeplitz, embed

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

function embed(To::Toeplitz{T}) where T
    return Circulant([To.c[div(end+1,2):end];To.c[1:div(end-1,2)]])
end
