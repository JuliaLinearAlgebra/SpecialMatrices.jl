export Hankel

struct Hankel{T} <: AbstractArray{T, 2}
    c :: Vector{T}
end
#Hankel{T}(c::Vector{T}) = length(c) % 2 == 1 ? Hankel{T}(c) : throw(ArgumentError(""))

getindex(H::Hankel, i::Int, j::Int) = H.c[i+j-1]
isassigned(H::Hankel, i::Int, j::Int) = isassigned(H.c, i+j-1)


size(H::Hankel, r::Int) = (r==1 || r==2) ? 1 + div(length(H.c),2) :
    throw(ArgumentError("Invalid dimension $r"))
size(H::Hankel) = size(H,1), size(H,2)

# Fast matrix x vector via fft()
function *(A::Hankel{T},x::Vector{T}) where T
    Toeplitz(A.c)*reverse(x)
end

function mul!(y::StridedVector{T},A::Hankel{T},x::StridedVector{T}) where T
    xx=reverse(x)
    return mul!(y,Toeplitz(A.c),xx)
end

function Matrix(H::Hankel{T}) where T
    n=size(H, 1)
    M=Array{T}(undef, n, n)
    for i=1:n
        M[:,i] = H.c[i:i+n-1]
    end
    M
end
