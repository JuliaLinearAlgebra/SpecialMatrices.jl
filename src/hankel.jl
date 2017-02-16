export Hankel

immutable Hankel{T<:Number} <: AbstractArray{T, 2}
    c::Vector{T}
    function (::Type{Hankel}){T}(c::Vector{T})
        if length(c) % 2 == 0
            throw(ArgumentError("Number of elements should be odd"))
        end
        new{T}(c)
    end
end

getindex(H::Hankel, i::Int, j::Int) = H.c[i+j-1]
isassigned(H::Hankel, i::Int, j::Int) = isassigned(H.c, i+j-1)

size(H::Hankel, r::Int) = (r==1 || r==2) ? 1 + div(length(H.c),2) :
    throw(ArgumentError("Invalid dimension $r"))
size(H::Hankel) = size(H,1), size(H,2)

# Fast matrix x vector via fft()
function *{T}(A::Hankel{T},x::Vector{T})
    Toeplitz(A.c)*reverse(x)
end

function A_mul_B!{T}(y::StridedVector{T},A::Hankel{T},x::StridedVector{T})
    xx=reverse(x)
    return A_mul_B!(y,Toeplitz(A.c),xx)
end

function full{T}(H::Hankel{T})
    n=size(H, 1)
    M=Array(T, n, n)
    for i=1:n
        M[:,i] = H.c[i:i+n-1]
    end
    M
end
