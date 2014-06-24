export Hankel

immutable Hankel{T} <: AbstractArray{T, 2}
    c :: Vector{T}
end
Hankel{T}(c::Vector{T}) = length(c) % 2 == 1 ? Hankel{T}(c) : throw(ArgumentError(""))

#XXX Inefficient but works
getindex(H::Hankel, i, j) = getindex(full(H), i, j)
isassigned(H::Hankel, i, j) = isassigned(full(H), i, j)

size(H::Hankel, r::Int) = (r==1 || r==2) ? 1 + div(length(H.c),2) :
    throw(ArgumentError("Invalid dimension $r"))
size(H::Hankel) = size(H,1), size(H,2)

function full{T}(H::Hankel{T})
    n=size(H, 1)
    M=Array(T, n, n)
    for i=1:n
        M[:,i] = H.c[i:i+n-1]
    end
    M
end
