#### Cauchy matrix

export Cauchy

immutable Cauchy{T<:Number} <: AbstractMatrix{T}
    x::Vector{T}
    y::Vector{T}
end

function Cauchy(x, y)
    cx = collect(x)
    cy = collect(y)
    T = promote_type(eltype(cx), eltype(cy))
    vx = Vector{T}(cx)
    vy = Vector{T}(cy)
    Cauchy(vx, vy)
end
function Cauchy(x)
    vx = collect(x)
    Cauchy(vx, vx)
end
Cauchy(k::Number) =  Cauchy(1:k)

size(A::Cauchy) = (size(A.x,1), size(A.y,1))

function getindex(A::Cauchy,i::Integer,j::Integer)
    return 1.0/(A.x[i]+A.y[j])
end

full(A::Cauchy) = [A[i,j] for i=1:size(A,1), j=1:size(A,2)]
