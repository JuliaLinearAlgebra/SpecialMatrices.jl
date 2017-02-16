#### Cauchy matrix

export Cauchy

immutable Cauchy{T<:Number} <: AbstractMatrix{T}
    x::Vector{T}
    y::Vector{T}
end

function Cauchy(x, y)
    T = promote_type(eltype(x), eltype(y))
    vx = Vector{T}(x)
    vy = Vector{T}(y)
    Cauchy(vx, vy)
end
function Cauchy(x)
    vx = Vector(x)
    Cauchy(vx, vx)
end
Cauchy(k::Number) =  Cauchy(1:k)

size(A::Cauchy) = (size(A.x,1), size(A.y,1))

function getindex(A::Cauchy,i::Integer,j::Integer)
    return 1.0/(A.x[i]+A.y[j])
end

full(A::Cauchy) = [A[i,j] for i=1:size(A,1), j=1:size(A,2)]
