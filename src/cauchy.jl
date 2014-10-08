#### Cauchy matrix

export Cauchy

immutable Cauchy{T} <: AbstractMatrix{T}
    x::Vector{T} #
    y::Vector{T} #
end # immutable

function Cauchy(k::Number)
         Cauchy([1:k],[1:k])
end

function Cauchy(x::Vector)
         Cauchy(x,x)
end

# Define its size

size(A::Cauchy, dim::Integer) = length(A.x)
size(A::Cauchy)= size(A,1), size(A,1)

# Index into a Cauchy
function getindex(A::Cauchy,i::Integer,j::Integer)
    return 1.0/(A.x[i]+A.y[j])
end # getindex

# Dense version of Cauchy
full(A::Cauchy) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]
