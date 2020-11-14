#### Cauchy matrix

export Cauchy

struct Cauchy{T} <: AbstractMatrix{T}
    x::Vector{T} #
    y::Vector{T} #
end # immutable

function Cauchy(k::Number)
         Cauchy(collect(1:k),collect(1:k))
end

function Cauchy(x::Vector)
         Cauchy(x,x)
end

# Define its size

size(A::Cauchy, dim::Integer) = dim==1 ? length(A.x) : dim==2 ? length(A.y) : throw(ArgumentError("Invalid dimension $dim"))
size(A::Cauchy)= size(A,1), size(A,1)

# Index into a Cauchy
function getindex(A::Cauchy,i::Integer,j::Integer)
    return 1.0/(A.x[i]+A.y[j])
end # getindex

# Dense version of Cauchy
Matrix(A::Cauchy) = [A[i,j] for i=1:size(A,1), j=1:size(A,2)]
