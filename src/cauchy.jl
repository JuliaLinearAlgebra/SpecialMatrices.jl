#### Cauchy matrix

export Cauchy
"""
[`Cauchy` matrix](http://en.wikipedia.org/wiki/Cauchy_matrix)
```julia
Cauchy(x,y)[i,j]=1/(x[i]+y[j])
Cauchy(x)=Cauchy(x,x)
cauchy(k::Number)=Cauchy(collect(1:k))

julia> Cauchy([1,2,3],[3,4,5])
3x3 Cauchy{Int64}:
 0.25      0.2       0.166667
 0.2       0.166667  0.142857
 0.166667  0.142857  0.125

julia> Cauchy([1,2,3])
3x3 Cauchy{Int64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667

julia> Cauchy(pi)
3x3 Cauchy{Float64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667
```
"""
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
size(A::Cauchy)= length(A.x), length(A.y)

# Index into a Cauchy
function getindex(A::Cauchy,i::Integer,j::Integer)
    return 1.0/(A.x[i]+A.y[j])
end # getindex

# Dense version of Cauchy
Matrix(A::Cauchy) = [A[i,j] for i=1:size(A,1), j=1:size(A,2)]
