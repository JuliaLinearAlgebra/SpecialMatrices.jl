#### Cauchy matrix

export Cauchy

"""
    [`Cauchy` matrix](http://en.wikipedia.org/wiki/Cauchy_matrix)

* `Cauchy(x,y)[i,j] = 1/(x[i]+y[j])`
* `Cauchy(x) = Cauchy(x,x)`
* `Cauchy(k::Int) = Cauchy(1:k)`

```julia
julia> Cauchy(1:3, 3:5)
3x3 Cauchy{Int64}:
 0.25      0.2       0.166667
 0.2       0.166667  0.142857
 0.166667  0.142857  0.125

julia> Cauchy(1:3)
3x3 Cauchy{Int64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667

julia> Cauchy(3)
3x3 Cauchy{Float64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667
```
"""
struct Cauchy{T} <: AbstractMatrix{T}
    x::Vector{T}
    y::Vector{T}
end # immutable

function Cauchy(x::AbstractVector{Tx}, y::AbstractVector{Ty}) where {Tx <: Number, Ty <: Number}
    T = promote_type(Tx, Ty)
    T = eltype(1f0 * one(T)) # ensure at least Float32
    vx = Vector{T}(x)
    vy = Vector{T}(y)
    Cauchy(vx, vy)
end

function Cauchy(x)
    vx = collect(x)
    Cauchy(vx, vx)
end

Cauchy(k::Int) = Cauchy(Rational{Int}.(1:k))

size(A::Cauchy) = (size(A.x,1), size(A.y,1))

function getindex(A::Cauchy{T}, i::Integer, j::Integer) where {T}
    @boundscheck checkbounds(A, i, j)
    @inbounds return 1 / (A.x[i] + A.y[j])
end
