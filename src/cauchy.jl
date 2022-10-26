#### Cauchy matrix

export Cauchy

"""
    [`Cauchy` matrix](http://en.wikipedia.org/wiki/Cauchy_matrix)

* `Cauchy(x,y)[i,j] = 1 / (x[i] + y[j])`
* `Cauchy(x) = Cauchy(x,x)`
* `Cauchy(k::Int) = Cauchy(1:k)`

Both `x` and `y` can be any iterable (typically vectors),
but all elements of `x` must have the same type; likewise for `y`.

```julia
julia> Cauchy([2.0 1], (0, 1, 2))
2×3 Cauchy{Float64, Matrix{Float64}, Tuple{Int64, Int64, Int64}}:
 0.5  0.333333  0.25
 1.0  0.5       0.333333

julia> Cauchy(1:3)
3×3 Cauchy{Rational{Int64}, UnitRange{Int64}, UnitRange{Int64}}:
 1//2  1//3  1//4
 1//3  1//4  1//5
 1//4  1//5  1//6

julia> Cauchy(3)
3×3 Cauchy{Rational{Int64}, UnitRange{Int64}, UnitRange{Int64}}:
 1//2  1//3  1//4
 1//3  1//4  1//5
 1//4  1//5  1//6
```
"""
struct Cauchy{T,X,Y} <: AbstractMatrix{T}
    x::X
    y::Y

    # The element type T corresponds to the reciprocal
    # of the sum of pairs of elements of x and y
    function Cauchy(x::X, y::Y) where {X, Y}
        Tx = eltype(first(x))
        Ty = eltype(first(y))
        all(==(Tx), eltype.(x)) || throw(ArgumentError("inconsistent x element types"))
        all(==(Ty), eltype.(y)) || throw(ArgumentError("inconsistent y element types"))
        T = eltype(one(Tx) + one(Ty)) 
        T = T <: Integer ? Rational{T} : eltype(1 / one(T))
        return new{T,X,Y}(x, y)
    end
end

Cauchy(x) = Cauchy(x, x)
Cauchy(k::Integer) = Cauchy(1:k)

size(A::Cauchy) = (length(A.x), length(A.y))

@inline Base.@propagate_inbounds function getindex(
    A::Cauchy{<:Rational},
    i::Integer,
    j::Integer,
)
    @boundscheck checkbounds(A, i, j)
    @inbounds return 1 // (A.x[i] + A.y[j])
end

@inline Base.@propagate_inbounds function getindex(
    A::Cauchy,
    i::Integer,
    j::Integer,
)
    @boundscheck checkbounds(A, i, j)
    @inbounds return 1 / (A.x[i] + A.y[j])
end
