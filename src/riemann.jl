#### Riemann Matrix
export Riemann

import Base: size, getindex

"""
    Riemann(N::Int)

Construct `N × N` `Riemann` matrix,
defined as `A = B[2:N+1, 2:N+1]`, where
`B[i,j] = i-1` if `i` divides `j`, and `-1` otherwise.
The
[Riemann hypothesis](http://en.wikipedia.org/wiki/Riemann_hypothesis) holds
if and only if `det(A) = O( N! N^(-1/2+epsilon))` for every `epsilon > 0`.

```jldoctest
julia> Riemann(7)
7x7 Riemann{Int64}:
  1  -1   1  -1   1  -1   1
 -1   2  -1  -1   2  -1  -1
 -1  -1   3  -1  -1  -1   3
 -1  -1  -1   4  -1  -1  -1
 -1  -1  -1  -1   5  -1  -1
 -1  -1  -1  -1  -1   6  -1
 -1  -1  -1  -1  -1  -1   7
```

For more details see Friedrich Roesler, Riemann's hypothesis as an eigenvalue problem,
Linear Algebra and its Applications, Vol. 81, (1986)
http://www.sciencedirect.com/science/article/pii/0024379586902557
"""
struct Riemann{Int} <: AbstractMatrix{Int}
    n::Int
end # immutable

size(A::Riemann) = (A.n, A.n)

# Index into a Riemann
@inline Base.@propagate_inbounds function getindex(A::Riemann, i::Integer, j::Integer)
#    return (i+1)%(j+1)==0 ? i : -1
    @boundscheck checkbounds(A, i, j)
    if i ≤ j && (j+1) % (i+1) == 0
        return i
    else
        return -1
    end
end
