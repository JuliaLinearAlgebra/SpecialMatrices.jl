#### Riemann Matrix
import Base: size
export Riemann
"""
`Riemann` matrix

Riemann matrix is defined as `A = B[2:N+1, 2:N+1]`, where
`B[i,j] = i-1` if `i` divides `j`, and `-1` otherwise.
[Riemann hypothesis](http://en.wikipedia.org/wiki/Riemann_hypothesis) holds
if and only if `det(A) = O( N! N^(-1/2+epsilon))` for every `epsilon > 0`.

```julia
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

# Define its size

size(A::Riemann, dim::Integer) = A.n
size(A::Riemann)= size(A,1), size(A,1)

# Index into a Riemann
function getindex(A::Riemann,i::Integer,j::Integer)
#    return (i+1)%(j+1)==0 ? i : -1
    if i<=j && (j+1)%(i+1)==0
        return i
    else
        return -1
    end
end # getindex

# Dense version of Riemann
Matrix(A::Riemann) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]
