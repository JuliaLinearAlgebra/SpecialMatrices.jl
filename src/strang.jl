export Strang
"""
`Strang` matrix

A special `SymTridiagonal` matrix named after Gilbert Strang

```julia
julia> Strang(6)
6x6 Strang{Float64}:
  2.0  -1.0   0.0   0.0   0.0   0.0
 -1.0   2.0  -1.0   0.0   0.0   0.0
  0.0  -1.0   2.0  -1.0   0.0   0.0
  0.0   0.0  -1.0   2.0  -1.0   0.0
  0.0   0.0   0.0  -1.0   2.0  -1.0
  0.0   0.0   0.0   0.0  -1.0   2.0
```
"""
struct Strang{T} <: AbstractArray{T, 2}
    n :: Int
end
Strang(n::Int) = Strang{Float64}(n)
strang(T, n)= n >1 ? SymTridiagonal(2ones(T, n),-ones(T, n-1)) :
              n==1 ? Diagonal([2one(T)]) : error("Invalid dimension ", n)
function getindex(S::Strang{T}, i, j) where T
    i == j && return 2
    abs(i - j) == 1 && return -1
    0
end
getindex(S::Strang{T}, I...) where T = getindex(S,ind2sub(size(S),I...)...)
size(S::Strang, r::Int) = r==1 || r==2 ? S.n : throw(ArgumentError("Invalid dimension $r"))
size(S::Strang) = S.n, S.n
Matrix(S::Strang{T}) where T = Matrix(strang(T, S.n))
*(A::VecOrMat{T}, S::Strang{T}) where T = A*Matrix(S)
*(S::Strang{T}, A::VecOrMat{T}) where T = Matrix(S)*A
