export Strang

using LinearAlgebra: SymTridiagonal

"""
    Strang{T}([T=Float64,] n::Int)

Construct `Strang` matrix with elements of type `T`,
a special `SymTridiagonal` matrix named after Gilbert Strang.

```jldoctest
julia> Strang(6)
6×6 SymTridiagonal{Float64, Vector{Float64}}:
  2.0  -1.0    ⋅     ⋅     ⋅     ⋅ 
 -1.0   2.0  -1.0    ⋅     ⋅     ⋅ 
   ⋅   -1.0   2.0  -1.0    ⋅     ⋅ 
   ⋅     ⋅   -1.0   2.0  -1.0    ⋅ 
   ⋅     ⋅     ⋅   -1.0   2.0  -1.0
   ⋅     ⋅     ⋅     ⋅   -1.0   2.0

julia> Strang(Int16, 3)
3×3 SymTridiagonal{Int16, Vector{Int16}}:
  2  -1   0
 -1   2  -1
  0  -1   2
```
"""
function Strang(::Type{T}, n::Int) where {T <: Number}
    n > 0 || throw(ArgumentError("$n ≤ 0"))
    return SymTridiagonal(T(2)*ones(T, n), -ones(T, n-1))
end

Strang(n::Int) = Strang(Float64, n)
