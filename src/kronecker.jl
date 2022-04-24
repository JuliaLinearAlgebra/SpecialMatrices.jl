export Kronecker


"""
    Kronecker(A::AbstractVecOrMat,B::AbstractVecOrMat)

Create a "lazy" M = A ⊗ B
[`Kronecker` product matrix](https://en.wikipedia.org/wiki/Kronecker_product)

The `transpose` and `adjoint` operations are also lazy.

```jldoctest kron
julia> A = Kronecker([1 2; 3 4],[5 6; 7 8])
4×4 Kronecker{Int64}:
  5   6  10  12
  7   8  14  16
 15  18  20  24
 21  24  28  32

julia> A'
4×4 Kronecker{Int64}:
  5   7  15  21
  6   8  18  24
 10  14  20  28
 12  16  24  32
```

The properties of the Kronecker product are utilized to efficiently
compute without forming it.

```jldoctest kron
julia> n = 50; A,B = rand(n,n), rand(n,n); K = Kronecker(A,B); M = kron(A,B);

julia> @time inv(K) ≈ @time inv(M)
  1.142912 seconds (6 allocations: 48.924 MiB)
 42.812522 seconds (411.07 M allocations: 11.528 GiB, 3.72% gc time)
true

```
"""
struct Kronecker{T} <: AbstractMatrix{T}
  A :: AbstractMatrix{T}
  B :: AbstractMatrix{T}
end

"""
    krontype(T1::Type, T2::Type)
Determine the return type of `Kronecker(A::T1,B::T2)`.
"""
function krontype(Ta::Type, Tb::Type)
    T = typeof(oneunit(Ta) * oneunit(Tb))
    return T
end

# function Kronecker(A::AbstractMatrix{Ta},B::AbstractMatrix{Tb}) where {Ta,Tb}
#   T = krontype(Ta::Type, Tb::Type)
#   return Kronecker{T}(convert(AbstractMatrix{T}, A),convert(AbstractMatrix{T}, B))
# end

function Kronecker(A,B)
  T = krontype(eltype(A), eltype(B))
  return Kronecker{T}(A,B)
end

@inline Base.@propagate_inbounds function getindex(M::Kronecker, i::Int, j::Int,)
  @boundscheck checkbounds(M, i, j)

  mB,nB = size(M.B)

  qi,ri = divrem(i - 1, mB)
  qj,rj = divrem(j - 1, nB)

  return (@inbounds M.A[qi + 1, qj + 1] * M.B[ri + 1, rj + 1])
end

size(M::Kronecker) = size(M.A) .* size(M.B)

# y = (A \otimes B) * V
function *(M::Kronecker{T1}, V::AbstractMatrix{T2}) where {T1, T2}

  T = krontype(T1,T2)

  mA,nA = size(M.A)
  mB,nB = size(M.B)

  mV,nV = size(V)

  y = Array{T}(undef, mA*mB, nV)

  for j=1:nV
    y[:,j] .= reshape(M.B * reshape(V[:,j], nB, mV ÷ nB) * transpose(M.A), mA*mB,1)
  end

  return y
end

# y = (A \otimes B) * v
function *(M::Kronecker{T1}, v::AbstractVector{T2}) where {T1, T2}

  T = krontype(T1,T2)

  mA,nA = size(M.A)
  mB,nB = size(M.B)

  mv = length(v)

  y = reshape(M.B * reshape(v, nB, mv ÷ nB) * transpose(M.A), mA*mB,1)

  return y
end

# y = (A \otimes B) \ v
function \(M::Kronecker{T1}, v::AbstractVector{T2}) where {T1, T2}

  T = krontype(T1,T2)

  mA,nA = size(M.A)
  mB,nB = size(M.B)

  mv = length(v)

  y = reshape(M.B \ reshape(v, nB, mv ÷ nB) / transpose(M.A), mA*mB,1)

  return y
end

# y = (A \otimes B) \ v
function \(M::Kronecker{T1}, V::AbstractMatrix{T2}) where {T1, T2}

  T = krontype(T1,T2)

  mA,nA = size(M.A)
  mB,nB = size(M.B)

  mV,nV = size(V)

  Y = Array{T}(undef, mA*mB, nV)

  for j=1:nV
    Y[:,j] .= reshape(M.B \ reshape(V[:,j], nB, mV ÷ nB) / transpose(M.A), mA*mB,1)
  end

  return Y
end

function inv(M::Kronecker)
  Kronecker(inv(M.A),inv(M.B))
end

transpose(M::Kronecker) = Kronecker(transpose(M.A),transpose(M.B))

adjoint(M::Kronecker) = Kronecker(adjoint(M.A),adjoint(M.B))

tr(M::Kronecker) = tr(M.A) * tr(M.B)

zlt(c,z) = real(c) == real(z) ? imag(c)<imag(z) : real(c) < real(z)

eigvals(M::Kronecker) = sort(kron(eigvals(M.A), eigvals(M.B)); lt=zlt)
