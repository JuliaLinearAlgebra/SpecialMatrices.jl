import Base: A_mul_B!, full, getindex, inv, isassigned, size, *

export FrobeniusMatrix

#Frobenius matrices or Gaussian elimination matrices
#
#Matrices of the form
#[ 1 0 ...     0 ]
#[ 0 1 ...     0 ]
#[ .........     ]
#[ ... 1 ...     ]
#[ ... z1 1 ...  ]
#[ ... z2 0 1 ...]
#[ ............. ]
#[ ... zk ...   1]
#
#i.e. an identity matrix with nonzero subdiagonal elements along a single
#column.
#
#In this implementation, the subdiagonal of the nonzero column is stored as a
#dense vector, so that the size can be inferred automatically as j+k where j is
#the index of the column and k is the number of subdiagonal elements.

type FrobeniusMatrix{T} <: AbstractArray{T, 2}
    colidx :: Int
    z ::Vector{T}
end

#Basic property computations
size(F::FrobeniusMatrix, r::Int) = (r==1 || r==2) ? F.colidx + length(F.z) : 
    throw(ArgumentError("Frobenius matrix is of rank 2"))

function size(F::FrobeniusMatrix)
    n = F.colidx + length(F.z)
    n, n
end

#XXX Inefficient but works
getindex(F::FrobeniusMatrix, i, j) = getindex(full(F), i, j)
isassigned(F::FrobeniusMatrix, i, j) = isassigned(full(F), i, j)

function full{T}(F::FrobeniusMatrix{T})
    M = eye(T, size(F, 1))
    M[F.colidx+1:end, F.colidx] = F.z
    M
end

#Linear algebra stuff
function A_mul_B!{T}(F::FrobeniusMatrix{T}, b::Vector{T})
    (n = size(F, 2)) == length(b) || throw(DimensionMismatch("$n $(length(b))"))
    for i=F.colidx+1:n
        b[i] += F.z[i-F.colidx] * b[F.colidx]
    end
    b
end
*{T}(F::FrobeniusMatrix{T}, b::Vector{T}) = A_mul_B!(F, copy(b))

function *{T}(F::FrobeniusMatrix{T}, G::FrobeniusMatrix{T})
    (n = size(F, 2)) == size(G, 2) || throw(DimensionMismatch(""))
    if F.colidx == G.colidx #Answer is still expressible as a FrobeniusMatrix
        return FrobeniusMatrix(F.colidx, F.z + G.z)
    else
        M⧀, M⧁ = F.colidx > G.colidx ? (G, F) : (F, G)
        M = full(M⧁)
        M[M⧀.colidx+1:end, M⧀.colidx] = FrobeniusMatrix(M⧁.colidx-M⧀.colidx, F.z) * G.z
        return M
    end
end

inv(F::FrobeniusMatrix) = FrobeniusMatrix(F.colidx, -F.z)

