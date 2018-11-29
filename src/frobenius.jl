export Frobenius

#Frobenius matrices or Gaussian elimination matrices
#
#Matrices of the form
#[ 1 0 ...     0 ]
#[ 0 1 ...     0 ]
#[ .........     ]
#[ ... 1 ...     ]
#[ ... c1 1 ...  ]
#[ ... c2 0 1 ...]
#[ ............. ]
#[ ... ck ...   1]
#
#i.e. an identity matrix with nonzero subdiagonal elements along a single
#column.
#
#In this implementation, the subdiagonal of the nonzero column is stored as a
#dense vector, so that the size can be inferred automatically as j+k where j is
#the index of the column and k is the number of subdiagonal elements.

struct Frobenius{T} <: AbstractArray{T, 2}
    colidx :: Int
    c :: Vector{T}
end

#Basic property computations
size(F::Frobenius, r::Int) = (r==1 || r==2) ? F.colidx + length(F.c) : 
    throw(ArgumentError("Frobenius matrix is of rank 2"))

function size(F::Frobenius)
    n = F.colidx + length(F.c)
    n, n
end

#XXX Inefficient but works
getindex(F::Frobenius, i, j) = getindex(Matrix(F), i, j)
isassigned(F::Frobenius, i::Integer, j::Integer) = isassigned(Matrix(F), i, j)

function Matrix(F::Frobenius{T}) where T
    M = Matrix{T}(I, size(F))
    M[F.colidx+1:end, F.colidx] = F.c
    M
end

#Linear algebra stuff
function mul!(F::Frobenius{T}, b::Vector{T}) where T
    (n = size(F, 2)) == length(b) || throw(DimensionMismatch("$n $(length(b))"))
    for i=F.colidx+1:n
        b[i] += F.c[i-F.colidx] * b[F.colidx]
    end
    b
end
*(F::Frobenius{T}, b::Vector{T}) where T = mul!(F, copy(b))

function *(F::Frobenius{T}, G::Frobenius{T}) where T
    (n = size(F, 2)) == size(G, 2) || throw(DimensionMismatch(""))
    if F.colidx == G.colidx #Answer is still expressible as a Frobenius
        return Frobenius(F.colidx, F.c + G.c)
    else
        M = Matrix(F)
        M[G.colidx+1:end, G.colidx] = F.colidx < G.colidx ? G.c :
            Frobenius(F.colidx-G.colidx, F.c) * G.c
        return M
    end
end

inv(F::Frobenius) = Frobenius(F.colidx, -F.c)
