export Vandermonde

struct Vandermonde{T} <: AbstractArray{T, 2}
    c :: Vector{T}
end

getindex(V::Vandermonde, i::Int, j::Int) = V.c[i]^(j-1)
isassigned(V::Vandermonde, i::Int, j::Int) = isassigned(V.c, i)

size(V::Vandermonde, r::Int) = (r==1 || r==2) ? length(V.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(V::Vandermonde) = length(V.c), length(V.c)

function Matrix(V::Vandermonde{T}) where T
    n=size(V, 1)
    M=Array{T}(undef, n, n)
    M[:,1] .= 1
    for j=2:n
        for i=1:n
	    M[i,j] = M[i,j-1]*V.c[i]
        end
    end
    M
end

function Matrix(V::Adjoint{T,Vandermonde{T}}) where T
    n=size(V, 1)
    M=Array{T}(undef, n, n)
    M[1,:] .= 1
    for j=1:n
        for i=2:n
	    M[i,j] = M[i-1,j]*adjoint(V.parent.c[j])
        end
    end
    M
end

function Matrix(V::Transpose{T,Vandermonde{T}}) where T
    n=size(V, 1)
    M=Array{T}(undef, n, n)
    M[1,:] .= 1
    for j=1:n
        for i=2:n
	    M[i,j] = M[i-1,j]*V.parent.c[j]
        end
    end
    M
end

function \(V::Adjoint{T1,Vandermonde{T1}}, y::AbstractVecOrMat{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, size(y))
    copyto!(x, y)
    pvand!(adjoint(V.parent.c), x)
    return x
end

function \(V::Transpose{T1,Vandermonde{T1}}, y::AbstractVecOrMat{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, size(y))
    copyto!(x, y)
    pvand!(V.parent.c, x)
    return x
end

function \(V::Vandermonde{T1}, y::AbstractVecOrMat{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, size(y))
    copyto!(x, y)
    dvand!(V.c, x)
    return x
end


function vandtype(T1::Type, T2::Type)
    # Figure out the return type of Vandermonde{T1} \ Vector{T2}
    T = promote_type(T1, T2)
    S = typeof(oneunit(T)/oneunit(T1))
    return S
end

"""
    pvand!(a, b) -> b

Solves system ``A^T*x = b`` in-place.

``A^T`` is transpose of Vandermonde matrix ``A_{ij} = a_i^{j-1}``.

Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function pvand!(alpha, B)
    n = length(alpha);
    if n != size(B,1)
        throw(DimensionMismatch("matrix has dimensions ($n,$n) but right hand side has $(size(B,1)) rows"))
    end
    nrhs = size(B,2)
    @inbounds begin
        for j=1:nrhs
            for k=1:n-1
                for i=n:-1:k+1
                    B[i,j] = B[i,j]-alpha[k]*B[i-1,j]
                end
            end
            for k=n-1:-1:1
                for i=k+1:n
                    B[i,j] = B[i,j]/(alpha[i]-alpha[i-k])
                end
                for i=k:n-1
                    B[i,j] = B[i,j]-B[i+1,j]
                end
            end
        end
    end
end
"""
    dvand!(a, b) -> b

Solves system ``A*x = b`` in-place.

``A`` is Vandermonde matrix ``A_{ij} = a_i^{j-1}``.

Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function dvand!(alpha, B)
    n = length(alpha)
    if n != size(B,1)
        throw(DimensionMismatch("matrix has dimensions ($n,$n) but right hand side has $(size(B,1)) rows"))
    end
    nrhs = size(B,2)
    @inbounds begin
        for j=1:nrhs
            for k=1:n-1
                for i=n:-1:k+1
                    B[i,j] = (B[i,j]-B[i-1,j])/(alpha[i]-alpha[i-k])
                end
            end
            for k=n-1:-1:1
                for i=k:n-1
                    B[i,j] = B[i,j]-alpha[k]*B[i+1,j]
                end
            end
        end
    end
end
