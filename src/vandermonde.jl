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

if VERSION >= v"0.7.0" # Only use Adjoint and Transpose for 0.7 and up

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

function \(V::Adjoint{T1,Vandermonde{T1}}, y::Vector{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, length(y))
    pvand!(x, adjoint(V.parent.c), y)
end    

function \(V::Transpose{T1,Vandermonde{T1}}, y::Vector{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, length(y))
    pvand!(x, V.parent.c, y)
end

end # End version check

function \(V::Vandermonde{T1}, y::Vector{T2}) where T1 where T2
    T = vandtype(T1,T2)
    x = Array{T}(undef, length(y))    
    dvand!(x, V.c, y)
end


function vandtype(T1::Type, T2::Type)
    # Figure out the return type of Vandermonde{T1} \ Vector{T2}
    T = promote_type(T1, T2)
    S = typeof(oneunit(T)/oneunit(T1))
    return S
end

"""
    pvand!(x, alpha, b) -> x

Solves transposed system ``A^T*x = b`` \\
A is Vandermonde matrix, \\
``A_{ij} = α_i^{j-1}``

Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function pvand!(x, alpha, b)
    n = length(alpha);
    copyto!(x,b)
    for k=1:n
        for j=n:-1:k+1
            x[j] = x[j]-alpha[k]*x[j-1]
        end
    end
    for k=n-1:-1:1
        for j=k+1:n
            x[j] = x[j]/(alpha[j]-alpha[j-k])
        end
        for j=k:n-1
            x[j] = x[j]-x[j+1]
        end
    end
    return x
end

"""
    dvand(x, alpha, b) -> x

Solves system ``A*x = b`` \\
A is Vandermonde matrix, \\
``A_{ij} = α_i^{j-1}``

Algorithm by Bjorck & Pereyra,
Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903,
https://doi.org/10.2307/2004623
"""
function dvand!(x, alpha, b)
    n = length(alpha)
    copyto!(x,b)
    for k=1:n-1
        for j=n:-1:k+1
            x[j] = (x[j]-x[j-1])/(alpha[j]-alpha[j-k])
        end
    end
    for k=n-1:-1:1
        for j=k:n-1
            x[j] = x[j]-alpha[k]*x[j+1]
        end
    end
    return x
end

