export Vandermonde

"""
    V = Vandermonde(c::AbstractVector)

Create a "lazy" `n × n`
[`Vandermonde` matrix](http://en.wikipedia.org/wiki/Vandermonde_matrix)
but requiring only `O(n)` storage for the vector `c`.

```jldoctest van
julia> a = 1:5; A = Vandermonde(a)
5×5 Vandermonde{Int64}:
 1  1   1    1    1
 1  2   4    8   16
 1  3   9   27   81
 1  4  16   64  256
 1  5  25  125  625
```

Adjoint Vandermonde:
```jldoctest van
julia> A'
5×5 adjoint(::Vandermonde{Int64}) with eltype Int64:
 1   1   1    1    1
 1   2   3    4    5
 1   4   9   16   25
 1   8  27   64  125
 1  16  81  256  625
```

The backslash operator `\\` is overloaded to solve Vandermonde and adjoint
Vandermonde systems in `O(n^2)` time using the algorithm of
Björck & Pereyra (1970), https://doi.org/10.2307/2004623.

```jldoctest van
julia> A \\ a
5-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0
```

```jldoctest van
julia> A' \\ A[2,:]
5-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0
```
"""
struct Vandermonde{T} <: AbstractMatrix{T}
    c :: AbstractVector{T}

    function Vandermonde(c::AbstractVector{T}) where T
        axes(c,1) isa Base.OneTo || throw(ArgumentError("must be OneTo"))
        new{T}(c)
    end
end

@inline Base.@propagate_inbounds function getindex(
    V::Vandermonde,
    i::Int,
    j::Int,
)
    @boundscheck checkbounds(V, i, j)
    return (@inbounds V.c[i])^(j-1)
end

# not needed: comes from size(V)
#isassigned(V::Vandermonde, i::Int, j::Int) = isassigned(V.c, i)

# not needed: comes from size(V)
#size(V::Vandermonde, r::Int) = (r==1 || r==2) ? length(V.c) :
#   throw(ArgumentError("Invalid dimension $r"))
size(V::Vandermonde) = (length(V.c), length(V.c))

#=
# not needed: comes from getindex(V)
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
=#

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
    S = typeof(oneunit(T) / oneunit(T1))
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
    n == size(B,1) ||
        throw(DimensionMismatch("matrix has dimensions ($n,$n) but right hand side has $(size(B,1)) rows"))
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
