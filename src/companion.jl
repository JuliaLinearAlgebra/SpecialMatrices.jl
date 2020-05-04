export Companion

struct Companion{T} <: AbstractArray{T, 2}
    c :: Vector{T}
end

# Generate companion matrix from a polynomial

function Companion(P::Polynomial{T}) where T
   n = length(P)
   c = Array{T}(undef, n-1)
   d=P.a[n]
   for i=1:n-1
       c[i]=P.a[i]/d
   end
   Companion(c)   
end

#Basic property computations
size(C::Companion, r::Int) = (r==1 || r==2) ? length(C.c) : 
    throw(ArgumentError("Companion is of rank 2"))

function size(C::Companion)
    n = length(C.c)
    n, n
end

#XXX Inefficient but works
# getindex(C::Companion, i, j) = getindex(Matrix(C), i, j)
# isassigned(C::Companion, i, j) = isassigned(Matrix(C), i, j)
getindex(C::Companion{T}, i::Int, j::Int) where T = (j==length(C.c)) ? -C.c[i] : (i==j+1 ? one(T) : zero(T) )
isassigned(C::Companion, i::Int, j::Int) = (j==length(C.c)) ? isassigned(C.c,i) : true


function Matrix(C::Companion{T}) where T
    M = zeros(T, size(C)...)
    M[:,end]=-C.c
    for i=1:size(C,1)-1
    	M[i+1, i] = one(T)
	end
    M
end

#Linear algebra stuff
function mul!(C::Companion{T}, b::Vector{T}) where T
	x = b[end]
	y = -C.c[1]*x
	b[2:end] = b[1:end-1]-C.c[2:end]*x
	b[1] = y
    b
end
*(C::Companion{T}, b::Vector{T}) where T = mul!(C, copy(b))

function mul!(A::Matrix{T}, C::Companion{T}) where T
	v = Array{T}(undef, size(A,1))
	for i=1:size(A,1)
		v[i] =dot(vec(A[i,:]),-C.c)
	end
	for i=1:size(A,1), j=1:size(A,2)-1
		A[i,j] = A[i,j+1]
	end
	A[:,end] = v
	A
end
*(A::Matrix{T}, C::Companion{T}) where T = mul!(copy(A), C)

function inv(C::Companion{T}) where T
	M = zeros(T, size(C)...)
    for i=1:size(C,1)-1
    	M[i, i+1] = one(T)
	end
	d = M[end, 1] = -one(T)/C.c[1]
	M[1:end-1, 1] = d*C.c[2:end]
    M
end
