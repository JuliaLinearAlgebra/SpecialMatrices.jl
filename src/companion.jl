export CompanionMatrix

immutable CompanionMatrix{T} <: AbstractArray{T, 2}
    c :: Vector{T}
end

#Basic property computations
size(C::CompanionMatrix, r::Int) = (r==1 || r==2) ? length(C.c) : 
    throw(ArgumentError("CompanionMatrix is of rank 2"))

function size(C::CompanionMatrix)
    n = length(C.c)
    n, n
end

#XXX Inefficient but works
getindex(C::CompanionMatrix, i, j) = getindex(full(C), i, j)
isassigned(C::CompanionMatrix, i, j) = isassigned(full(C), i, j)

function full{T}(C::CompanionMatrix{T})
    M = zeros(T, size(C)...)
    M[:,end]=-C.c
    for i=1:size(C,1)-1
    	M[i+1, i] = one(T)
	end
    M
end

#Linear algebra stuff
function A_mul_B!{T}(C::CompanionMatrix{T}, b::Vector{T})
	x = b[end]
	y = -C.c[1]*x
	b[2:end] = b[1:end-1]-C.c[2:end]*x
	b[1] = y
    b
end
*{T}(C::CompanionMatrix{T}, b::Vector{T}) = A_mul_B!(C, copy(b))

function A_mul_B!{T}(A::Matrix{T}, C::CompanionMatrix{T})
	v = Array(T, size(A,1))
	for i=1:size(A,1)
		v[i] =(A[i,:]*-C.c)[1]
	end
	for i=1:size(A,1), j=1:size(A,2)-1
		A[i,j] = A[i,j+1]
	end
	A[:,end] = v
	A
end
*{T}(A::Matrix{T}, C::CompanionMatrix{T}) = A_mul_B!(copy(A), C)

function inv{T}(C::CompanionMatrix{T})
	M = zeros(T, size(C)...)
    for i=1:size(C,1)-1
    	M[i, i+1] = one(T)
	end
	d = M[end, 1] = -one(T)/C.c[1]
	M[1:end-1, 1] = d*C.c[2:end]
    M
end
