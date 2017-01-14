export Circulant

immutable Circulant{T} <: AbstractMatrix{T}
    c::Vector{T}
end
typealias Circ Circulant

function full{T}(C::Circ{T})
    n = size(C, 1)
    M = Array(T, n, n)
    for i=1:n
        M[i:n,i] = C.c[1:n-i+1]
        M[1:i-1,i] = C.c[n-i+2:n]
    end
    M
end

getindex(C::Circ, i::Int, j::Int) = C.c[mod(i-j,length(C.c))+1]
isassigned(C::Circ, i::Int, j::Int) = isassigned(C.c,mod(i-j,length(C.c))+1)
size(C::Circ, r::Int) = (r==1 || r==2) ? length(C.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(C::Circ) = size(C,1), size(C,2)

Base.copy(C::Circ) = Circ(C.c)
Base.conj(C::Circ) = Circ(conj(C.c))
Base.conj!(C::Circ) = (conj!(C.c); C)

Base.transpose(C::Circ) = Circ(circshift(reverse(C.c), 1))
Base.ctranspose(C::Circ) = conj!(transpose(C))

+(C::Circ, a::Number) = Circ(C.c + a)
+(a::Number, C::Circ) = Circ(a + C.c)
+(C::Circ, D::Circ) = Circ(C.c + D.c)
-(C::Circ, a::Number) = Circ(C.c - a)
-(a::Number, C::Circ) = Circ(a - C.c)
-(C::Circ, D::Circ) = Circ(C.c - D.c)
-(C::Circ) = Circ(-C.c)

# # Return the eigen-factorisation of C. Constituents are efficient to compute
# # and the eigenvectors are represented efficiently using a UnitaryDFT object.
# function Base.eigfact{T}(C::Circ{T})
#     n = length(C.c)
#     Base.LinAlg.Eigen(DFT{T}(n) * C.c, UnitaryDFT{T}(n))
# end

eig{T}(C::Circ{T}) = (N = length(C.c); (DFT{T}(N) * C.c, UnitaryDFT{T}(N)))

# Generic code for multiplication by a Circulant matrix. Dispatch is used later
# on to ensure that the types returned are the correct ones. The idea is that
# if both inputs are Integer, then the output should also be integer.
# Similarly, if the inputs are both real, then the output should also be real.
function circ_dot{T}(C::Circ{T}, X::StridedVecOrMat)
    if size(X, 1) != size(C, 2)
        throw(ArgumentError("C and X are not conformal."))
    end
    γ, U = eig(C)
    return Ac_mul_B(U, (Diagonal(γ) * (U * X)))
end
function circ_dot(X::StridedVecOrMat, C::Circ)
    if size(X, 2) != size(C, 1)
        println(size(X, 2))
        println(size(C, 1))
        throw(ArgumentError("C and X are not conformal."))
    end
    γ, U = eig(C)
    return (A_mul_Bc(X, U) * Diagonal(γ)) * U
end

# In general, the return type of the * operation should be a complex float.
*(C::Circ, X::SVM) = circ_dot(C, X)
*(X::SVM, C::Circ) = circ_dot(X, C)

# If both arguments contain reals, then the result should really contain reals.
*{T<:Real, V<:Real}(C::Circ{T}, X::SVM{V}) = real!.(circ_dot(C, X))
*{T<:Real, V<:Real}(X::SVM{V}, C::Circ{T}) = real!.(circ_dot(X, C))

# TODO: How attached are people to this? This should probably change. Also, this doesn't
# really get you anything in terms of performance. There's not reduction in memory use as far as I can see
# as the in-place fourier transform is not used.
function A_mul_B!{T}(y::StridedVector{T},C::Circulant{T},x::StridedVector{T})
    xt=fft(x)
    vt=fft(C.c)
    yt=ifft(vt.*xt)
    if T<: Int
        map!(round,y,yt) 
    elseif T<: Real
        map!(real,y,yt)
    else
        copy!(y,yt)
    end
    return y
end

# TODO: Implement all of the other functionality that you would expect to find with a matrix.
