export Circulant, CircEig

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

# Return the eigen-factorisation of C. Constituents are efficient to compute
# and the eigenvectors are represented efficiently using a UnitaryDFT object.
eig{T<:Real}(C::Circ{T}) = (N = length(C.c); (DFT{Complex{T}}(N) * C.c, UDFT{Complex{T}}(N)))
eig{T<:Complex}(C::Circ{T}) = (N = length(C.c); (DFT{T}(N) * C.c, UDFT{T}(N)))
eigfact(C::Circ) = Eigen(eig(C)...)
typealias CircEig{T<:Complex} Eigen{T, T, UDFT{T}, Vector{T}}

# Check for conformal arguments.
function conf{T}(C::CircEig{T}, X::SVM)
    if size(X, 1) != size(C.vectors, 2)
        throw(ArgumentError("C and X are not conformal."))
    end
end
function conf{T}(X::SVM, C::CircEig{T})
    if size(X, 2) != size(C.vectors, 1)
        throw(ArgumentError("X and C are not conformal."))
    end
end

# TODO: How do I correctly do the in-place bits of this? Do some profiling
# before moving forwards with more functionality.

# Perform the multiplication A * B in place, overwriting B.
# TODO: TEST THIS!
function A_mul_B!(A::Diagonal, B::SVM)
    M, N = size(B)
    if size(A, 2) != M
        throw(ArgumentError("A and B are not conformal."))
    end
    d = D.diag
    @show M, N
    for j in 1:N
        for i in 1:M
            B[i, j] *= d[i]
        end
    end
end

# Always compute the eigen-factorisation to compute the multiplication. If
# this operation is to be performed multiple times, then we should just cache
# the factorisation to save re-computing it for each multiplication operation.
function *{T}(C::CircEig{T}, X::SVM)
    conf(C, X)
    Γ, U = Diagonal(C.values), C.vectors
    Ac_mul_B(U, (Γ * (U * X)))
    tmp = U * X
    bfft!(A_mul_B!(tmp, Γ, tmp), 1) / U.sqrtN
    @show Ac_mul_B(U, (Γ * (U * X)))
    tmp = U * X
    @show bfft!(A_mul_B!(Γ, tmp), 1) / U.sqrtN
end
function *{T}(X::SVM, C::CircEig{T})
    conf(X, C)
    Γ, U = Diagonal(C.values), C.vectors
    (A_mul_Bc(X, U) * Γ) * U
end
*(C::Circ, X::SVM) = eigfact(C) * X
*(X::SVM, C::Circ) = X * eigfact(C)

# Use the eigen-factorisation to compute the inv(C) * X. The inverse
# factorisation should be cached if the operation is required multiple times.
function \{T}(C::CircEig{T}, X::SVM)
    conf(C, X)
    Γ, U = Diagonal(C.values), C.vectors
    Ac_mul_B(U, (Γ \ (U * X)))
end
function /{T}(X::SVM, C::CircEig{T})
    conf(X, C)
    Γ, U = Diagonal(C.values), C.vectors
    (A_mul_Bc(X, U) / Γ) * U
end
\(C::Circ, X::SVM) = eigfact(C) \ X
/(X::SVM, C::Circ) = X / eigfact(C)

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
