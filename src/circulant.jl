export Circulant, CircEig, tocirc

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

# getindex(C::Circ, i::Int, j::Int) = C.c[mod(i-j,length(C.c))+1]
isassigned(C::Circ, i::Int, j::Int) = isassigned(C.c,mod(i-j,length(C.c))+1)
size(C::Circ, r::Int) = (r==1 || r==2) ? length(C.c) :
    throw(ArgumentError("Invalid dimension $r"))
size(C::Circ) = size(C,1), size(C,2)

copy(C::Circ) = Circ(C.c)
conj(C::Circ) = Circ(conj(C.c))
conj!(C::Circ) = (conj!(C.c); C)

transpose(C::Circ) = Circ(circshift(reverse(C.c), 1))
ctranspose(C::Circ) = conj!(transpose(C))

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

# Utility function to convert back to usual circulant formulation.
tocirc{T}(C::CircEig{T}) = (N = length(C.values); Circ(DFT{T}(N)'C.values / N))

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

# Always compute the eigen-factorisation to compute the multiplication. If
# this operation is to be performed multiple times, then we should just cache
# the factorisation to save re-computing it for each multiplication operation.
function *{T}(C::CircEig{T}, X::SVM)
    conf(C, X)
    Γ, U = Diagonal(C.values), C.vectors
    Ac_mul_B(U, (Γ * (U * X)))
    Y = U * X
    Ac_mul_B!(U, Γ * (U * X))
end
function *{T}(X::SVM, C::CircEig{T})
    conf(X, C)
    Γ, U = Diagonal(C.values), C.vectors
    (A_mul_Bc(X, U) * Γ) * U
    Y = A_mul_Bc(X, U)
    A_mul_B!(A_mul_B!(Y, Y, Γ), U)
end
*(C::Circ, X::SVM) = eigfact(C) * X
*(X::SVM, C::Circ) = X * eigfact(C)

# "In-place" multiplication.
function A_mul_B!{T<:Number, V<:Complex}(Y::SVM{V}, C::CircEig{T}, X::SVM{V})
    conf(C, X) 
    Γ, U = Diagonal(C.values), C.vectors
    Ac_mul_B!(U, A_mul_B!(Y, Γ, A_mul_B!(U, X)))
end
function A_mul_B!{T<:Number, V<:Complex}(Y::SVM{V}, X::SVM{V}, C::CircEig{T})
    conf(X, C)
    Γ, U = Diagonal(C.values), C.vectors
    A_mul_B!(A_mul_B!(Y, A_mul_Bc!(X, U), Γ), U)
end
A_mul_B!(Y::SVM, C::Circ, X::SVM) = A_mul_B!(Y, eigfact(C), X)
A_mul_B!(Y::SVM, X::SVM, C::Circ)  = A_mul_B!(Y, X, eigfact(C))

inv{T}(C::CircEig{T}) = Eigen(1 ./ C.values, UDFT{T}(length(C.values)))
inv(C::Circ) = tocirc(inv(eigfact(C)))

det{T}(C::CircEig{T}) = prod(C.values)
det(C::Circ) = det(eigfact(C))

logdet{T}(C::CircEig{T}) = sum(log(C.values))
logdet(C::Circ) = logdet(eigfact(C))

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

# Helper functionality for casting.
real(C::Circ) = Circ(real(C.c))
convert{T}(::Type{Circ{T}}, C::Circ) = Circ(convert(Vector{T}, C.c))

round(C::Circ) = Circ(round(C.c))
round{T}(::Type{Circ{T}}, C::Circ) = Circ(round(Vector{T}, C.c))
