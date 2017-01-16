export DFT, UnitaryDFT, UDFT

# Non-unitary DFT matrix.
abstract AbstractDFT{T} <: AbstractMatrix{T}
immutable DFT{T<:Complex} <: AbstractDFT{T}
    N::Int
end

# Unitary DFT Matrix.
immutable UnitaryDFT{T<:Complex} <: AbstractDFT{T}
    N::Int
    sqrtN::Float64
    UnitaryDFT(N) = new(N, sqrt(N))
end
typealias UDFT UnitaryDFT

Base.show(io::IO, U::DFT) = print(io, "DFT matrix of size $U.N.")
Base.show(io::IO, U::UDFT) = print(io, "Unitary DFT matrix of size $U.N.")

size(M::AbstractDFT) = M.N, M.N
length(M::AbstractDFT) = M.N^2

# getindex(M::DFT, m::Int, n::Int) =
#     exp(-2π * im * (m-1) * (n-1) / M.N)
# getindex(M::UnitaryDFT, m::Int, n::Int) =1
#     exp(-2π * im * (m-1) * (n-1) / M.N) / M.sqrtN
# getindex(M::AbstractDFT, m::Int) = getindex(M, rem(m, M.N), div(m, M.N))

# isassigned(M::AbstractDFT, i::Int) = 1 <= i <= M.N ? true : false

# Helper functions to check whether stuff is conformal or not.
function conf(U::AbstractDFT, x::SV)
    if U.N != length(x)
        throw(ArgumentError("U and x are not conformal (in this order)."))
    end
end
# Always throw an error as it does not make sense to multiply a column Vector
# by a matrix in this order.
function conf(x::SV, U::AbstractDFT)
    throw(ArgumentError("x and U are not conformal (in this order)."))
end
function conf(U::AbstractDFT, X::SM)
    if U.N != size(X, 1)
        throw(ArgumentError("U and X are not conformal (in this order)."))
    end
end
function conf(X::SM, U::AbstractDFT)
    if size(X, 2) != U.N
        throw(ArgumentError("X and U are not conformal (in this order)."))
    end
end

# Multiplication from the left and right by the DFT matrices.
*(U::DFT, X::SVM) = (conf(U, X); fft(X, 1))
*(X::SVM, U::DFT) = (conf(X, U); fft(X, 2))
*(U::UDFT, X::SVM) = (conf(U, X); fft(X, 1) / U.sqrtN)
*(X::SVM, U::UDFT) = (conf(X, U); fft(X, 2) / U.sqrtN)

# In-place functionality. Note that diagonal dispatch is used as both arguments should be of the same type.
A_mul_B!(Y::SV, U::DFT, X::SV) = (conf(U, X); copy!(Y, X); fft!(Y))
A_mul_B!(Y::SM, U::DFT, X::SM) = (conf(U, X); copy!(Y, X); fft!(Y, 1))
A_mul_B!(Y::SV, X::SV, U::DFT) = (conf(X, U); copy!(Y, X); fft!(Y))
A_mul_B!(Y::SM, X::SM, U::DFT) = (conf(X, U); copy!(Y, X); fft!(Y, 2))
A_mul_B!(Y::SV, U::UDFT, X::SV) = (conf(U, X); copy!(Y, X); fft!(Y); y ./= U.sqrtN)
A_mul_B!(Y::SM, U::UDFT, X::SM) = (conf(U, X); copy!(Y, X); fft!(Y, 1); y ./= U.sqrtN)
A_mul_B!(Y::SV, X::SV, U::UDFT) = (conf(X, U); copy!(Y, X); fft!(Y, 2); y ./= U.sqrtN)
A_mul_B!(Y::SM, X::SM, U::UDFT) = (conf(X, U); copy!(Y, X); fft!(Y, 2); y ./= U.sqrtN)

# Conjugate transpose of U'X and X*U'.
Ac_mul_B(U::DFT, X::SVM) = (conf(U, X); bfft(X, 1))
A_mul_Bc(X::SVM, U::DFT) = (conf(X, U); bfft(X, 2))
Ac_mul_B(U::UDFT, X::SVM) = (conf(U, X); bfft(X, 1) / U.sqrtN)
A_mul_Bc(X::SVM, U::UDFT) = (conf(X, U); bfft(X, 2) / U.sqrtN)

# Explicit implementation of X'U and U*X' to avoid fallback in Base being called.
Ac_mul_B(X::SVM, U::AbstractDFT) = (Xc = ctranspose(X); conf(Xc, U); A_mul_B!(Xc, Xc, U))
A_mul_Bc(U::AbstractDFT, X::SVM) = (Xc = ctranspose(X); conf(U, Xc); A_mul_B!(Xc, U, Xc))

# All DFT matrices are symmetric.
At_mul_B(U::AbstractDFT, X::SVM) = U * X
A_mul_Bt(X::SVM, U::AbstractDFT) = X * U
