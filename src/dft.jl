export DFT, IDFT, UnitaryDFT, UnitaryIDFT

# Implementation of the DFT matrix. Not really an AbstractArray as it doesn't
# make sense to write to this matrix. Should probably resolve this.

# Non-unitary DFT matrices.
abstract AbstractDFT
immutable DFT <: AbstractDFT
    N::Int
end
immutable IDFT <: AbstractDFT
    N::Int
end

# Unitary DFT matrices.
abstract AbstractUnitaryDFT <: AbstractDFT
immutable UnitaryDFT <: AbstractUnitaryDFT
    N::Int
    sqrtN::Float64
    UnitaryDFT(N) = new(N, sqrt(N))
end
immutable UnitaryIDFT <: AbstractUnitaryDFT
    N::Int
    sqrtN::Float64
    UnitaryIDFT(N) = new(N, sqrt(N))
end

size(M::AbstractDFT) = (M.N, M.N)
length(M::AbstractDFT) = M.N^2

getindex(M::DFT, m::Int, n::Int) =
    exp(-2π * im * (m-1) * (n-1) / M.N)
getindex(N::IDFT, m::Int, n::Int) =
    exp(2π * im * (m-1) * (n-1) / M.N) / M.N
getindex(M::UnitaryDFT, m::Int, n::Int) =
    exp(-2π * im * (m-1) * (n-1) / M.N) / sqrt(M.N)
getindex(M::UnitaryIDFT, m::Int, n::Int) =
    exp(2π * im * (m-1) * (n-1) / M.N) / sqrt(M.N)
getindex(M::AbstractDFT, m::Int) = getindex(M, rem(m, M.N), div(m, M.N))

isassigned(M::AbstractDFT, i::Int) = 1 <= i <= M.N ? true : false

# Helper functions to check whether stuff is conformal or not.
function conformal(U::AbstractDFT, x::Vector)
    if U.N != length(x)
        throw(ArgumentError("U and x are not conformal."))
    end
end
function conformal(U::AbstractDFT, X::Matrix)
    if U.N != size(X, 1)
        throw(ArgumentError("U and X are not conformal."))
    end
end
function conformal(X::Matrix, U::AbstractDFT)
    if size(X, 2) != U.N
        throw(ArgumentError("X and U are not conformal."))
    end
end

# Multiplication from the left by the DFT matrices.
*(U::DFT, x::Vector) = (conformal(U, x); fft(x))
*(U::IDFT, x::Vector) = (conformal(U, x); ifft(x))
*(U::UnitaryDFT, x::Vector) = (conformal(U, x); fft(x) / U.sqrtN)
*(U::UnitaryIDFT, x::Vector) = (conformal(U, x); bfft(x) / U.sqrtN)

*(U::DFT, X::Matrix) = (conformal(U, X); fft(X, 1))
*(U::IDFT, X::Matrix) = (conformal(U, X); ifft(X, 1))
*(U::UnitaryDFT, X::Matrix) = (conformal(U, X); fft(X, 1) / U.sqrtN)
*(U::UnitaryIDFT, X::Matrix) = (conformal(U, X);  bfft(X, 1) / U.sqrtN)

# Multiplication from the right by the DFT matrices.
*(X::Matrix, U::DFT) = (conformal(X, U); fft(X, 2))
*(X::Matrix, U::IDFT) = (conformal(X, U); ifft(X, 2))
*(X::Matrix, U::UnitaryDFT) = (conformal(X, U); fft(X, 2) / U.sqrtN)
*(X::Matrix, U::UnitaryIDFT) = (conformal(X, U); bfft(X, 2) / U.sqrtN)

# In-place functionality.
A_mul_B!(y::Vector, U::DFT, x::Vector) =
    (conformal(U, x); copy!(y, x); fft!(y))
A_mul_B!(y::Vector, U::IDFT, x::Vector) =
    (conformal(U, x); copy!(y, x); ifft!(y))
A_mul_B!(y::Vector, U::UnitaryDFT, x::Vector) =
    (conformal(U, x); copy!(y, x); fft!(y); y ./= U.sqrtN)
A_mul_B!(y::Vector, U::UnitaryIDFT, x::Vector) =
    (conformal(U, x); copy!(y, x); bfft!(y); y ./= U.sqrtN)

A_mul_B!(Y::Matrix, U::DFT, X::Matrix) =
    (conformal(U, X); copy!(Y, X); fft!(Y, 1))
A_mul_B!(Y::Matrix, U::IDFT, X::Matrix) =
    (conformal(U, X); copy!(Y, X); ifft!(Y, 1))
A_mul_B!(Y::Matrix, U::UnitaryDFT, X::Matrix) =
    (conformal(U, X); copy!(Y, X); fft!(Y, 1); y ./= U.sqrtN)
A_mul_B!(Y::Matrix, U::UnitaryIDFT, X::Matrix) =
    (conformal(U, X); copy!(Y, X); bfft!(Y, 1); y ./= U.sqrtN)

A_mul_B!(Y::Matrix, X::Matrix, U::DFT) =
    (conformal(X, U); copy!(Y, X); fft!(Y, 2))
A_mul_B!(Y::Matrix, X::Matrix, U::IDFT) =
    (conformal(X, U); copy!(Y, X); ifft!(Y, 2))
A_mul_B!(Y::Matrix, X::Matrix, U::UnitaryDFT) =
    (conformal(X, U); copy!(Y, X); fft!(Y, 2); y ./= U.sqrtN)
A_mul_B!(Y::Matrix, X::Matrix, U::UnitaryIDFT) =
    (conformal(X, U); copy!(Y, X); bfft!(Y, 2); y ./= U.sqrtN)

# Conjugate transpose of DFT times X.
Ac_mul_B(U::DFT, X::Matrix) = bfft(X, 1)
Ac_mul_B(U::IDFT, X::Matrix) = fft(X, 1) / U.N
Ac_mul_B(U::UnitaryDFT, X::Matrix) = bfft(X, 1) / U.sqrtN
Ac_mul_B(U::UnitaryIDFT, X::Matrix) = fft(X, 1) / U.sqrtN

# All DFT matrices are symmetric.
At_mul_B(U::AbstractDFT, X::Matrix) = U * X

















