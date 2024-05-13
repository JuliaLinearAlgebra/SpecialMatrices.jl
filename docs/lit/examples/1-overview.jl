#=
# [SpecialMatrices overview](@id 1-overview)

This page illustrates the Julia package
[`SpecialMatrices`](https://github.com/JuliaLinearAlgebra/SpecialMatrices.jl).

This package extends the `LinearAlgebra` library
with support for special matrices that are used in linear algebra.
Every special matrix has its own type
and is stored efficiently.
Use `Matrix(A)` to access the full matrix if needed.

## Related packages

[ToeplitzMatrices.jl](https://github.com/JuliaLinearAlgebra/ToeplitzMatrices.jl)
supports
Toeplitz, Hankel, and circulant matrices.
=#

#srcURL


# ### Setup

# Packages needed here.

using SpecialMatrices
using Polynomials
using LinearAlgebra: factorize, Diagonal


# ## [`Cauchy` matrix](https://en.wikipedia.org/wiki/Cauchy_matrix)

Cauchy(1:3, 2:4)

#
Cauchy((1., 2.), 1:3)

#
Cauchy(3)


# ## [`Companion` matrix](https://en.wikipedia.org/wiki/Companion_matrix)

Companion(1:3)

# From a polynomial
p = Polynomial(4:-1:1)

#
Companion(p)


# ## [`Frobenius` matrix](https://en.wikipedia.org/wiki/Frobenius_matrix)

F = Frobenius(3, 2:4) # Specify subdiagonals of column 3

# Special form of inverse:
inv(F)

# Special form preserved if the same column has the subdiagonals
F * F

# Otherwise it promotes to a `Matrix`:
F * Frobenius(2, 2:5)

# Efficient matrix-vector multiplication:
F * (1:6)


# ## [`Hilbert` matrix](https://en.wikipedia.org/wiki/Hilbert_matrix)

H = Hilbert(5)

# Inverses are also integer matrices:
inv(H)


#=
## [`Kahan` matrix](https://math.nist.gov/MatrixMarket/data/MMDELI/kahan/kahan.html)

See [N. J. Higham (1987)](https://doi.org/10.1137/1029112).
=#


Kahan(5, 5, 1, 35)

#
Kahan(5, 3, 0.5, 0)

#
Kahan(3, 5, 0.5, 1e-3)


#=
## `Riemann` matrix

A Riemann matrix is defined as
`A = B[2:N+1, 2:N+1]`,
where
`B[i,j] = i-1` if `i` divides `j`, and `-1` otherwise.
The [Riemann hypothesis](https://en.wikipedia.org/wiki/Riemann_hypothesis)
holds if and only if
`det(A) = O( N! N^(-1/2+ϵ))` for every `ϵ > 0`.

See
[F. Roesler (1986)](https://doi.org/10.1016/0024-3795(86)90255-7).
=#

Riemann(5)


#=
## [`Strang` matrix](https://doi.org/10.1137/141000671)

A special symmetric, tridiagonal, Toeplitz matrix named after Gilbert Strang.
=#

S = Strang(5)

# The Strang matrix has a special ``L D L'`` factorization:
F = factorize(S)

# Here is a verification:
F.L * Diagonal(F.D) * F.L'


# ## [`Vandermonde` matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix)

a = 1:4
V = Vandermonde(a)

# Adjoint Vandermonde:

V'

#=
The backslash operator `\` is overloaded
to solve Vandermonde and adjoint Vandermonde systems
in ``O(n^2)`` time using the algorithm of
[Björck & Pereyra (1970)](https://doi.org/10.2307/2004623).
=#
V \ a

#
V' \ V[3,:]


include("../../../inc/reproduce.jl")
