# SpecialMatrices.jl

[![action status][action-img]][action-url]
[![pkgeval status][pkgeval-img]][pkgeval-url]
[![codecov][codecov-img]][codecov-url]
[![coveralls][coveralls-img]][coveralls-url]
[![license][license-img]][license-url]
[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
[![deps](https://juliahub.com/docs/SpecialMatrices/deps.svg)](https://juliahub.com/ui/Packages/SpecialMatrices)
[![version](https://juliahub.com/docs/SpecialMatrices/version.svg)](https://juliahub.com/ui/Packages/SpecialMatrices)
[![pkgeval](https://juliahub.com/docs/SpecialMatrices/pkgeval.svg)](https://juliahub.com/ui/Packages/SpecialMatrices)

https://github.com/JuliaLinearAlgebra/SpecialMatrices.jl

A [Julia](https://julialang.org) package for working with special matrix types.

This Julia package extends the `LinearAlgebra` library
with support for special matrices that are used in linear algebra.
Every special matrix has its own type
and is stored efficiently.
The full matrix is accessed by the command `Matrix(A)`.


## Installation

```julia
julia> ] add SpecialMatrices
```


## Related packages

[ToeplitzMatrices.jl](https://github.com/JuliaLinearAlgebra/ToeplitzMatrices.jl)
supports
Toeplitz, Hankel, and circulant matrices.


## Currently supported special matrices

### [`Cauchy` matrix](https://en.wikipedia.org/wiki/Cauchy_matrix)

```julia
Cauchy(x,y)[i,j] = 1/(x[i] + y[j])
Cauchy(x) = Cauchy(x,x)
Cauchy(k::Int) = Cauchy(1:k)

julia> Cauchy([1,2,3],[3,4,5])
3×3 Cauchy{Int64}:
 0.25      0.2       0.166667
 0.2       0.166667  0.142857
 0.166667  0.142857  0.125

julia> Cauchy([1,2,3])
3×3 Cauchy{Int64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667

julia> Cauchy(3)
3×3 Cauchy{Float64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667
```


### [`Companion` matrix](https://en.wikipedia.org/wiki/Companion_matrix)

```julia
julia> A=Companion([3,2,1])
3×3 Companion{Int64}:
 0  0  -3
 1  0  -2
 0  1  -1
```
Also, directly from a polynomial:

```julia
julia> using Polynomials

julia> P=Polynomial([2.0,3,4,5])
Polynomial(2 + 3x + 4x^2 + 5x^3)

julia> C=Companion(P)
3×3 Companion{Float64}:
 0.0  0.0  -0.4
 1.0  0.0  -0.6
 0.0  1.0  -0.8
```


### [`Frobenius` matrix](https://en.wikipedia.org/wiki/Frobenius_matrix)

Example

```julia
julia> using SpecialMatrices

julia> F=Frobenius(3, [1.0,2.0,3.0]) #Specify subdiagonals of column 3
6×6 Frobenius{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  0.0  0.0
 0.0  0.0  2.0  0.0  1.0  0.0
 0.0  0.0  3.0  0.0  0.0  1.0

julia> inv(F) #Special form of inverse
6×6 Frobenius{Float64}:
 1.0  0.0   0.0  0.0  0.0  0.0
 0.0  1.0   0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0
 0.0  0.0  -1.0  1.0  0.0  0.0
 0.0  0.0  -2.0  0.0  1.0  0.0
 0.0  0.0  -3.0  0.0  0.0  1.0

julia> F*F #Special form preserved if the same column has the subdiagonals
6×6 Frobenius{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  2.0  1.0  0.0  0.0
 0.0  0.0  4.0  0.0  1.0  0.0
 0.0  0.0  6.0  0.0  0.0  1.0

julia> F*Frobenius(2, [5.0,4.0,3.0,2.0]) #Promotes to Matrix
6×6 Matrix{Float64}:
 1.0   0.0  0.0  0.0  0.0  0.0
 0.0   1.0  0.0  0.0  0.0  0.0
 0.0   5.0  1.0  0.0  0.0  0.0
 0.0   9.0  1.0  1.0  0.0  0.0
 0.0  13.0  2.0  0.0  1.0  0.0
 0.0  17.0  3.0  0.0  0.0  1.0

julia> F*[10.0,20,30,40,50,60.0]
6-element Vector{Float64}:
  10.0
  20.0
  30.0
  70.0
 110.0
 150.0
```


### [`Hilbert` matrix](https://en.wikipedia.org/wiki/Hilbert_matrix)

```julia
julia> A=Hilbert(5)
5×5 Hilbert{Rational{Int64}}:
 1//1  1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9
```
Inverses are also integer matrices:

```julia
julia> inv(A)
5×5 InverseHilbert{Rational{Int64}}:
    25//1    -300//1     1050//1    -1400//1     630//1
  -300//1    4800//1   -18900//1    26880//1  -12600//1
  1050//1  -18900//1    79380//1  -117600//1   56700//1
 -1400//1   26880//1  -117600//1   179200//1  -88200//1
   630//1  -12600//1    56700//1   -88200//1   44100//1
```


### [`Kahan` matrix](https://math.nist.gov/MatrixMarket/data/MMDELI/kahan/kahan.html)

```julia
julia> Kahan(5,5,1,35)
5×5 Kahan{Int64,Int64}:
 1.0  -0.540302  -0.540302  -0.540302  -0.540302
 0.0   0.841471  -0.454649  -0.454649  -0.454649
 0.0   0.0        0.708073  -0.382574  -0.382574
 0.0   0.0        0.0        0.595823  -0.321925
 0.0   0.0        0.0        0.0        0.501368

julia> Kahan(5,3,0.5,0)
5×3 Kahan{Float64,Int64}:
 1.0  -0.877583  -0.877583
 0.0   0.479426  -0.420735
 0.0   0.0        0.229849
 0.0   0.0        0.0
 0.0   0.0        0.0

julia> Kahan(3,5,0.5,1e-3)
3×5 Kahan{Float64,Float64}:
 1.0  -0.877583  -0.877583  -0.877583  -0.877583
 0.0   0.479426  -0.420735  -0.420735  -0.420735
 0.0   0.0        0.229849  -0.201711  -0.201711
```

For more details see [N. J. Higham (1987)][Higham87].

[Higham87]: https://eprints.ma.man.ac.uk/695/01/covered/MIMS_ep2007_10.pdf "N. Higham, A Survey of Condition Number Estimation for Triangular Matrices, SIMAX, Vol. 29, No. 4, pp. 588, 1987"


### `Riemann` matrix

Riemann matrix is defined as `A = B[2:N+1, 2:N+1]`, where
`B[i,j] = i-1` if `i` divides `j`, and `-1` otherwise.
[Riemann hypothesis](https://en.wikipedia.org/wiki/Riemann_hypothesis) holds
if and only if `det(A) = O( N! N^(-1/2+ϵ))` for every `ϵ > 0`.

```julia
julia> Riemann(7)
7×7 Riemann{Int64}:
  1  -1   1  -1   1  -1   1
 -1   2  -1  -1   2  -1  -1
 -1  -1   3  -1  -1  -1   3
 -1  -1  -1   4  -1  -1  -1
 -1  -1  -1  -1   5  -1  -1
 -1  -1  -1  -1  -1   6  -1
 -1  -1  -1  -1  -1  -1   7
```

For more details see [F. Roesler (1986)][Roesler1986].

[Roesler1986]: https://doi.org/10.1016/0024-3795(86)90255-7 "Friedrich Roesler,
Riemann's hypothesis as an eigenvalue problem,
Linear Algebra and its Applications, Vol. 81, p.153-198, Sep. 1986"


### [`Strang` matrix](https://doi.org/10.1137/141000671)

A special symmetric, tridiagonal, Toeplitz matrix named after Gilbert Strang.

```julia
julia> Strang(6)
6×6 Strang{Int64}:
  2  -1   0   0   0   0
 -1   2  -1   0   0   0
  0  -1   2  -1   0   0
  0   0  -1   2  -1   0
  0   0   0  -1   2  -1
  0   0   0   0  -1   2
```


### [`Vandermonde` matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix)

```julia
julia> a = 1:5
julia> A = Vandermonde(a)
5×5 Vandermonde{Int64}:
 1  1   1    1    1
 1  2   4    8   16
 1  3   9   27   81
 1  4  16   64  256
 1  5  25  125  625
```

Adjoint Vandermonde:
```julia
julia> A'
5×5 adjoint(::Vandermonde{Int64}) with eltype Int64:
 1   1   1    1    1
 1   2   3    4    5
 1   4   9   16   25
 1   8  27   64  125
 1  16  81  256  625
```

The backslash operator `\` is overloaded
to solve Vandermonde and adjoint Vandermonde systems
in ``O(n^2)`` time using the algorithm of
[Björck & Pereyra (1970)](https://doi.org/10.2307/2004623).
```julia
julia> A \ a
5-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0

julia> A' \ A[2,:]
5-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0
```

<!-- URLs -->
[action-img]: https://github.com/JuliaLinearAlgebra/SpecialMatrices.jl/workflows/CI/badge.svg
[action-url]: https://github.com/JuliaLinearAlgebra/SpecialMatrices.jl/actions
[build-img]: https://github.com/JuliaLinearAlgebra/SpecialMatrices.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaLinearAlgebra/SpecialMatrices.jl/actions?query=workflow%3ACI+branch%3Amaster
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SpecialMatrices.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/SpecialMatrices.html
[code-blue-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[code-blue-url]: https://github.com/invenia/BlueStyle
[codecov-img]: https://codecov.io/github/JuliaLinearAlgebra/SpecialMatrices.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/JuliaLinearAlgebra/SpecialMatrices.jl?branch=master
[coveralls-img]: https://coveralls.io/repos/JuliaLinearAlgebra/SpecialMatrices.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaLinearAlgebra/SpecialMatrices.jl?branch=master
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaLinearAlgebra.github.io/SpecialMatrices.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaLinearAlgebra.github.io/SpecialMatrices.jl/dev
[license-img]: https://img.shields.io/badge/license-MIT-brightgreen.svg
[license-url]: LICENSE
