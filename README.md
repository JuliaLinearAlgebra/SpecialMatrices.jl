# Special Matrices

A [Julia](http://julialang.org) package for working with special matrix types.

This package extends the `LinearAlgebra` library with support for special
matrices which are used in linear algebra. Every special matrix has its own type.
The full matrix is accessed by the command `Matrix(A)`.

[![Build Status](https://travis-ci.org/JuliaMatrices/SpecialMatrices.jl.svg)](https://travis-ci.org/JuliaMatrices/SpecialMatrices.jl) [![Coverage Status](https://img.shields.io/coveralls/jiahao/SpecialMatrices.jl.svg)](https://coveralls.io/r/jiahao/SpecialMatrices.jl)

## Currently supported special matrices

## [`Cauchy`](http://en.wikipedia.org/wiki/Cauchy_matrix) matrix

```julia
julia> Cauchy([1,2,3],[3,4,5])
3x3 Cauchy{Int64}:
 0.25      0.2       0.166667
 0.2       0.166667  0.142857
 0.166667  0.142857  0.125

julia> Cauchy([1,2,3])
3x3 Cauchy{Int64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667

julia> Cauchy(pi)
3x3 Cauchy{Float64}:
 0.5       0.333333  0.25
 0.333333  0.25      0.2
 0.25      0.2       0.166667
```

## `Circulant` matrix

```julia
julia> Circulant([1,2,3,4])
4x4 Circulant{Int64}:
 1  4  3  2
 2  1  4  3
 3  2  1  4
 4  3  2  1
```

## [`Companion`](http://en.wikipedia.org/wiki/Companion_matrix) matrix

```julia
julia> A=Companion([3,2,1])
3x3 Companion{Int64}:
 0  0  -3
 1  0  -2
 0  1  -1
```
Also, directly from a polynomial:

```julia
julia> using Polynomials

julia> P=Poly([2.0,3,4,5])
Poly(2 + 3x + 4x^2 + 5x^3)

julia> C=Companion(P)
3×3 Companion{Float64}:
 0.0  0.0  -0.4
 1.0  0.0  -0.6
 0.0  1.0  -0.8
```

## [`Frobenius`](http://en.wikipedia.org/wiki/Frobenius_matrix) matrix

Example

```julia
julia> using SpecialMatrices

julia> F=Frobenius(3, [1.0,2.0,3.0]) #Specify subdiagonals of column 3
6x6 Frobenius{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  0.0  0.0
 0.0  0.0  2.0  0.0  1.0  0.0
 0.0  0.0  3.0  0.0  0.0  1.0

julia> inv(F) #Special form of inverse
6x6 Frobenius{Float64}:
 1.0  0.0   0.0  0.0  0.0  0.0
 0.0  1.0   0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0
 0.0  0.0  -1.0  1.0  0.0  0.0
 0.0  0.0  -2.0  0.0  1.0  0.0
 0.0  0.0  -3.0  0.0  0.0  1.0

julia> F*F #Special form preserved if the same column has the subdiagonals
6x6 Frobenius{Float64}:
 1.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  2.0  1.0  0.0  0.0
 0.0  0.0  4.0  0.0  1.0  0.0
 0.0  0.0  6.0  0.0  0.0  1.0

julia> F*Frobenius(2, [5.0,4.0,3.0,2.0]) #Promotes to Matrix
6x6 Array{Float64,2}:
 1.0   0.0  0.0  0.0  0.0  0.0
 0.0   1.0  0.0  0.0  0.0  0.0
 0.0   5.0  1.0  0.0  0.0  0.0
 0.0   9.0  1.0  1.0  0.0  0.0
 0.0  13.0  2.0  0.0  1.0  0.0
 0.0  17.0  3.0  0.0  0.0  1.0

julia> F*[10.0,20,30,40,50,60.0]
6-element Array{Float64,1}:
  10.0
  20.0
  30.0
  70.0
 110.0
 150.0
```

## [`Hankel`](http://en.wikipedia.org/wiki/Hankel_matrix) matrix

Input is a vector of odd length.

```julia
julia> Hankel(collect(-4:4))
5x5 Hankel{Int64}:
 -4  -3  -2  -1  0
 -3  -2  -1   0  1
 -2  -1   0   1  2
 -1   0   1   2  3
  0   1   2   3  4
```

## [`Hilbert`](http://en.wikipedia.org/wiki/Hilbert_matrix) matrix

```julia
julia> A=Hilbert(5)
Hilbert{Rational{Int64}}(5,5)

julia> Matrix(A)
5x5 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9

julia> Matrix(Hilbert(5))
5x5 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9
```
Inverses are also integer matrices:

```julia
julia> inv(A)
5x5 Array{Rational{Int64},2}:
    25//1    -300//1     1050//1    -1400//1     630//1
  -300//1    4800//1   -18900//1    26880//1  -12600//1
  1050//1  -18900//1    79380//1  -117600//1   56700//1
 -1400//1   26880//1  -117600//1   179200//1  -88200//1
   630//1  -12600//1    56700//1   -88200//1   44100//1
```

## [`Kahan`](http://math.nist.gov/MatrixMarket/data/MMDELI/kahan/kahan.html) matrix

```julia
julia> A=Kahan(5,5,1,35)
5x5 Kahan{Int64,Int64}:
 1.0  -0.540302  -0.540302  -0.540302  -0.540302
 0.0   0.841471  -0.454649  -0.454649  -0.454649
 0.0   0.0        0.708073  -0.382574  -0.382574
 0.0   0.0        0.0        0.595823  -0.321925
 0.0   0.0        0.0        0.0        0.501368

julia> A=Kahan(5,3,0.5,0)
5x3 Kahan{Float64,Int64}:
 1.0  -0.877583  -0.877583
 0.0   0.479426  -0.420735
 0.0   0.0        0.229849
 0.0   0.0        0.0
 0.0   0.0        0.0

julia> A=Kahan(3,5,0.5,1e-3)
3x5 Kahan{Float64,Float64}:
 1.0  -0.877583  -0.877583  -0.877583  -0.877583
 0.0   0.479426  -0.420735  -0.420735  -0.420735
 0.0   0.0        0.229849  -0.201711  -0.201711
```

For more details see [N. J. Higham (1987)][Higham87].

[Higham87]: http://eprints.ma.man.ac.uk/695/01/covered/MIMS_ep2007_10.pdf "N. Higham, A Survey of Condition Number Estimation for Triangular Matrices, SIMAX, Vol. 29, No. 4, pp. 588, 1987"

## `Riemann` matrix

Riemann matrix is defined as `A = B[2:N+1, 2:N+1]`, where
`B[i,j] = i-1` if `i` divides `j`, and `-1` otherwise.
[Riemann hypothesis](http://en.wikipedia.org/wiki/Riemann_hypothesis) holds
if and only if `det(A) = O( N! N^(-1/2+epsilon))` for every `epsilon > 0`.

```julia
julia> Riemann(7)
7x7 Riemann{Int64}:
  1  -1   1  -1   1  -1   1
 -1   2  -1  -1   2  -1  -1
 -1  -1   3  -1  -1  -1   3
 -1  -1  -1   4  -1  -1  -1
 -1  -1  -1  -1   5  -1  -1
 -1  -1  -1  -1  -1   6  -1
 -1  -1  -1  -1  -1  -1   7
```

For more details see [F. Roesler (1986)][Roesler1986].

[Roesler1986]: http://www.sciencedirect.com/science/article/pii/0024379586902557 "Friedrich Roesler, Riemann's hypothesis as an eigenvalue problem, Linear Algebra and its Applications, Vol. 81, (1986)"




## `Strang` matrix

A special `SymTridiagonal` matrix named after Gilbert Strang

```julia
julia> Strang(6)
6x6 Strang{Float64}:
  2.0  -1.0   0.0   0.0   0.0   0.0
 -1.0   2.0  -1.0   0.0   0.0   0.0
  0.0  -1.0   2.0  -1.0   0.0   0.0
  0.0   0.0  -1.0   2.0  -1.0   0.0
  0.0   0.0   0.0  -1.0   2.0  -1.0
  0.0   0.0   0.0   0.0  -1.0   2.0
```

## [`Toeplitz`](http://en.wikipedia.org/wiki/Toeplitz_matrix) matrix

Input is a vector of odd length.

```julia
julia> Toeplitz(collect(-4:4))
5x5 Toeplitz{Int64}:
 0  -1  -2  -3  -4
 1   0  -1  -2  -3
 2   1   0  -1  -2
 3   2   1   0  -1
 4   3   2   1   0
```
## [`Vandermonde`](http://en.wikipedia.org/wiki/Vandermonde_matrix) matrix

```julia
julia> a = collect(1.0:5.0)
julia> A = Vandermonde(a)
5×5 Vandermonde{Float64}:
 1.0  1.0   1.0    1.0    1.0
 1.0  2.0   4.0    8.0   16.0
 1.0  3.0   9.0   27.0   81.0
 1.0  4.0  16.0   64.0  256.0
 1.0  5.0  25.0  125.0  625.0
```

Adjoint Vandermonde:
```julia
julia> A'
5×5 LinearAlgebra.Adjoint{Float64,Vandermonde{Float64}}:
 1.0   1.0   1.0    1.0    1.0
 1.0   2.0   3.0    4.0    5.0
 1.0   4.0   9.0   16.0   25.0
 1.0   8.0  27.0   64.0  125.0
 1.0  16.0  81.0  256.0  625.0
```

The backslash overator `\` is overloaded to solve Vandermonde and adjoint Vandermonde systems in ``O(n^2)`` time using the algorithm of [Björck & Pereyra (1970)](https://doi.org/10.2307/2004623
),
```julia
julia> A\a
5-element Array{Float64,1}:
 0.0
 1.0
 0.0
 0.0
 0.0

julia> r2 = A[2,:]
julia> A'\r2
5-element Array{Float64,1}:
 0.0
 1.0
 0.0
 0.0
 0.0
```
