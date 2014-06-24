# Special Matrices

A [Julia](http://julialang.org) package for working with special matrix types.

This package extends the `Base.LinAlg` library with support for special
matrices which are used in linear algebra.

[![Build Status](https://travis-ci.org/jiahao/SpecialMatrices.jl.svg)](https://travis-ci.org/jiahao/SpecialMatrices.jl) [![Coverage Status](https://img.shields.io/coveralls/jiahao/SpecialMatrices.jl.svg)](https://coveralls.io/r/jiahao/SpecialMatrices.jl)


## Currently supported special matrices

## [`Frobenius`](http://en.wikipedia.org/wiki/Frobenius_matrix) matrix

Example

```julia
julia> using SpecialMatrices

julia> F=Frobenius(3, [1.0:3.0]) #Specify subdiagonals of column 3
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

julia> F*Frobenius(2, [5.0:-1.0:2.0]) #Promotes to Matrix
6x6 Array{Float64,2}:
 1.0   0.0  0.0  0.0  0.0  0.0
 0.0   1.0  0.0  0.0  0.0  0.0
 0.0   5.0  1.0  0.0  0.0  0.0
 0.0   9.0  1.0  1.0  0.0  0.0
 0.0  13.0  2.0  0.0  1.0  0.0
 0.0  17.0  3.0  0.0  0.0  1.0

julia> F*[10.0:10.0:60.0]
6-element Array{Float64,1}:
  10.0
  20.0
  30.0
  70.0
 110.0
 150.0
```

## [`Companion`](http://en.wikipedia.org/wiki/Companion_matrix) matrix

```julia
julia> A=Companion([1,2,1])
3x3 Companion{Int64}:
 0  0  -1
 1  0  -2
 0  1  -1
```

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

## `Hankel` matrix

```julia
julia> Hankel([-4:4])
5x5 Hankel{Int64}:
 -4  -3  -2  -1  0
 -3  -2  -1   0  1
 -2  -1   0   1  2
 -1   0   1   2  3
  0   1   2   3  4
```

## `Toeplitz` matrix

```julia
julia> Toeplitz([-4:4])
5x5 Toeplitz{Int64}:
 0  -1  -2  -3  -4
 1   0  -1  -2  -3
 2   1   0  -1  -2
 3   2   1   0  -1
 4   3   2   1   0
```

## `Circulant` matrix

```julia
julia> Circulant([1:4])
4x4 Circulant{Int64}:
 1  4  3  2
 2  1  4  3
 3  2  1  4
 4  3  2  1
```

## `Hilbert` matrix

```julia
julia> full(Hilbert(5))
5x5 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4  1//5
 1//2  1//3  1//4  1//5  1//6
 1//3  1//4  1//5  1//6  1//7
 1//4  1//5  1//6  1//7  1//8
 1//5  1//6  1//7  1//8  1//9
```

## `Vandermonde` matrix

```julia
julia> Vandermonde([1:5])
5x5 Vandermonde{Int64}:
 1  1   1    1    1
 1  2   4    8   16
 1  3   9   27   81
 1  4  16   64  256
 1  5  25  125  625
```
