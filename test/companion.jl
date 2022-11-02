# test/companion.jl

import LinearAlgebra: mul!

Z = @inferred Companion(1:3) # AbstractVector

n = 9
c = rand(n)
Z = @inferred Companion(c)

@test Z isa Companion{Float64}
@test (@inferred size(Z)) == (n,n)
@test size(Z,1) == n
@test size(Z,3) == 1

@test Z[2,1] == 1
@test (@inferred getindex(Z, 2, n)) == -c[2]

@test isassigned(Z, 1)
@test isassigned(Z, n^2)
@test !isassigned(Z, 0)
@test !isassigned(Z, n^2+1)

# Special properties
Zi = @inferred inv(Z)
Zm = @inferred Matrix(Z)
@test Matrix(Zi) ≈ inv(Zm)

# Matvec product
x = randn(n)
@test Z * x ≈ Zm * x

y = similar(x)
@inferred mul!(y, Z, x)
@test y ≈ Z * x


# matrix * companion
m = 8
B = randn(m, n)
@test B * Z ≈ B * Zm

A = copy(B)
@inferred mul!(A, B, Z)
@test A ≈ B * Zm


# Polynomial construction
using Polynomials
p = Polynomial([-1,0,1])
C = @inferred Companion(p)
v1 = [1,1]
v2 = [-1,1]
@test C * v1 == v1
@test C * v2 == -v2
