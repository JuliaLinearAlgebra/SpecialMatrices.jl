# test/frobenius.jl

using SpecialMatrices: Frobenius
using Test: @test, @test_throws, @inferred
using LinearAlgebra: I
import LinearAlgebra: mul!

n = 7
j = 4
function test_frobenius(col)
    A = @inferred Frobenius(j, col) # AbstractVector
    @test A isa Frobenius

    # size
    @test size(A) == (n,n)
    @test size(A,1) == n
    @test size(A,2) == n
    @test size(A,3) == 1

    # inv
    Ai = @inferred inv(A)
    @test Ai * A == I

    # getindex
    @test A[1,1] == 1
    @test (@inferred getindex(A, 2, 2)) == 1
    @test (@inferred getindex(A, j+2, j)) == col[2]
    @test_throws BoundsError getindex(A, 8, 8)

    # Matrix
    M = @inferred Matrix(A)
    @test M isa Matrix

    # mul! *
    y = Vector{Float32}(undef, n)
    x = 10 * (1:n)
    @inferred mul!(y, A, x)
    @test y == Matrix(A) * x
    @test y == A * x

    x = randn(n)
    @inferred mul!(y, A, x)
    @test y ≈ Matrix(A) * x
    @test y ≈ A * x

    # Matrix-Matrix
    B = Frobenius(j, randn(size(col))) # another one of the same column
#   C = @inferred *(A, B) # A * B (NOT type stable)
    C = *(A, B) # A * B
    @test C isa Frobenius
    @test Matrix(C) ≈ Matrix(A) * Matrix(B)

    B = Frobenius(j-1, randn(1 .+ size(col))) # different column
    C = *(A, B) # A * B
    @test C isa Matrix
    @test Matrix(C) ≈ Matrix(A) * Matrix(B)
end

for col in (1:3, 2.0:4.0) # test both Int and Float32 types
   test_frobenius(col)
end

m = rand(2:10)
n = rand(1:10)
Z = Frobenius(m, randn(n)) # size m+n

#Special properties
@test Matrix(inv(Z)) ≈ inv(Matrix(Z))

#Matvec product
b = randn(m+n)
@test Z*b ≈ Matrix(Z)*b

#Matmul product
Y1 = Frobenius(m, randn(n)) #Another one of the same column
@test Matrix(Y1*Z) ≈ Matrix(Y1)*Matrix(Z)

n2 = rand(1:(m+n-1))
Y2 = Frobenius(m+n-n2, randn(n2)) #Probably not the same column
@test Matrix(Y2*Z) ≈ Matrix(Y2)*Matrix(Z)
