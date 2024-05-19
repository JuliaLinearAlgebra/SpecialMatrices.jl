# test/kahan.jl

using SpecialMatrices: Kahan
using Test: @test, @inferred

m, n = 3, 5
A = Kahan(m, n, 0.5, 1e-3)

@test (@inferred size(A)) == (m,n)
@test (@inferred size(A,1)) == m
@test (@inferred size(A,2)) == n

@test A[1,1] == 1

M = @inferred Matrix(A)
@test M isa Matrix
@test size(M) == size(A)
