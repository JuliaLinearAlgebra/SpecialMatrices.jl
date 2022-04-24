# kronecker.jl tests

using SpecialMatrices: Kronecker
import SpecialMatrices
using Test: @test, @testset, @test_throws, @inferred

# abstract vector construction
# V = @inferred Vandermonde(5:7)
# @test V == (5:7) .^ (0:2)'

A = rand(2,3) + 1im*rand(2,3)
B = rand(4,5)
v = rand(15,6)
K = @inferred Kronecker(A,B)
M = @inferred Matrix(K)

# Test element
@testset "basics" begin
  @test K[8,15] == A[2,3]*B[4,5]

  @test isassigned(K, 8*15)
  @test !isassigned(K, 8*15+1)
  @test size(K) == (8, 15)
end

# Test full matrix conversions
@testset "convert" begin
  @test Matrix(K') == M'
  @test Matrix(transpose(K)) == transpose(M)
end

#
@testset "vector multiplication" begin
  A = rand(2,3)
  B = rand(4,5)
  for nv in [1 3]
    v = rand(15,nv)
    @test Kronecker(A,B)*v ≈ kron(A,B)*v
  end
end

@testset "inversion" begin
  n = 6;
  A = rand(n,n); B = rand(2*n,2*n)
  @test inv(Kronecker(A,B)) ≈ inv(kron(A,B))
end

@testset "transpose" begin
  @test transpose(K) ≈ transpose(M)
end

@testset "adjoint" begin
  @test adjoint(K) ≈ adjoint(M)
end

@testset "trace" begin
  n = 6;
  A = rand(n,n); B = rand(2*n,2*n)

  @test tr(Kronecker(A,B)) ≈ tr(kron(A,B))
end

#
@testset "solver" begin
  n = 6;
  A = rand(n,n); B = rand(n,n);

  for nv in [1 3]
    v = rand(n^2,nv)
    x = Kronecker(A,B) \ v
    y = kron(A,B) \ v
    @test x ≈ y
  end
end

@testset "eigenvalues" begin
  n = 4;
  A = rand(n,n); B = rand(n,n)

  @test eigvals(Kronecker(A,B)) ≈ eigvals(kron(A,B))
end
