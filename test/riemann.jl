# riemann.jl tests

using SpecialMatrices: Riemann
using Test: @test, @testset, @test_throws, @inferred

@testset "riemann" begin
    n = 5
    R = @inferred Riemann(n)
    @test size(R) == (n, n)
    @test R[1,1] == 1
    @test R[n,n] == n

    M = @inferred Matrix(R)
    @test Matrix(R') == M'
    @test Matrix(transpose(R)) == transpose(M)
end
