# strang.jl

using SpecialMatrices: Strang
using LinearAlgebra: diagm
using Test: @test, @testset, @test_throws, @inferred

@testset "strang" begin
    @test_throws ArgumentError Strang(0)

    A = @inferred Strang(1)
    @test Matrix(A) == reshape([2.0], (1,1))

    n = 5
    A = @inferred Strang(Int16, 5)
    @test A isa Strang{Int16}
    @test A == diagm(0 => 2ones(n), 1 => -ones(n-1), -1 => -ones(n-1))

    x = rand(n)
    @test A * x ≈ Matrix(A) * x
end
