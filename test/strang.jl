# strang.jl

using SpecialMatrices: Strang
using LinearAlgebra: diagm, SymTridiagonal
using Test: @test, @testset, @test_throws, @inferred

@testset "strang" begin
    @test_throws ArgumentError Strang(0)

    Z = @inferred Strang(1)
    @test Matrix(Z) == reshape([2.0],(1,1))

    n = 5
    Z = @inferred Strang(Int64, 5)
    @test Z == diagm(0 => 2ones(n), 1 => -ones(n-1), -1 => -ones(n-1))

    @test Z isa SymTridiagonal
end

# no further tests are needed because SymTridiagonal is tested elsewhere
