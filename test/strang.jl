# strang.jl

using SpecialMatrices: Strang
using LinearAlgebra: diagm, issymmetric, ishermitian, isposdef
using LinearAlgebra: factorize, SymTridiagonal, Diagonal
using Test: @test, @testset, @test_throws, @inferred

@testset "strang" begin
    @test_throws ArgumentError Strang(0)

    A = @inferred Strang(1)
    @test Matrix(A) == reshape([2.0], (1,1))

    n = 5
    A = @inferred Strang(Int16, 5)
    @test A isa Strang{Int16}
    @test A == diagm(0 => 2ones(n), 1 => -ones(n-1), -1 => -ones(n-1))

    @test (@inferred getindex(A, 1, 2)) == -1
    @test A[begin] == 2
    @test A[end] == 2
    @test A[1,end] == 0

    x = rand(n)
    y = @inferred *(A, x)
    @test y ≈ Matrix(A) * x

    @test transpose(A) === A
    @test adjoint(A) === A
    @test issymmetric(A)
    @test ishermitian(A)
    @test isposdef(A)
end


if VERSION >= v"1.8"

@testset "strang-factor" begin
    # factorize:
    n = 5
    A = Strang(5)
    M = Matrix(A)
    T = SymTridiagonal(M)
    F = factorize(T)
    G = factorize(A)
    H = G.L * Diagonal(G.D) * G.L'
    @test H ≈ M
    @test F.L ≈ G.L
    @test F.D ≈ G.D
end

end
