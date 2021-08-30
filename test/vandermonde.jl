# vandermonde.jl tests

using SpecialMatrices: Vandermonde
import SpecialMatrices # dvand!, pvand!
using Test: @test, @testset, @test_throws, @inferred

# abstract vector construction
V = @inferred Vandermonde(5:7)
@test V == (5:7) .^ (0:2)'

a = [1,2,3+1im,8,5]
V = @inferred Vandermonde(a)
M = @inferred Matrix(V)

# Test element
@testset "basics" begin
    @test V[4,5] == a[4]^4

    @test isassigned(V, 2)
    @test !isassigned(V, 26)
    @test size(V) == (5, 5)
end

# Test full matrix conversions
@testset "convert" begin
    @test Matrix(V') == M'
    @test Matrix(transpose(V)) == transpose(M)
end

@testset "solve" begin
    # Test solving with vector and matrix rhs
    y = [1im,1,5,0,2]
    Y = [y 2*y]
    for rhs in [y, Y]
        # Test solution
        @test isapprox(V \ rhs, M \ rhs)
        @test isapprox(V' \ rhs, M' \ rhs)
        @test isapprox(rhs' / V, rhs' / M)
        @test isapprox(transpose(V) \ rhs, transpose(M) \ rhs)

        # Check that overloading works
        x = zero(M \ rhs)
        copyto!(x, rhs)
        SpecialMatrices.dvand!(a, x)
#       @which V \ rhs
        @test V \ rhs == x

        copyto!(x, rhs)
        SpecialMatrices.pvand!(a', x)
#       @which V' \ rhs
        @test V' \ rhs == x

#       @which rhs' / V
        @test rhs' / V == x'

        copyto!(x, rhs)
        SpecialMatrices.pvand!(a, x)
#       @which transpose(V) \ rhs
        @test transpose(V) \ rhs == x
    end
end

# Test dimension errors
@testset "dims" begin
    rhs = zeros(2,2)
    @test_throws DimensionMismatch V \ rhs
    @test_throws DimensionMismatch V' \ rhs
end
