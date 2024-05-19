using SpecialMatrices: SpecialMatrices
using Test: @test, @test_broken, @testset, detect_ambiguities

include("aqua.jl")

const files = (
    "cauchy",
    "companion",
    "frobenius",
    "hilbert",
    "kahan",
    "riemann",
    "strang",
    "vandermonde",
)

@testset "SpecialMatrices.jl" begin
    @testset "$(titlecase(f)) matrix" for f in files
        include("$f.jl")
    end
end

@testset "ambiguities" begin
if VERSION < v"1.11"
    @test_broken isempty(detect_ambiguities(SpecialMatrices)) # see aqua.jl
else
    @test isempty(detect_ambiguities(SpecialMatrices))
end
end
