using SpecialMatrices: SpecialMatrices
using Test: @test_broken, @testset, detect_ambiguities

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
    @test_broken isempty(detect_ambiguities(SpecialMatrices)) # see aqua.jl
end
