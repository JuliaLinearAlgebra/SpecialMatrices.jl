using SpecialMatrices
using Test
using LinearAlgebra

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
