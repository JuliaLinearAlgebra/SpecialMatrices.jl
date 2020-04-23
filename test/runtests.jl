using SpecialMatrices
using Test
using LinearAlgebra

const files = (
    "companion",
    "frobenius",
    "hankel",
    "hilbert",
    "strang",
    "toeplitz",
    "vandermonde",
)

@testset "SpecialMatrices.jl" begin
    @testset "$(titlecase(f)) matrix" for f in files
        include("$f.jl")
    end
end
