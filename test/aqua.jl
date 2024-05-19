using SpecialMatrices: SpecialMatrices
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(SpecialMatrices; ambiguities = false)
end

#=
todo: 1 ambiguity
mul!(A::Matrix, B::AbstractMatrix, C::SpecialMatrices.Companion)
 SpecialMatrices/src/companion.jl:89
mul!(C::AbstractMatrix, A::LinearAlgebra.AbstractTriangular, B::AbstractMatrix)
 LinearAlgebra/src/triangular.jl:691
=#
