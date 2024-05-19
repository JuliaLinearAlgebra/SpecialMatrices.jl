using SpecialMatrices: SpecialMatrices
import Aqua
using Test: @testset

@testset "aqua" begin
  if VERSION < v"1.11"
    Aqua.test_all(SpecialMatrices; ambiguities = false)
  else
    Aqua.test_all(SpecialMatrices)
  end
end

#=
todo: 1 ambiguity in julia v1.10
mul!(A::Matrix, B::AbstractMatrix, C::SpecialMatrices.Companion)
 SpecialMatrices/src/companion.jl:89
mul!(C::AbstractMatrix, A::LinearAlgebra.AbstractTriangular, B::AbstractMatrix)
 LinearAlgebra/src/triangular.jl:691
=#
