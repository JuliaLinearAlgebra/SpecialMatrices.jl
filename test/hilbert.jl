# test/hilbert.jl

using SpecialMatrices: Hilbert, InverseHilbert
using LinearAlgebra: I, det, ishermitian, isposdef
using Test: @test, @test_throws, @inferred

n = 10
H = @inferred Hilbert(n)
@test H isa Hilbert
@test H[n,n] == 1//(2n-1)

@test size(H) == (n,n)
@test (@inferred ishermitian(H))
@test (@inferred isposdef(H))
@test (@inferred Matrix(H)) isa Matrix

Hi = @inferred inv(H)
@test Hi isa InverseHilbert
@test Hi == InverseHilbert(n)
@test Hi[1,1] == n^2//1

@test size(Hi) == (n,n)
@test (@inferred ishermitian(Hi))
@test (@inferred isposdef(Hi))
@test (@inferred Matrix(Hi)) isa Matrix

Hii = @inferred inv(Hi)
@test Hii == H
@test Hi * H == I

# For "large" n the determinant overflows for Int type:
@test_throws OverflowError det(H)
# But it works fine for Float64 type:
Hf = @inferred Hilbert(Float64, n)
@test (@inferred det(Hf)) isa Float64

Hs = Hilbert(5)
@test (@inferred det(Hs)) ≈ det(Matrix(Hs))

Hb = @inferred InverseHilbert{BigInt}(n)
@test (@inferred getindex(Hb,n,n)) isa BigInt
@test Hb[n,n] isa BigInt
