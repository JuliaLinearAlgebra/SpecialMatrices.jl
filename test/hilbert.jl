# test/hilbert.jl

n = 10
H = @inferred Hilbert(n)
@test H isa Hilbert
@test H[n,n] == 1//(2n-1)

@test size(H) == (n,n)
@test ishermitian(H)
@test (@inferred Matrix(H)) isa Matrix

Hi = @inferred inv(H)
@test Hi isa InverseHilbert
@test Hi == InverseHilbert(n)
@test Hi[1,1] == n^2//1

@test size(Hi) == (n,n)
@test ishermitian(Hi)
@test (@inferred Matrix(Hi)) isa Matrix

if VERSION >= v"1.6" # This test fails with 1.0 so exclude it there.
   @test isposdef(H)
   @test isposdef(Hi)
end

Hii = @inferred inv(Hi)
@test Hii == H
@test Hi * H == I

# For "large" n the determinant overflows for Int type:
@test_throws OverflowError det(H)
# But it works fine for Float64 type:
Hf = Hilbert(Float64, n)
@test (@inferred det(Hf)) isa Float64

Hs = Hilbert(5)
@test (@inferred det(Hs)) ≈ det(Matrix(Hs))

Hb = @inferred Hilbert(BigInt, n)
@test Hb[n,n] isa Rational{BigInt}
