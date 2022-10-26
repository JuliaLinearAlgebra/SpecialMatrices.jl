C = @inferred Cauchy(3)
@test 1 ./ C ≈ [2 3 4; 3 4 5; 4 5 6]
@test size(C) == (3,3)
Cf = @inferred Matrix(C)
@test Cf isa Matrix{<:Rational}
@inferred getindex(C, 1, 1)
@test C[1,1] == 1//2

x = (1, 2, 3) # test iterator
y = [2im 10f0] # test row vector
C = @inferred Cauchy(x, y)
Cf = @inferred Matrix(C)
@test 1 ./ C ≈ [1+2im 11; 2+2im 12; 3+2im 13]
@inferred getindex(C, 1, 1)
@test C[1,1] ≈ 1 / (1 + 2im)

x = (1, 2.0) # inconsistent element types
y = 1:3
@test_throws ArgumentError Cauchy(x, y)
