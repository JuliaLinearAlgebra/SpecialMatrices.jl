C = @inferred Cauchy(3)
@test 1 ./ C ≈ [2 3 4; 3 4 5; 4 5 6]
@test size(C) == (3,3)
Cf = @inferred Matrix(C)
@test Cf isa Matrix

@test C[1,1] == 1
@inferred getindex(C, 1, 1)

x = 1:3
y = [2im 10f0]
C = @inferred Cauchy(x, y)
@test 1 ./ C ≈ [1+2im 11; 2+2im 12; 3+2im 13]
