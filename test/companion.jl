n = rand(1:10)
Z = Companion(randn(n))

#Special properties
@test Matrix(inv(Z)) ≈ inv(Matrix(Z))

#Matvec product
b = randn(n)
@test Z*b ≈ Matrix(Z)*b

m = rand(1:10)
A = randn(m, n)
@test A*Z ≈ A*Matrix(Z)

# Polynomial construction
using Polynomials
p = Polynomial([-1,0,1])
p_c = Companion(p)
v1 = [1,1]
v2 = [-1,1]
@test p_c*v1 == v1
@test p_c*v2 == -v2
