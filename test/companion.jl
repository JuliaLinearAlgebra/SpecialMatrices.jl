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
