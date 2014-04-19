n = rand(1:10)
Z = CompanionMatrix(randn(n))

#Special properties
@test_approx_eq full(inv(Z)) inv(full(Z))

#Matvec product
b = randn(n)
@test_approx_eq Z*b full(Z)*b
