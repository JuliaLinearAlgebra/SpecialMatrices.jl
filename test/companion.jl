n = rand(1:10)
Z = Companion(randn(n))

#Special properties
@test_approx_eq full(inv(Z)) inv(full(Z))

#Matvec product
b = randn(n)
@test_approx_eq Z*b full(Z)*b

m = rand(1:10)
A = randn(m, n)
@test_approx_eq A*Z A*full(Z)
