Z = Strang(1)
@test full(Z) == reshape([2.0],(1,1))

n = rand(1:10)
Z = Strang(n)

#Matvec product
b = randn(n)
@test_approx_eq Z*b full(Z)*b

m = rand(1:10)
A = randn(m, n)
@test_approx_eq A*Z A*full(Z)
