m=rand(1:10)
n=rand(1:10)
Z = FrobeniusMatrix(m, randn(n))
@test_approx_eq full(inv(Z)) inv(full(Z))

Y1 = FrobeniusMatrix(m, randn(n)) #Another one of the same column
@test_approx_eq full(Y1*Z) full(Y1)*full(Z)

n2=rand(1:10)
Y2 = FrobeniusMatrix(m+n-n2, randn(n2)) #Probably not the same column
@test_approx_eq full(Y2*Z) full(Y2)*full(Z) #XXX NOT WORKING
