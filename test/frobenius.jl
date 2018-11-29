m=rand(2:10)
n=rand(1:10)
Z = Frobenius(m, randn(n)) #size m+n

#Special properties
@test Matrix(inv(Z)) ≈ inv(Matrix(Z))

#Matvec product
b = randn(m+n)
@test Z*b ≈ Matrix(Z)*b

#Matmul product
Y1 = Frobenius(m, randn(n)) #Another one of the same column
@test Matrix(Y1*Z) ≈ Matrix(Y1)*Matrix(Z)

n2=rand(1:(m+n-1))
Y2 = Frobenius(m+n-n2, randn(n2)) #Probably not the same column
@test Matrix(Y2*Z) ≈ Matrix(Y2)*Matrix(Z)
