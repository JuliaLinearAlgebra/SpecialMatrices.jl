n=rand(1:10)
H = Hankel(collect(-n:n))
@test diag(Matrix(H)) == collect(-n:2:n)
