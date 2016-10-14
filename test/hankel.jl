n=rand(1:10)
H = Hankel(collect(-n:n))
@test diag(full(H)) == collect(-n:2:n)
