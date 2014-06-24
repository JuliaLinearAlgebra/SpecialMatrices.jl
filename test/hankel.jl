n=rand(1:10)
H = Hankel([-n:n])
@test diag(full(H)) == [-n:2:n]
