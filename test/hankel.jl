n = 9
H = Hankel(collect(-n:n))
@test diag(full(H)) == collect(-n:2:n)

@test full(Hankel([1,2,3,4,5])) == [1 2 3;2 3 4; 3 4 5]

@test_throws ArgumentError Hankel([1,2,3,4])
