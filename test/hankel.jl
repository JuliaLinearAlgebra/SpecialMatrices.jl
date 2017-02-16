n = 9
H = Hankel(collect(-n:n))
@test diag(full(H)) == collect(-n:2:n)

@test full(Hankel([1,2,3,4,5])) == [1 2 3;2 3 4; 3 4 5]

@test_throws ArgumentError Hankel([1,2,3,4])

H = Hankel([11,12,13,14,15])
@test all(map(x -> Hankel([11,12,13,14,15])[x] == [11 12 13;12 13 14;13 14 15][x], 1:9))

@test size(H) == (3,3)
@test size(H,1) == 3
@test size(H,2) == 3
@test size(H,22) == 1
