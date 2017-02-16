n = 9
H = Hankel(collect(-n:n))
@test diag(full(H)) == collect(-n:2:n)

@test full(Hankel([1,2,3,4,5])) == [1 2 3; 2 3 4; 3 4 5]

@test_throws ArgumentError Hankel([1,2,3,4])

H = Hankel([11,12,13,14,15])
@test all(map(x -> Hankel([11,12,13,14,15])[x] == [11 12 13;12 13 14;13 14 15][x], 1:9))

@test size(H) == (3,3)
@test size(H,1) == 3
@test size(H,2) == 3
@test size(H,22) == 1

@test isassigned(H,2,2) == true
@test isassigned(H,22,22) == false

H = Hankel([1,3.2,-0.2,9.1,3.14])
V = Vector([2.1,0.1,-9.9])
@test maxabs(H*V - full(H)*V) < 1e-12

H = Hankel([1,3.2im,-0.2,9.1,3.14])
V = Vector([2.1,0.1,-9.9])
@test maxabs(H*V - full(H)*V) < 1e-12

H = Hankel([1,3.2,-0.2,9.1,3.14])
V = Vector([2.1,0.1im,-9.9])
@test maxabs(H*V - full(H)*V) < 1e-12

H = Hankel([1,3.2im,-0.2,9.1,3.14])
V = Vector([2.1,0.1im,-9.9])
@test maxabs(H*V - full(H)*V) < 1e-12

H = Hankel([1,3.2,-0.2,9.1,3.14])
V = Vector([2.1,0.1,-9.9])
y = zeros(3)
A_mul_B!(y, H, V)
@test maxabs(y - full(H)*V) < 1e-12
