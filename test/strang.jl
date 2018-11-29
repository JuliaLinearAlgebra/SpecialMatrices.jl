Z = Strang(1)
@test Matrix(Z) == reshape([2.0],(1,1))

n = rand(1:10)
Z = Strang(n)

for i in 1:n, j in 1:n
    i==j && @test Z[i,j] == 2
    abs(i-j)==1 && @test Z[i,j] == -1
    abs(i-j)>1 && @test Z[i,j] == 0
end

A = Strang(10)
u = ones(10)
@test A*u == [1.0;zeros(8);1.0]

#Matvec product
b = randn(n)
@test Z*b ≈ Matrix(Z)*b

m = rand(1:10)
A = randn(m, n)
@test A*Z ≈ A*Matrix(Z)
