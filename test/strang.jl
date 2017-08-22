Z = Strang(1)
@test full(Z) == reshape([2.0],(1,1))

n = rand(1:10)
Z = Strang(n)

for i in 1:n, j in 1:n
    i==j && @test Z[i,j] == 2
    abs(i-j)==1 && @test Z[i,j] == -1
    abs(i-j)>1 && @test Z[i,j] == 0
end

A = Strang(10)
u = ones(10)
v = similar(u)
A_mul_B!(v,A,u)
@test v == [1.0;zeros(8);1.0]

#Matvec product
b = randn(n)
@test_approx_eq Z*b full(Z)*b

m = rand(1:10)
A = randn(m, n)
@test_approx_eq A*Z A*full(Z)
