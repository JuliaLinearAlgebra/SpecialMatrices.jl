 n=rand(1:10)
 @test diag(Matrix(Toeplitz(collect(-n:n)))) == zeros(n+1)
