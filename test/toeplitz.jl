 n=rand(1:10)
 @test diag(full(Toeplitz(collect(-n:n)))) == zeros(n+1)
