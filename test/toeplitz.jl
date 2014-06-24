 n=rand(1:10)
 @test diag(full(Toeplitz([-n:n]))) == zeros(n+1)
