function check_vector()
    N = 5
    x = randn(N)
    U = DFT(N)
    all((U * x) .== (U * reshape(x, N, 1)))
end
@test check_vector()