N = 5
x = randn(N)
X = reshape(x, N, 1)
U = DFT{Float64}(N)

function check_dot(N=N, x=x, X=X, U=U)
    all((U * x) .== (U * X))
end
@test check_dot()

function check_inplace_dot(N=N, x=x, X=X, U=U)
    x = convert(Vector{Complex{Float64}}, x)
    X = convert(Matrix{Complex{Float64}}, X)
    y = zeros(Complex{Float64}, N)
    Y = zeros(Complex{Float64}, N, 1)
    all(A_mul_B!(y, U, x) == U * x) && all(A_mul_B!(Y, U, X) == U * X) &&
        all(A_mul_B!(Y', X', U) == X'U)
end
@test check_inplace_dot()