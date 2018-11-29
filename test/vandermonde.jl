using Compat
using Compat.Test
using SpecialMatrices

a = [1,2,3+1im,8,5]
V=Vandermonde(a)

# Test element
@test V[4,5] == a[4]^4

# Test full matrix conversions
M = Matrix(V)
@test Matrix(V') == M'
@test Matrix(transpose(V)) == transpose(M)

# Test solving with vector and matrix rhs
y = [1im,1,5,0,2]
Y = [y 2*y]
for rhs=[y, Y]
    # Test solution
    @test isapprox(V\rhs, M\rhs)
    @test isapprox(V'\rhs, M'\rhs)
    @test isapprox(rhs'/V, rhs'/M)
    @test isapprox(transpose(V)\rhs, transpose(M)\rhs)
    # Check that overloading works
    x = zero(M\rhs)
    copyto!(x, rhs)
    SpecialMatrices.dvand!(a, x)
    @test V\rhs==x
    if VERSION >= v"0.7.0"
        copyto!(x, rhs)
        SpecialMatrices.pvand!(a', x)
        @test V'\rhs==x
        @test rhs'/V==x'
        copyto!(x, rhs)
        SpecialMatrices.pvand!(a, x)
        @test transpose(V)\rhs==x
    end
end

# Test dimension errors
rhs = zeros(2,2)
try
    V\rhs
catch e
    @test typeof(e)==DimensionMismatch
end
try
    V'\rhs
catch e
    @test typeof(e)==DimensionMismatch
end
