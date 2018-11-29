using SpecialMatrices
using Compat.Test

a = [1,2,3+1im,8,5]
V=Vandermonde(a)

# Test element
@test V[4,5] == a[4]^4

# Test full matrix conversions
M = Matrix(V)
@test Matrix(V') == M'
@test Matrix(transpose(V)) == transpose(M)

# Test solving
y = [1im,1,5,0,2]
@test isapprox(V\y, M\y)
@test isapprox(V'\y, M'\y)
@test isapprox(y'/V, y'/M)
@test isapprox(transpose(V)\y, transpose(M)\y)

# Test that overloading works
x2 = zero(M\y)
SpecialMatrices.dvand!(x2, a, y)
@test V\y==x2
if VERSION >= v"0.7.0"
    SpecialMatrices.pvand!(x2, a', y)
    @test V'\y==x2
    @test y'/V==x2'
    SpecialMatrices.pvand!(x2, a, y)
    @test transpose(V)\y==x2    
end

