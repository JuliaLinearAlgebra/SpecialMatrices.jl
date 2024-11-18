var documenterSearchIndex = {"docs":
[{"location":"methods/#Methods-list","page":"Methods","title":"Methods list","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"","category":"page"},{"location":"methods/#Methods-usage","page":"Methods","title":"Methods usage","text":"","category":"section"},{"location":"methods/","page":"Methods","title":"Methods","text":"Modules = [SpecialMatrices]","category":"page"},{"location":"methods/#SpecialMatrices.Cauchy","page":"Methods","title":"SpecialMatrices.Cauchy","text":"Cauchy matrix\n\nCauchy(x,y)[i,j]=1/(x[i]+y[j])\nCauchy(x)=Cauchy(x,x)\ncauchy(k::Number)=Cauchy(collect(1:k))\n\njulia> Cauchy([1,2,3],[3,4,5])\n3x3 Cauchy{Int64}:\n 0.25      0.2       0.166667\n 0.2       0.166667  0.142857\n 0.166667  0.142857  0.125\n\njulia> Cauchy([1,2,3])\n3x3 Cauchy{Int64}:\n 0.5       0.333333  0.25\n 0.333333  0.25      0.2\n 0.25      0.2       0.166667\n\njulia> Cauchy(pi)\n3x3 Cauchy{Float64}:\n 0.5       0.333333  0.25\n 0.333333  0.25      0.2\n 0.25      0.2       0.166667\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Companion","page":"Methods","title":"SpecialMatrices.Companion","text":"Companion matrix\n\njulia> A=Companion([3,2,1])\n3x3 Companion{Int64}:\n 0  0  -3\n 1  0  -2\n 0  1  -1\n\nAlso, directly from a polynomial:\n\njulia> using Polynomials\n\njulia> P=Polynomial([2.0,3,4,5])\nPolynomial(2 + 3x + 4x^2 + 5x^3)\n\njulia> C=Companion(P)\n3×3 Companion{Float64}:\n 0.0  0.0  -0.4\n 1.0  0.0  -0.6\n 0.0  1.0  -0.8\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Frobenius","page":"Methods","title":"SpecialMatrices.Frobenius","text":"Frobenius matrix\n\nFrobenius matrices or Gaussian elimination matrices of the form\n\n[ 1 0 ...     0 ]\n[ 0 1 ...     0 ]\n[ .........     ]\n[ ... 1 ...     ]\n[ ... c1 1 ...  ]\n[ ... c2 0 1 ...]\n[ ............. ]\n[ ... ck ...   1]\n\ni.e. an identity matrix with nonzero subdiagonal elements along a single column.\n\n In this implementation, the subdiagonal of the nonzero column is stored as a\n dense vector, so that the size can be inferred automatically as j+k where j is\n the index of the column and k is the number of subdiagonal elements.\n\njulia> using SpecialMatrices\n\njulia> F=Frobenius(3, [1.0,2.0,3.0]) #Specify subdiagonals of column 3\n6x6 Frobenius{Float64}:\n 1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  1.0  1.0  0.0  0.0\n 0.0  0.0  2.0  0.0  1.0  0.0\n 0.0  0.0  3.0  0.0  0.0  1.0\n\njulia> inv(F) #Special form of inverse\n6x6 Frobenius{Float64}:\n 1.0  0.0   0.0  0.0  0.0  0.0\n 0.0  1.0   0.0  0.0  0.0  0.0\n 0.0  0.0   1.0  0.0  0.0  0.0\n 0.0  0.0  -1.0  1.0  0.0  0.0\n 0.0  0.0  -2.0  0.0  1.0  0.0\n 0.0  0.0  -3.0  0.0  0.0  1.0\n\njulia> F*F #Special form preserved if the same column has the subdiagonals\n6x6 Frobenius{Float64}:\n 1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  1.0  0.0  0.0  0.0  0.0\n 0.0  0.0  1.0  0.0  0.0  0.0\n 0.0  0.0  2.0  1.0  0.0  0.0\n 0.0  0.0  4.0  0.0  1.0  0.0\n 0.0  0.0  6.0  0.0  0.0  1.0\n\njulia> F*Frobenius(2, [5.0,4.0,3.0,2.0]) #Promotes to Matrix\n6x6 Array{Float64,2}:\n 1.0   0.0  0.0  0.0  0.0  0.0\n 0.0   1.0  0.0  0.0  0.0  0.0\n 0.0   5.0  1.0  0.0  0.0  0.0\n 0.0   9.0  1.0  1.0  0.0  0.0\n 0.0  13.0  2.0  0.0  1.0  0.0\n 0.0  17.0  3.0  0.0  0.0  1.0\n\njulia> F*[10.0,20,30,40,50,60.0]\n6-element Array{Float64,1}:\n  10.0\n  20.0\n  30.0\n  70.0\n 110.0\n 150.0\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Hilbert","page":"Methods","title":"SpecialMatrices.Hilbert","text":"Hilbert matrix\n\njulia> A=Hilbert(5)\nHilbert{Rational{Int64}}(5,5)\n\njulia> Matrix(A)\n5x5 Array{Rational{Int64},2}:\n 1//1  1//2  1//3  1//4  1//5\n 1//2  1//3  1//4  1//5  1//6\n 1//3  1//4  1//5  1//6  1//7\n 1//4  1//5  1//6  1//7  1//8\n 1//5  1//6  1//7  1//8  1//9\n\njulia> Matrix(Hilbert(5))\n5x5 Array{Rational{Int64},2}:\n 1//1  1//2  1//3  1//4  1//5\n 1//2  1//3  1//4  1//5  1//6\n 1//3  1//4  1//5  1//6  1//7\n 1//4  1//5  1//6  1//7  1//8\n 1//5  1//6  1//7  1//8  1//9\n\nInverses are also integer matrices:\n\njulia> inv(A)\n5x5 Array{Rational{Int64},2}:\n    25//1    -300//1     1050//1    -1400//1     630//1\n  -300//1    4800//1   -18900//1    26880//1  -12600//1\n  1050//1  -18900//1    79380//1  -117600//1   56700//1\n -1400//1   26880//1  -117600//1   179200//1  -88200//1\n   630//1  -12600//1    56700//1   -88200//1   44100//1\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Kahan","page":"Methods","title":"SpecialMatrices.Kahan","text":"Kahan matrix\n\njulia> A=Kahan(5,5,1,35)\n5x5 Kahan{Int64,Int64}:\n 1.0  -0.540302  -0.540302  -0.540302  -0.540302\n 0.0   0.841471  -0.454649  -0.454649  -0.454649\n 0.0   0.0        0.708073  -0.382574  -0.382574\n 0.0   0.0        0.0        0.595823  -0.321925\n 0.0   0.0        0.0        0.0        0.501368\n\njulia> A=Kahan(5,3,0.5,0)\n5x3 Kahan{Float64,Int64}:\n 1.0  -0.877583  -0.877583\n 0.0   0.479426  -0.420735\n 0.0   0.0        0.229849\n 0.0   0.0        0.0\n 0.0   0.0        0.0\n\njulia> A=Kahan(3,5,0.5,1e-3)\n3x5 Kahan{Float64,Float64}:\n 1.0  -0.877583  -0.877583  -0.877583  -0.877583\n 0.0   0.479426  -0.420735  -0.420735  -0.420735\n 0.0   0.0        0.229849  -0.201711  -0.201711\n\nFor more details see: N. Higham, A Survey of Condition Number Estimation for Triangular Matrices, SIMAX, Vol. 29, No. 4, pp. 588, 1987, http://eprints.ma.man.ac.uk/695/01/covered/MIMSep200710.pdf\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Riemann","page":"Methods","title":"SpecialMatrices.Riemann","text":"Riemann(N::Int)\n\nConstruct N × N Riemann matrix, defined as A = B[2:N+1, 2:N+1], where B[i,j] = i-1 if i divides j, and -1 otherwise. The Riemann hypothesis holds if and only if det(A) = O( N! N^(-1/2+epsilon)) for every epsilon > 0.\n\njulia> Riemann(7)\n7x7 Riemann{Int64}:\n  1  -1   1  -1   1  -1   1\n -1   2  -1  -1   2  -1  -1\n -1  -1   3  -1  -1  -1   3\n -1  -1  -1   4  -1  -1  -1\n -1  -1  -1  -1   5  -1  -1\n -1  -1  -1  -1  -1   6  -1\n -1  -1  -1  -1  -1  -1   7\n\nFor more details see Friedrich Roesler, Riemann's hypothesis as an eigenvalue problem, Linear Algebra and its Applications, Vol. 81, (1986) http://www.sciencedirect.com/science/article/pii/0024379586902557\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Strang","page":"Methods","title":"SpecialMatrices.Strang","text":"Strang matrix\n\nA special SymTridiagonal matrix named after Gilbert Strang\n\njulia> Strang(6)\n6x6 Strang{Float64}:\n  2.0  -1.0   0.0   0.0   0.0   0.0\n -1.0   2.0  -1.0   0.0   0.0   0.0\n  0.0  -1.0   2.0  -1.0   0.0   0.0\n  0.0   0.0  -1.0   2.0  -1.0   0.0\n  0.0   0.0   0.0  -1.0   2.0  -1.0\n  0.0   0.0   0.0   0.0  -1.0   2.0\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.Vandermonde","page":"Methods","title":"SpecialMatrices.Vandermonde","text":"Vandermonde(c::AbstractVector)\n\nCreate a \"lazy\" n × n Vandermonde matrix where A_ij = c_i^j-1, requiring only O(n) storage for the vector c.\n\nThe transpose and adjoint operations are also lazy.\n\njulia> a = 1:5; A = Vandermonde(a)\n5×5 Vandermonde{Int64}:\n 1  1   1    1    1\n 1  2   4    8   16\n 1  3   9   27   81\n 1  4  16   64  256\n 1  5  25  125  625\n\njulia> A'\n5×5 adjoint(::Vandermonde{Int64}) with eltype Int64:\n 1   1   1    1    1\n 1   2   3    4    5\n 1   4   9   16   25\n 1   8  27   64  125\n 1  16  81  256  625\n\nThe backslash operator \\ is overloaded to solve Vandermonde and adjoint Vandermonde systems in O(n^2) time using the algorithm of Björck & Pereyra (1970), https://doi.org/10.2307/2004623.\n\njulia> A \\ a\n5-element Vector{Float64}:\n 0.0\n 1.0\n 0.0\n 0.0\n 0.0\n\njulia> A' \\ A[2,:]\n5-element Vector{Float64}:\n 0.0\n 1.0\n 0.0\n 0.0\n 0.0\n\n\n\n\n\n","category":"type"},{"location":"methods/#SpecialMatrices.dvand!-Tuple{Any, Any}","page":"Methods","title":"SpecialMatrices.dvand!","text":"dvand!(a, b) -> b\n\nSolve system A*x = b in-place.\n\nA is Vandermonde matrix A_ij = a_i^j-1.\n\nAlgorithm by Bjorck & Pereyra, Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903, https://doi.org/10.2307/2004623\n\n\n\n\n\n","category":"method"},{"location":"methods/#SpecialMatrices.pvand!-Tuple{Any, Any}","page":"Methods","title":"SpecialMatrices.pvand!","text":"pvand!(a, b) -> b\n\nSolve system A^T*x = b in-place.\n\nA^T is transpose of Vandermonde matrix A_ij = a_i^j-1.\n\nAlgorithm by Bjorck & Pereyra, Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903, https://doi.org/10.2307/2004623\n\n\n\n\n\n","category":"method"},{"location":"methods/#SpecialMatrices.vandtype-Tuple{Type, Type}","page":"Methods","title":"SpecialMatrices.vandtype","text":"vandtype(T1::Type, T2::Type)\n\nDetermine the return type of Vandermonde{T1} \\ Vector{T2}.\n\n\n\n\n\n","category":"method"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"EditURL = \"https://github.com/JuliaMatrices/SpecialMatrices.jl/blob/master/docs/lit/examples/1-overview.jl\"","category":"page"},{"location":"examples/1-overview/#overview","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"","category":"section"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"This page illustrates the Julia package SpecialMatrices.","category":"page"},{"location":"examples/1-overview/#Setup","page":"SpecialMatrices overview","title":"Setup","text":"","category":"section"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"Packages needed here.","category":"page"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"using SpecialMatrices","category":"page"},{"location":"examples/1-overview/#Cauchy","page":"SpecialMatrices overview","title":"Cauchy","text":"","category":"section"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"Cauchy(collect(1:3), collect(2:4))","category":"page"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"","category":"page"},{"location":"examples/1-overview/","page":"SpecialMatrices overview","title":"SpecialMatrices overview","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SpecialMatrices","category":"page"},{"location":"#SpecialMatrices.jl-Documentation","page":"Home","title":"SpecialMatrices.jl Documentation","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This Julia module exports methods for defining special matrices that are used in linear algebra. Every special matrix has its own type. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(This documentation is under construction).","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the package README for details.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Examples.","category":"page"}]
}
