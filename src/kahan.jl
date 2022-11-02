#### Kahan matrix

export Kahan

"""
[`Kahan` matrix](http://math.nist.gov/MatrixMarket/data/MMDELI/kahan/kahan.html)

```jldoctest
julia> A = Kahan(5,5,1,35)
5×5 Kahan{Int64, Int64}:
 1.0  -0.540302  -0.540302  -0.540302  -0.540302
 0.0   0.841471  -0.454649  -0.454649  -0.454649
 0.0   0.0        0.708073  -0.382574  -0.382574
 0.0   0.0        0.0        0.595823  -0.321925
 0.0   0.0        0.0        0.0        0.501368

julia> A = Kahan(5,3,0.5,0)
5×3 Kahan{Float64, Int64}:
 1.0  -0.877583  -0.877583
 0.0   0.479426  -0.420735
 0.0   0.0        0.229849
 0.0   0.0        0.0
 0.0   0.0        0.0

julia> A = Kahan(3,5,0.5,1e-3)
3×5 Kahan{Float64, Float64}:
 1.0  -0.877583  -0.877583  -0.877583  -0.877583
 0.0   0.479426  -0.420735  -0.420735  -0.420735
 0.0   0.0        0.229849  -0.201711  -0.201711
```

For more details see:
N. Higham,
"A Survey of Condition Number Estimation for Triangular Matrices,"
SIMAX, Vol. 29, No. 4, pp. 588, 1987,
https://doi.org/10.1137/1029112
"""
struct Kahan{T<:Number,T1<:Number} <: AbstractMatrix{T}
    m::Int # dimension
    n::Int # dimension
    theta::T # angle
    pert::T1 # perturbation is pert*eps()
end # immutable

# Define its size
size(A::Kahan, r::Int) = r==1 ? A.m : r==2 ? A.n : 1
size(A::Kahan) = A.m, A.n

# Index into a Kahan
function getindex(A::Kahan,i::Integer,j::Integer)
    m=minimum(size(A))
    t=tan(A.theta)
    c=1.0/sqrt(1.0+t^2)
    s=t*c
    if i>m; return 0.0
    elseif i>j; return 0.0
    elseif i==j; return s^(i-1)+A.pert*eps()*(m-i+1)
    else return -c*s^(i-1)
    end
end # getindex

# Dense version of Kahan (in Core)
#Matrix(A::Kahan) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]
