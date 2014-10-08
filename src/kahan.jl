#### Kahan matrix

export Kahan

immutable Kahan{T<:Number,T1<:Number} <: AbstractMatrix{T}
    m::Int # dimension
    n::Int # dimension
    theta::T # angle
    pert::T1 # perturbation is pert*eps()
end # immutable

# Define its size
size(A::Kahan, r::Int) = r==1 ? A.m : A.n
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

# Dense version of Kahan
full(A::Kahan) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]



