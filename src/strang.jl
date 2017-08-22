export Strang

immutable Strang{T} <: AbstractArray{T, 2}
    n :: Int
end
Strang(n::Int) = Strang{Float64}(n)
strang(T, n)= n >1 ? SymTridiagonal(2ones(T, n),-ones(T, n-1)) :
              n==1 ? Diagonal([2one(T)]) : error("Invalid dimension ", n)
function getindex{T}(S::Strang{T}, i, j)
    i == j && return 2
    abs(i - j) == 1 && return -1
    0
end
getindex{T}(S::Strang{T}, I...) = getindex(S,ind2sub(size(S),I...)...)
size(S::Strang, r::Int) = r==1 || r==2 ? S.n : throw(ArgumentError("Invalid dimension $r"))
size(S::Strang) = S.n, S.n
full{T}(S::Strang{T}) = full(strang(T, S.n))
*{T}(A::VecOrMat{T}, S::Strang{T}) = A*full(strang(T, S.n))
*{T}(S::Strang{T}, A::VecOrMat{T}) = full(strang(T, S.n))*A
