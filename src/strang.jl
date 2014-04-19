export StrangMatrix

immutable StrangMatrix{T} <: AbstractArray{T, 2}
	n :: Int
end
StrangMatrix(n::Int) = StrangMatrix{Float64}(n)
strang(T, n)=SymTridiagonal(2ones(T, n),-ones(T, n-1))

getindex{T}(S::StrangMatrix{T}, i, j) = getindex(strang(T, S.n), i, j)
size(S::StrangMatrix, r::Int) = r==1 || r==2 ? S.n : throw(ArgumentError("Invalid dimension $r"))
size(S::StrangMatrix) = S.n, S.n
full{T}(S::StrangMatrix{T}) = full(strang(T, S.n))
*{T}(A::VecOrMat{T}, S::StrangMatrix{T}) = A*full(strang(T, S.n))
*{T}(S::StrangMatrix{T}, A::VecOrMat{T}) = full(strang(T, S.n))*A
