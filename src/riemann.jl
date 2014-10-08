#### Riemann Matrix 
import Base: size
export Riemann

immutable Riemann{Int} <: AbstractMatrix{Int}
    n::Int
end # immutable

# Define its size

size(A::Riemann, dim::Integer) = A.n
size(A::Riemann)= size(A,1), size(A,1)

# Index into a Riemann
function getindex(A::Riemann,i::Integer,j::Integer)
#    return (i+1)%(j+1)==0 ? i : -1
    if i<=j && (j+1)%(i+1)==0
        return i
    else
        return -1
    end
end # getindex

# Dense version of Riemann
full(A::Riemann) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]

