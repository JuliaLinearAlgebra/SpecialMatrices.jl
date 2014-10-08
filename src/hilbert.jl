#The m-by-n Hilbert matrix has matrix elements
# H_{ij} = 1/(i+j-1)
import Base: full
export Hilbert

type Hilbert{T}
    m :: Int
    n :: Int
end
Hilbert(T::Type, m::Integer, n::Integer) = Hilbert{T}(m, n)
Hilbert(m::Integer, n::Integer) = Hilbert{Rational{Int}}(m, n)
Hilbert(n::Integer) = Hilbert(n, n)

function full{T}(H::Hilbert{T})
    Hf = zeros(T, H.m, H.n)
    for i=1:H.m, j=1:H.n
        Hf[i, j] = one(T)/(i+j-1)
    end
    Hf
end

function inv{T}(H::Hilbert{T})
    invH=zeros(T,H.m,H.m)
    for i=1:H.m, j=1:H.m
        invH[i,j]=(-1)^(i+j)*(i+j-1)*binomial(H.m+i-1,H.m-j)*
                  binomial(H.m+j-1,H.m-i)*binomial(i+j-2,i-1)^2
    end
    invH
end

