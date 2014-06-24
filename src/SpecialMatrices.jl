module SpecialMatrices

import Base: A_mul_B!, full, getindex, inv, isassigned, size, *

include("companion.jl") #Companion matrix
include("frobenius.jl") #Frobenius matrix
include("strang.jl") #Strang matrix
include("hankel.jl") #Hankel matrix
include("hilbert.jl") #Hilbert matrix
include("toeplitz.jl") #Toeplitz matrix
include("vandermonde.jl") #Vandermonde matrix

end # module
