module SpecialMatrices

import Base: A_mul_B!, full, getindex, inv, isassigned, size, *

include("companion.jl") #Companion matrix
include("frobenius.jl") #Frobenius matrix
include("strang.jl") #Strang matrix

end # module
