module SpecialMatrices

using Compat
using LinearAlgebra
using FFTW
using Polynomials

import LinearAlgebra: Matrix, inv, det, ishermitian, isposdef, mul!

import Base: getindex, isassigned, size, *
import Base.\

include("cauchy.jl") #Cauchy matrix
include("companion.jl") #Companion matrix
include("frobenius.jl") #Frobenius matrix
include("hankel.jl") #Hankel matrix
include("hilbert.jl") #Hilbert matrix
include("kahan.jl") #Kahan matrix
include("riemann.jl") #Riemann matrix
include("strang.jl") #Strang matrix
include("toeplitz.jl") #Toeplitz matrix, Circulant matrix
include("vandermonde.jl") #Vandermonde matrix

end # module
