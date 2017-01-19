module SpecialMatrices

import Base: A_mul_B!, full, getindex, inv, isassigned, size, *, length, +, -,
    Ac_mul_B, A_mul_Bc, At_mul_B, A_mul_Bt, eig, eigfact, Eigen, \, /,
    transpose, ctranspose, copy, conj, conj!, inv, det, logdet, real, convert,
    round

typealias SV StridedVector
typealias SM StridedMatrix
typealias SVM StridedVecOrMat

include("dft.jl") # Discrete Fourier Transform matrices.

include("cauchy.jl") # Cauchy matrix
include("circulant.jl") # Circulant matrix.
include("companion.jl") # Companion matrix
include("frobenius.jl") # Frobenius matrix
include("hankel.jl") # Hankel matrix
include("hilbert.jl") # Hilbert matrix
include("kahan.jl") # Kahan matrix
include("riemann.jl") # Riemann matrix
include("strang.jl") # Strang matrix
include("toeplitz.jl") # Toeplitz matrix
include("vandermonde.jl") # Vandermonde matrix

end # module
