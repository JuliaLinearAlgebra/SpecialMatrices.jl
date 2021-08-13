#---------------------------------------------------------
# # [SpecialMatrices overview](@id 1-overview)
#---------------------------------------------------------

# This page illustrates the Julia package
# [`SpecialMatrices`](https://github.com/JuliaMatrices/SpecialMatrices.jl).

# ### Setup

# Packages needed here.

using SpecialMatrices

# ### Cauchy

Cauchy(x::AbstractArray, y::AbstractArray) = Cauchy(collect(x), collect(y))

Cauchy(1:3, 2:4)
