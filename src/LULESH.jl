module LULESH

using CEnum
using CUDA
using MPI

abstract type AbstractDomain end
# Device vector types
const VD{T}  = Union{Vector{T}, CuVector{T}} where T

# Index type
const IndexT = Int

include("types.jl")
include("domain.jl")

end
