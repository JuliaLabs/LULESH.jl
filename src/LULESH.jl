module LULESH

using CEnum
using MPI
using Parameters
using Printf
using Random

abstract type AbstractDomain end
# Device vector types
# const VD{T} = Union{Vector{T}, CuVector{T}} where T
const VD{T} = Vector{T} where T

# Index type
const IndexT = Int
const MAX_FIELDS_PER_MPI_COMM = 6

struct LuleshProblem
    num_iters::Int
    structured::Bool
    nx::Int
    nr::Int
    balance::Int
    cost::Int
    col::Int
    row::Int
    plane::Int
    side::Int
    devicetype
    floattype
    comm::Union{MPI.Comm, Nothing}
    function LuleshProblem(num_iters, structured, nx, nr, balance, cost, devicetype, floattype, comm)
        col, row, plane, side = InitMeshDecomp(comm)
        return new(num_iters, structured, nx, nr, balance, cost, col, row, plane, side, devicetype, floattype, comm)
    end
end

include("bc.jl")
include("types.jl")
include("domain.jl")
include("mpi.jl")
include("utils.jl")

export printUsage, IndexT, Domain, LuleshProblem, getMyRank, getNumRanks, getWtime,
       lagrangeLeapFrog, comm_max, timeIncrement!

function InitMeshDecomp(comm)
    # Assume cube processor layout for now
    numRanks = getNumRanks(comm)
    myRank = getMyRank(comm)
    testProcs = floor(cbrt(numRanks+0.5))
    @show testProcs, numRanks
    if (testProcs*testProcs*testProcs != numRanks)
        error("Num processors must be a cube of an integer (1, 8, 27, ...)")
    end
    @show typeof(testProcs)
    # TODO This is not good for the padding
    if (MAX_FIELDS_PER_MPI_COMM > getCacheCoherencePad(testProcs))
        error("corner element comm buffers too small. MAX_FIELDS_PER_MPI_COMM > CACHE_COHERENCE_PAD_REAL ($MAX_FIELDS_PER_MPI_COMM > $(getCacheCoherencePad(testProcs))")
    end
    testProcs = convert(Int, testProcs)
    dx = convert(Int, testProcs)
    dy = convert(Int, testProcs)
    dz = convert(Int, testProcs)

    # temporary test
    if dx*dy*dz != numRanks
        error("error -- must have as many domains as procs\n") ;
    end
    remainder = (dx*dy*dz) % numRanks
    myDom = Int(0)
    if (myRank < remainder)
        myDom = myRank*( 1+ (div(dx*dy*dz, numRanks)))
    else
      myDom = remainder*( 1+ (div(dx*dy*dz, numRanks))) +
         (myRank - remainder)*(div(dx*dy*dz, numRanks))
    end
   @show typeof(myDom)
   col = myDom % dx
   row = div(myDom, dx) % dy
   plane = div(myDom, (dx*dy))
   side = testProcs

   return col, row, plane, side
end
export InitMeshDecomp
end # module
