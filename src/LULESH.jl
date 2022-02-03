module LULESH

using CEnum
using MPI
using Parameters
using Printf
using Random
using StaticArrays
using LinearAlgebra

abstract type AbstractDomain end
# Device vector types
# const VD{T} = Union{Vector{T}, CuVector{T}} where T
const VD{T} = Vector{T} where T

# Index type
const IndexT = Int
# MPI Message Tags
const MSG_COMM_SBN      = 1024
const MSG_SYNC_POS_VEL  = 2048
const MSG_MONOQ         = 3072
const MAX_FIELDS_PER_MPI_COMM = 6
# Assume 128 byte coherence
# Assume Float64 is an "integral power of 2" bytes wide
const CACHE_COHERENCE_PAD_REAL = div(128, sizeof(Float64))

# Only supported configuration
const SEDOV_SYNC_POS_VEL_EARLY = true
const ALLOW_UNPACKED_PLANE = false

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
        col, row, plane, side = initMeshDecomp(comm)
        return new(num_iters, structured, nx, nr, balance, cost, col, row, plane, side, devicetype, floattype, comm)
    end
end

include("bc.jl")
include("types.jl")
include("domain.jl")
include("mpi.jl")
include("comm.jl")
include("utils.jl")

export printUsage, IndexT, Domain, LuleshProblem, getMyRank, getNumRanks, getWtime,
       lagrangeLeapFrog, comm_max, timeIncrement!, nodalMass, commRecv, MSG_COMM_SBN, verifyAndWriteFinalOutput
export commSend, commRecv, commSBN

function initMeshDecomp(comm)
    # Assume cube processor layout for now
    numRanks = getNumRanks(comm)
    myRank = getMyRank(comm)
    testProcs = floor(cbrt(numRanks+0.5))
    if (testProcs*testProcs*testProcs != numRanks)
        error("Num processors must be a cube of an integer (1, 8, 27, ...)")
    end
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
    col = myDom % dx
    row = div(myDom, dx) % dy
    plane = div(myDom, (dx*dy))
    side = testProcs

    return col, row, plane, side
end
export InitMeshDecomp

function printNormField(domain, fields)
    if getMyRank(domain.comm) == 0
        for field in fields
            # println("$field[", length(getfield(domain, field)), "]: ", norm(getfield(domain, field)))
        end
    end
end

function printField(domain, fields)
    if getMyRank(domain.comm) == 0
        for field in fields
            # println("$field[", length(getfield(domain, field)), "]: ", getfield(domain, field))
        end
    end
end

function printNormAllFields(domain, location="")
    if getMyRank(domain.comm) == 0
        # println("Location: ", location)
    end
    # fields = [:x, :xd, :fx, :nodalMass, :symmX, :dxx, :delv_xi]
    # fields = [:x, :commDataSend, :commDataRecv]
    # printNormField(domain, fields)
end


function printAllFields(domain, location="")
    if getMyRank(domain.comm) == 0
        # println("Location: ", location)
    end
    # fields = [:x, :xd, :fx, :nodalMass, :symmX, :nodelist]
    # fields = [:commDataSend, :commDataRecv]
    # printField(domain, fields)
end

export printField, printNormField, printAllFields, printNormAllFields
end # module
