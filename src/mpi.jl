SetupCommBuffers(edgeNodes::Int, dom::AbstractDomain) = nothing
BuildMesh(nx::Int, edgeNodes::Int, edgeElems::Int, domNodes::Int, padded_domElems::Int, x_h::Vector{FT}, y_h::Vector{FT}, z_h::Vector{FT}, nodelist_h::Vector{Int}, dom::AbstractDomain) where {FT} = nothing

getMyRank(comm::MPI.Comm) = MPI.Comm_rank(comm)
getMyRank(::Nothing) = 0

getNumRanks(comm::MPI.Comm) = MPI.Comm_size(comm)
getNumRanks(::Nothing) = 1

getWtime(::MPI.Comm) = MPI.Wtime()
getWtime(::Nothing) = time()

comm_max(data::Float64, comm::MPI.Comm) = MPI.Allreduce(data, MPI.MAX, comm)
comm_max(data::Float64, ::Nothing) = data
