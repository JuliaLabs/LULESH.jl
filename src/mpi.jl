function setupCommBuffers!(domain::Domain, edgeNodes::Int)
    function cache_align_real(n::Int64)
        n = Int32(n)
        (n + CACHE_COHERENCE_PAD_REAL%Int32 - one(Int32)) & ~(CACHE_COHERENCE_PAD_REAL%Int32 - one(Int32))
    end
    @unpack_Domain domain
    # allocate a buffer large enough for nodal ghost data
    maxEdgeSize = max(domain.sizeX, max(domain.sizeY, domain.sizeZ))+1

    m_maxPlaneSize = cache_align_real(maxEdgeSize*maxEdgeSize)
    m_maxEdgeSize = cache_align_real(maxEdgeSize)
    maxPlaneSize = m_maxPlaneSize
    maxEdgeSize = m_maxEdgeSize

    # assume communication to 6 neighbors by default
    m_rowMin = (m_rowLoc == 0)        ? 0 : 1
    m_rowMax = (m_rowLoc == m_tp-1)     ? 0 : 1
    m_colMin = (m_colLoc == 0)        ? 0 : 1
    m_colMax = (m_colLoc == m_tp-1)     ? 0 : 1
    m_planeMin = (m_planeLoc == 0)    ? 0 : 1
    m_planeMax = (m_planeLoc == m_tp-1) ? 0 : 1

    # account for face communication
    comBufSize =(
        (m_rowMin + m_rowMax + m_colMin + m_colMax + m_planeMin + m_planeMax) *
        m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM
    )

    # account for edge communication
    comBufSize += (
        ((m_rowMin & m_colMin) + (m_rowMin & m_planeMin) + (m_colMin & m_planeMin) +
        (m_rowMax & m_colMax) + (m_rowMax & m_planeMax) + (m_colMax & m_planeMax) +
        (m_rowMax & m_colMin) + (m_rowMin & m_planeMax) + (m_colMin & m_planeMax) +
        (m_rowMin & m_colMax) + (m_rowMax & m_planeMin) + (m_colMax & m_planeMin)) *
        m_maxEdgeSize * MAX_FIELDS_PER_MPI_COMM
    )

    # account for corner communication
    # factor of 16 is so each buffer has its own cache line
    comBufSize += (((m_rowMin & m_colMin & m_planeMin) +
            (m_rowMin & m_colMin & m_planeMax) +
            (m_rowMin & m_colMax & m_planeMin) +
            (m_rowMin & m_colMax & m_planeMax) +
            (m_rowMax & m_colMin & m_planeMin) +
            (m_rowMax & m_colMin & m_planeMax) +
            (m_rowMax & m_colMax & m_planeMin) +
            (m_rowMax & m_colMax & m_planeMax)) * CACHE_COHERENCE_PAD_REAL
            )

    commDataSend = Vector{Float64}(undef, comBufSize)
    commDataRecv = Vector{Float64}(undef, comBufSize)
    # prevent floating point exceptions
    fill!(commDataSend, 0)
    fill!(commDataRecv, 0)

    # Boundary nodesets
    if (m_colLoc == 0)
        resize!(symmX, edgeNodes*edgeNodes)
    end
    if (m_rowLoc == 0)
        resize!(symmY, edgeNodes*edgeNodes)
    end
    if (m_planeLoc == 0)
        resize!(symmZ, edgeNodes*edgeNodes)
    end
    @pack_Domain! domain
end

getMyRank(comm::MPI.Comm) = MPI.Comm_rank(comm)
getMyRank(::Nothing) = 0

getNumRanks(comm::MPI.Comm) = MPI.Comm_size(comm)
getNumRanks(::Nothing) = 1

getWtime(::MPI.Comm) = MPI.Wtime()
getWtime(::Nothing) = time()

comm_max(data::Float64, comm::MPI.Comm) = MPI.Allreduce(data, MPI.MAX, comm)
comm_max(data::Float64, ::Nothing) = data

comm_min(data::Float64, comm::MPI.Comm) = MPI.Allreduce(data, MPI.MIN, comm)
comm_min(data::Float64, ::Nothing) = data

comm_barrier(comm::MPI.Comm) = MPI.Barrier(comm)
comm_barrier(::Nothing) = nothing
