function setupCommBuffers!(domain::Domain, edgeNodes::Int)
    @unpack_Domain domain
    # allocate a buffer large enough for nodal ghost data
    maxEdgeSize = max(domain.sizeX, max(domain.sizeY, domain.sizeZ))+1
    m_maxPlaneSize = maxEdgeSize*maxEdgeSize
    m_maxEdgeSize = maxEdgeSize

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
            (m_rowMax & m_colMax & m_planeMax))
            )

    domain.commDataSend = Vector{Float64}(undef, comBufSize)
    domain.commDataRecv = Vector{Float64}(undef, comBufSize)
    # prevent floating point exceptions
    fill!(domain.commDataSend, 0)
    fill!(domain.commDataRecv, 0)

    # Boundary nodesets
    if (m_colLoc == 0)
        resize!(domain.symmX, edgeNodes*edgeNodes)
    end
    if (m_rowLoc == 0)
        resize!(domain.symmY, edgeNodes*edgeNodes)
    end
    if (m_planeLoc == 0)
        resize!(domain.symmZ, edgeNodes*edgeNodes)
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

# Assume 128 byte coherence
# Assume Float64 is an "integral power of 2" bytes wide
const CACHE_COHERENCE_PAD_REAL = div(128, sizeof(Float64))

function commRecv(domain::Domain, msgType, xferFields, dx, dy, dz, doRecv, planeOnly)
    if domain.comm === nothing
        return
    end

    # post recieve buffers for all incoming messages
    maxPlaneComm = xferFields * domain.maxPlaneSize
    maxEdgeComm  = xferFields * domain.maxEdgeSize
    pmsg = 0 # plane comm msg
    emsg = 0 # edge comm msg
    cmsg = 0 # corner comm msg

    baseType = MPI.Datatype(Float64) # TODO support Float32

    # assume communication to 6 neighbors by default
    rowMin = rowMax = colMin = colMax = planeMin = planeMax = true
    if domain.rowLoc == 0
      rowMin = false
    end
    if domain.rowLoc == domain.tp-1
        rowMax = false
    end
    if domain.colLoc == 0
        colMin = false
    end
    if domain.colLoc == domain.tp-1
        colMax = false
    end
    if domain.planeLoc == 0
        planeMin = false
    end
    if domain.planeLoc == domain.tp-1
        planeMax = false
    end

    # XXX: empty recvRequests?
    # for (Index_t i=0; i<26; ++i) {
    #   domain->recvRequest[i] = MPI_REQUEST_NULL ;
    # }

    myRank = MPI.Comm_rank(domain.comm)

    # post receives
    function irecv!(fromProc, offset, recvCount)
        idx = offset + 1
        data = view(domain.commDataRecv, idx:(idx+recvCount))
        return MPI.Irecv!(data, fromProc, msgType + fromProc, MPI.MPI_COMM_WORLD)
    end

    # receive data from neighboring domain faces
    if planeMin && doRecv
        # contiguous memory
        fromProc = myRank - domain.tp^2
        recvCount = dx * dy * xferFields
        req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
        domain.recvRequest[pmsg+1] = req
        pmsg += 1
    end

    if planeMax
        # contiguous memory
        fromProc = myRank + domain.tp^2
        recvCount = dx * dy * xferFields
        req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
        domain.recvRequest[pmsg+1] = req
        pmsg += 1

    end

    if rowMin && doRecv
        # semi-contiguous memory
        fromProc = myRank - domain.tp
        recvCount = dx * dz * xferFields
        req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
        domain.recvRequest[pmsg+1] = req
        pmsg += 1
    end

    if rowMax
        # semi-contiguous memory
        fromProc = myRank + domain.tp
        recvCount = dx * dz * xferFields
        req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
        domain.recvRequest[pmsg+1] = req
        pmsg += 1
    end

    if colMin && doRecv
        # scattered memory
        fromProc = myRank - 1
        recvCount = dy * dz * xferFields
        req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
        domain.recvRequest[pmsg+1] = req
        pmsg += 1
    end

    if colMax
        # scattered memory
        fromProc = myRank + 1
        recvCount = dy * dz * xferFields
        req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
        domain.recvRequest[pmsg+1] = req
        pmsg += 1
    end

    if !planeOnly
        # receive data from domains connected only by an edge
        if rowMin && colMin && doRecv
            fromProc = myRank - domain.tp - 1
            recvCount = dz * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMin && colMax && doRecv
            fromProc = myRank - domain.tp^2 - domain.tp
            recvCount = dx * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if colMin && planeMin && doRecv
            fromProc = myRank - domain.tp^2 - 1
            recvCount = dy * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMax && colMax
            fromProc = myRank + domain.tp + 1
            recvCount = dz * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMax && planeMax
            fromProc = myRank + domain.tp^2 + domain.tp
            recvCount = dx * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if colMin && planeMax
            fromProc = myRank + domain.tp^2 + 1
            recvCount = dy * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMax && colMin
            fromProc = myRank + domain.tp + 1
            recvCount = dz * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMin && planeMax
            fromProc = myRank + domain.tp^2 - domain.tp
            recvCount = dx * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if colMin && planeMax
            fromProc = myRank + domain.tp^2 - 1
            recvCount = dy * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMax && colMax && doRecv
            fromProc = myRank - domain.tp + 1
            recvCount = dz * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if rowMax && planeMin && doRecv
            fromProc = myRank - domain.tp^2 + domain.tp
            recvCount = dx * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end

        if colMax && planeMin && doRecv
            fromProc = myRank - domain.tp^2 + 1
            recvCount = dy * xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+1] = req
            emsg += 1
        end


        # receive data from domains connected only by a corner
        if rowMin && colMin && planeMin && doRecv
            # corner at domain logical coord (0, 0, 0)
            fromProc = myRank - domain.tp^2 - domain.tp - 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMin && colMin && planeMax
            # corner at domain logical coord (0, 0, 1)
            fromProc = myRank + domain.tp^2 - domain.tp - 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMin && colMax && planeMin && doRecv
            # corner at domain logical coord (1, 0, 0)
            fromProc = myRank - domain.tp^2 - domain.tp + 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMin && colMax && planeMax
            # corner at domain logical coord (1, 0, 1)
            fromProc = myRank + domain.tp^2 - domain.tp + 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMax && colMin && planeMin && doRecv
            # corner at domain logical coord (0, 1, 0)
            fromProc = myRank - domain.tp^2 + domain.tp - 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMax && colMin && planeMax
            # corner at domain logical coord (0, 1, 1)
            fromProc = myRank + domain.tp^2 + domain.tp - 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMax && colMax && planeMin && doRecv
            # corner at domain logical coord (1, 1, 0)
            fromProc = myRank - domain.tp^2 + domain.tp + 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end

        if rowMax && colMax && planeMax
            # corner at domain logical coord (1, 1, 1)
            fromProc = myRank + domain.tp^2 + domain.tp + 1
            recvCount = xferFields
            offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
            req = irecv!(fromProc, offset, recvCount)
            domain.recvRequest[pmsg+emsg+cmsg+1] = req
            cmsg += 1
        end
    end
end

function commSend(domain::Domain, msgType, fields,
                  dx, dy, dz, doSend, planeOnly)

    if domain.comm === nothing
        return
    end

    xferFields = length(fields)
    error("not implemented")
end

function commSBN(domain::Domain, fields)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end

function commMonoQ(domain::Domain)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end

function commSyncPosVel(domain::Domain)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end
