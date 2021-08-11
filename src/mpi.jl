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

function connect(domain::Domain)
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
    return rowMin, rowMax, colMin, colMax, planeMin, planeMax
end

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

    # assume communication to 6 neighbors by default
    rowMin, rowMax, colMin, colMax, planeMin, planeMax = connect(domain)

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

#     # post recieve buffers for all incoming messages
#     maxPlaneComm = xferFields * domain->maxPlaneSize ;
#     maxEdgeComm  = xferFields * domain->maxEdgeSize ;
#     pmsg = 0 # plane comm msg
#     emsg = 0 # edge comm msg
#     cmsg = 0 # corner comm msg

#     # MPI_Status status[26] ;
#     # Real_t *destAddr ;
#     # bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
#     # bool packable ;

#     # assume communication to 6 neighbors by default
#     rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
#     if domain.rowLoc == 0
#        rowMin = false
#     end
#     if domain.rowLoc == (domain.tp-1)
#        rowMax = false
#     end
#     if domain.colLoc == 0
#        colMin = false
#     }
#     if domain.colLoc == (domain.tp-1)
#        colMax = false
#     end
#     if domain.planeLoc == 0
#        planeMin = false
#     end
#     if domain.planeLoc == (domain.tp-1)
#        planeMax = false
#     end

#     packable = true
#     # TODO: What is this checking for?
#     # Probably that fieldData is all the same
#     for i = 1:(xferFields - 2)
#         if fieldData[i+1] - fieldData[i] !=
#            fieldData[i+2] - fieldData[i+1]
#             packable = false
#             break
#         end
#     end

#     # XXX:
#     # for (Index_t i=0; i<26; ++i) {
#     #    domain->sendRequest[i] = MPI_REQUEST_NULL ;
#     # }

#     myRank = MPI.Comm_rank(domain.comm)

#     # post sends
#     if planeMin || planeMax
#         # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
#         # static MPI_Datatype msgTypePlane ;
#         # static bool packPlane ;
#         sendCount = dx * dy

#         # if (msgTypePlane == 0) {
#         #    /* Create an MPI_struct for field data */
#         #    if (ALLOW_UNPACKED_PLANE && packable) {

#         #       MPI_Type_vector(xferFields, sendCount,
#         #                       (fieldData[1] - fieldData[0]),
#         #                       baseType, &msgTypePlane) ;
#         #       MPI_Type_commit(&msgTypePlane) ;
#         #       packPlane = false ;
#         #    }
#         #    else {
#         #       msgTypePlane = baseType ;
#         #       packPlane = true ;
#         #    }
#         # }
#         @assert !ALLOW_UNPACKED_PLANE
#         packPlane = true

#         function stage!(offset, sendcount)
#             for fi in 1:xferFields
#                 idx = offset + 1
#                 range = idx:(idx+sendcount)
#                 domain.commDataSend[range] = fieldData[fi] # ???
#                 offset += sendcount
#             end
#         end


#         if planeMin
#             # contiguous memory
#             offset = pmsg * maxPlaneComm
#             if packPlane
#                 stage!(offset, sendCount)
#             else
#                 error("unimplemented")
#                 # destAddr = fieldData[0] ;
#             end

#           MPI_Isend(destAddr,
#                     (packPlane ? xferFields*sendCount : 1),
#                      msgTypePlane, myRank - domain->tp*domain->tp,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg]) ;
#           ++pmsg ;
#        }
#        if (planeMax && doSend) {
#           /* contiguous memory */
#           Index_t offset = dx*dy*(dz - 1) ;
#           if (packPlane) {
#              destAddr = &domain->commDataSend[pmsg * maxPlaneComm] ;
#              for (Index_t fi=0 ; fi<xferFields; ++fi) {
#                 Real_t *srcAddr = &fieldData[fi][offset] ;
#                 memcpy(destAddr, srcAddr, sendCount*sizeof(Real_t)) ;
#                 destAddr += sendCount ;
#              }
#              destAddr -= xferFields*sendCount ;
#           }
#           else {
#              destAddr = &fieldData[0][offset] ;
#           }

#           MPI_Isend(destAddr,
#                     (packPlane ? xferFields*sendCount : 1),
#                      msgTypePlane, myRank + domain->tp*domain->tp,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg]) ;
#           ++pmsg ;
#        }
#     }
#     if (rowMin | rowMax) {
#        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
#        static MPI_Datatype msgTypeRow ;
#        static bool packRow ;
#        int sendCount = dx * dz ;

#        if (msgTypeRow == 0) {
#           /* Create an MPI_struct for field data */
#           if (ALLOW_UNPACKED_ROW && packable) {

#              static MPI_Datatype msgTypePencil ;

#              /* dz pencils per plane */
#              MPI_Type_vector(dz, dx, dx * dy, baseType, &msgTypePencil) ;
#              MPI_Type_commit(&msgTypePencil) ;

#              MPI_Type_vector(xferFields, 1, (fieldData[1] - fieldData[0]),
#                              msgTypePencil, &msgTypeRow) ;
#              MPI_Type_commit(&msgTypeRow) ;
#              packRow = false ;
#           }
#           else {
#              msgTypeRow = baseType ;
#              packRow = true ;
#           }
#        }

#        if (rowMin) {
#           /* contiguous memory */
#           if (packRow) {
#              destAddr = &domain->commDataSend[pmsg * maxPlaneComm] ;
#              for (Index_t fi=0; fi<xferFields; ++fi) {
#                 Real_t *srcAddr = fieldData[fi] ;
#                 for (Index_t i=0; i<dz; ++i) {
#                    memcpy(&destAddr[i*dx], &srcAddr[i*dx*dy],
#                           dx*sizeof(Real_t)) ;
#                 }
#                 destAddr += sendCount ;
#              }
#              destAddr -= xferFields*sendCount ;
#           }
#           else {
#              destAddr = fieldData[0] ;
#           }

#           MPI_Isend(destAddr,
#                     (packRow ? xferFields*sendCount : 1),
#                      msgTypeRow, myRank - domain->tp,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg]) ;
#           ++pmsg ;
#        }
#        if (rowMax && doSend) {
#           /* contiguous memory */
#           Index_t offset = dx*(dy - 1) ;
#           if (packRow) {
#              destAddr = &domain->commDataSend[pmsg * maxPlaneComm] ;
#              for (Index_t fi=0; fi<xferFields; ++fi) {
#                 Real_t *srcAddr = &fieldData[fi][offset] ;
#                 for (Index_t i=0; i<dz; ++i) {
#                    memcpy(&destAddr[i*dx], &srcAddr[i*dx*dy],
#                           dx*sizeof(Real_t)) ;
#                 }
#                 destAddr += sendCount ;
#              }
#              destAddr -= xferFields*sendCount ;
#           }
#           else {
#              destAddr = &fieldData[0][offset] ;
#           }

#           MPI_Isend(destAddr,
#                     (packRow ? xferFields*sendCount : 1),
#                      msgTypeRow, myRank + domain->tp,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg]) ;
#           ++pmsg ;
#        }
#     }
#     if (colMin | colMax) {
#        /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
#        static MPI_Datatype msgTypeCol ;
#        static bool packCol ;
#        int sendCount = dy * dz ;

#        if (msgTypeCol == 0) {
#           /* Create an MPI_struct for field data */
#           if (ALLOW_UNPACKED_COL && packable) {

#              static MPI_Datatype msgTypePoint ;
#              static MPI_Datatype msgTypePencil ;

#              /* dy points per pencil */
#              MPI_Type_vector(dy, 1, dx, baseType, &msgTypePoint) ;
#              MPI_Type_commit(&msgTypePoint) ;

#              /* dz pencils per plane */
#              MPI_Type_vector(dz, 1, dx*dy, msgTypePoint, &msgTypePencil) ;
#              MPI_Type_commit(&msgTypePencil) ;

#              MPI_Type_vector(xferFields, 1, (fieldData[1] - fieldData[0]),
#                              msgTypePencil, &msgTypeCol) ;
#              MPI_Type_commit(&msgTypeCol) ;
#              packCol = false ;
#           }
#           else {
#              msgTypeCol = baseType ;
#              packCol = true ;
#           }
#        }

#        if (colMin) {
#           /* contiguous memory */
#           if (packCol) {
#              destAddr = &domain->commDataSend[pmsg * maxPlaneComm] ;
#              for (Index_t fi=0; fi<xferFields; ++fi) {
#                 for (Index_t i=0; i<dz; ++i) {
#                    Real_t *srcAddr = &fieldData[fi][i*dx*dy] ;
#                    for (Index_t j=0; j<dy; ++j) {
#                       destAddr[i*dy + j] = srcAddr[j*dx] ;
#                    }
#                 }
#                 destAddr += sendCount ;
#              }
#              destAddr -= xferFields*sendCount ;
#           }
#           else {
#              destAddr = fieldData[0] ;
#           }

#           MPI_Isend(destAddr,
#                     (packCol ? xferFields*sendCount : 1),
#                      msgTypeCol, myRank - 1,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg]) ;
#           ++pmsg ;
#        }
#        if (colMax && doSend) {
#           /* contiguous memory */
#           Index_t offset = dx - 1 ;
#           if (packCol) {
#              destAddr = &domain->commDataSend[pmsg * maxPlaneComm] ;
#              for (Index_t fi=0; fi<xferFields; ++fi) {
#                 for (Index_t i=0; i<dz; ++i) {
#                    Real_t *srcAddr = &fieldData[fi][i*dx*dy + offset] ;
#                    for (Index_t j=0; j<dy; ++j) {
#                       destAddr[i*dy + j] = srcAddr[j*dx] ;
#                    }
#                 }
#                 destAddr += sendCount ;
#              }
#              destAddr -= xferFields*sendCount ;
#           }
#           else {
#              destAddr = &fieldData[0][offset] ;
#           }

#           MPI_Isend(destAddr,
#                     (packCol ? xferFields*sendCount : 1),
#                      msgTypeCol, myRank + 1,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg]) ;
#           ++pmsg ;
#        }
#     }

#     if (!planeOnly) {
#        if (rowMin && colMin) {
#           int toProc = myRank - domain->tp - 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = fieldData[fi] ;
#              for (Index_t i=0; i<dz; ++i) {
#                 destAddr[i] = srcAddr[i*dx*dy] ;
#              }
#              destAddr += dz ;
#           }
#           destAddr -= xferFields*dz ;
#           MPI_Isend(destAddr, xferFields*dz, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMin && planeMin) {
#           int toProc = myRank - domain->tp*domain->tp - domain->tp ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = fieldData[fi] ;
#              for (Index_t i=0; i<dx; ++i) {
#                 destAddr[i] = srcAddr[i] ;
#              }
#              destAddr += dx ;
#           }
#           destAddr -= xferFields*dx ;
#           MPI_Isend(destAddr, xferFields*dx, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (colMin && planeMin) {
#           int toProc = myRank - domain->tp*domain->tp - 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = fieldData[fi] ;
#              for (Index_t i=0; i<dy; ++i) {
#                 destAddr[i] = srcAddr[i*dx] ;
#              }
#              destAddr += dy ;
#           }
#           destAddr -= xferFields*dy ;
#           MPI_Isend(destAddr, xferFields*dy, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMax && colMax && doSend) {
#           int toProc = myRank + domain->tp + 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*dy - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dz; ++i) {
#                 destAddr[i] = srcAddr[i*dx*dy] ;
#              }
#              destAddr += dz ;
#           }
#           destAddr -= xferFields*dz ;
#           MPI_Isend(destAddr, xferFields*dz, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMax && planeMax && doSend) {
#           int toProc = myRank + domain->tp*domain->tp + domain->tp ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*(dy-1) + dx*dy*(dz-1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dx; ++i) {
#                destAddr[i] = srcAddr[i] ;
#              }
#              destAddr += dx ;
#           }
#           destAddr -= xferFields*dx ;
#           MPI_Isend(destAddr, xferFields*dx, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (colMax && planeMax && doSend) {
#           int toProc = myRank + domain->tp*domain->tp + 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*dy*(dz-1) + dx - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dy; ++i) {
#                 destAddr[i] = srcAddr[i*dx] ;
#              }
#              destAddr += dy ;
#           }
#           destAddr -= xferFields*dy ;
#           MPI_Isend(destAddr, xferFields*dy, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMax && colMin && doSend) {
#           int toProc = myRank + domain->tp - 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*(dy-1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dz; ++i) {
#                 destAddr[i] = srcAddr[i*dx*dy] ;
#              }
#              destAddr += dz ;
#           }
#           destAddr -= xferFields*dz ;
#           MPI_Isend(destAddr, xferFields*dz, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMin && planeMax && doSend) {
#           int toProc = myRank + domain->tp*domain->tp - domain->tp ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*dy*(dz-1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dx; ++i) {
#                 destAddr[i] = srcAddr[i] ;
#              }
#              destAddr += dx ;
#           }
#           destAddr -= xferFields*dx ;
#           MPI_Isend(destAddr, xferFields*dx, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (colMin && planeMax && doSend) {
#           int toProc = myRank + domain->tp*domain->tp - 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*dy*(dz-1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dy; ++i) {
#                 destAddr[i] = srcAddr[i*dx] ;
#              }
#              destAddr += dy ;
#           }
#           destAddr -= xferFields*dy ;
#           MPI_Isend(destAddr, xferFields*dy, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMin && colMax) {
#           int toProc = myRank - domain->tp + 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dz; ++i) {
#                 destAddr[i] = srcAddr[i*dx*dy] ;
#              }
#              destAddr += dz ;
#           }
#           destAddr -= xferFields*dz ;
#           MPI_Isend(destAddr, xferFields*dz, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (rowMax && planeMin) {
#           int toProc = myRank - domain->tp*domain->tp + domain->tp ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx*(dy - 1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dx; ++i) {
#                 destAddr[i] = srcAddr[i] ;
#              }
#              destAddr += dx ;
#           }
#           destAddr -= xferFields*dx ;
#           MPI_Isend(destAddr, xferFields*dx, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }

#        if (colMax && planeMin) {
#           int toProc = myRank - domain->tp*domain->tp + 1 ;
#           destAddr = &domain->commDataSend[pmsg * maxPlaneComm +
#                                            emsg * maxEdgeComm] ;
#           Index_t offset = dx - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              Real_t *srcAddr = &fieldData[fi][offset] ;
#              for (Index_t i=0; i<dy; ++i) {
#                 destAddr[i] = srcAddr[i*dx] ;
#              }
#              destAddr += dy ;
#           }
#           destAddr -= xferFields*dy ;
#           MPI_Isend(destAddr, xferFields*dy, baseType, toProc,
#                     (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg]) ;
#           ++emsg ;
#        }


#        if (rowMin && colMin && planeMin) {
#           /* corner at domain logical coord (0, 0, 0) */
#           int toProc = myRank - domain->tp*domain->tp - domain->tp - 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                        cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][0] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMin && colMin && planeMax && doSend) {
#           /* corner at domain logical coord (0, 0, 1) */
#           int toProc = myRank + domain->tp*domain->tp - domain->tp - 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx*dy*(dz - 1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMin && colMax && planeMin) {
#           /* corner at domain logical coord (1, 0, 0) */
#           int toProc = myRank - domain->tp*domain->tp - domain->tp + 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMin && colMax && planeMax && doSend) {
#           /* corner at domain logical coord (1, 0, 1) */
#           int toProc = myRank + domain->tp*domain->tp - domain->tp + 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMax && colMin && planeMin) {
#           /* corner at domain logical coord (0, 1, 0) */
#           int toProc = myRank - domain->tp*domain->tp + domain->tp - 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx*(dy - 1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMax && colMin && planeMax && doSend) {
#           /* corner at domain logical coord (0, 1, 1) */
#           int toProc = myRank + domain->tp*domain->tp + domain->tp - 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMax && colMax && planeMin) {
#           /* corner at domain logical coord (1, 1, 0) */
#           int toProc = myRank - domain->tp*domain->tp + domain->tp + 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx*dy - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#        if (rowMax && colMax && planeMax && doSend) {
#           /* corner at domain logical coord (1, 1, 1) */
#           int toProc = myRank + domain->tp*domain->tp + domain->tp + 1 ;
#           Real_t *comBuf = &domain->commDataSend[pmsg * maxPlaneComm +
#                                                  emsg * maxEdgeComm +
#                                           cmsg * CACHE_COHERENCE_PAD_REAL] ;
#           Index_t idx = dx*dy*dz - 1 ;
#           for (Index_t fi=0; fi<xferFields; ++fi) {
#              comBuf[fi] = fieldData[fi][idx] ;
#           }
#           MPI_Isend(comBuf,
#                     xferFields, baseType,
#                     toProc, (msgType + myRank) , MPI_COMM_WORLD,
#                     &domain->sendRequest[pmsg+emsg+cmsg]) ;
#           ++cmsg ;
#        }
#    }

#    MPI_Waitall(26, domain->sendRequest, status) ;
# }

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

    xferFields = 3 # delv_xi, delv_eta, delv_zeta
    fields = (domain.delv_xi, domain.delv_eta, domain.delv_zeta)
    maxPlaneComm = xferFields * domain.maxPlaneSize

    pmsg = 0 # plane comm msg
    dx = domain.sizeX
    dy = domain.sizeY
    dz = domain.sizeZ

    # assume communication to 6 neighbors by default
    rowMin, rowMax, colMin, colMax, planeMin, planeMax = connect(domain)

    numElem = domain.numElem

    # point into ghost data area
    fieldOffsets = [numElem, numElem, numElem]

    myRank = MPI.Comm_rank(domain.comm)

    if planeMin || planeMax
        # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
        opCount = dx * dy

        if planeMin
            # contiguous memory
            MPI.Wait!(domain.recvRequest[pmsg+1])

            offset = pmsg * maxPlaneComm
            for (fi, field) in enumerate(fields)
                fOffset = fieldOffsets[fi]
                copyto!(field, fOffset, domain.commDataRecv, offset, opCount)

                offset += opCount
                fieldOffsets[fi] += opCount
            end
            pmsg += 1
        end
        if planeMax
            # contiguous memory
            MPI.Wait!(domain.recvRequest[pmsg+1])
            offset = pmsg * maxPlaneComm
            for (fi, field) in enumerate(fields)
                fOffset = fieldOffsets[fi]
                copyto!(field, fOffset, domain.commDataRecv, offset, opCount)

                offset += opCount
                fieldOffsets[fi] += opCount
            end
            pmsg += 1
        end
    end

    error("unimplemented")
#     if (rowMin | rowMax) {
#       /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
#       Index_t opCount = dx * dz ;

#       if (rowMin) {
#          /* contiguous memory */
#          srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
#          MPI_Wait(&domain->recvRequest[pmsg], &status) ;
#          for (Index_t fi=0 ; fi<xferFields; ++fi) {
#             Real_t *destAddr = fieldData[fi] ;
#             for (Index_t i=0; i<opCount; ++i) {
#                destAddr[i] = srcAddr[i] ;
#             }
#             srcAddr += opCount ;
#             fieldData[fi] += opCount ;
#          }
#          ++pmsg ;
#       }
#       if (rowMax) {
#          /* contiguous memory */
#          srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
#          MPI_Wait(&domain->recvRequest[pmsg], &status) ;
#          for (Index_t fi=0 ; fi<xferFields; ++fi) {
#             Real_t *destAddr = fieldData[fi] ;
#             for (Index_t i=0; i<opCount; ++i) {
#                destAddr[i] = srcAddr[i] ;
#             }
#             srcAddr += opCount ;
#             fieldData[fi] += opCount ;
#          }
#          ++pmsg ;
#       }
#    }
#    if (colMin | colMax) {
#       /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
#       Index_t opCount = dy * dz ;

#       if (colMin) {
#          /* contiguous memory */
#          srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
#          MPI_Wait(&domain->recvRequest[pmsg], &status) ;
#          for (Index_t fi=0 ; fi<xferFields; ++fi) {
#             Real_t *destAddr = fieldData[fi] ;
#             for (Index_t i=0; i<opCount; ++i) {
#                destAddr[i] = srcAddr[i] ;
#             }
#             srcAddr += opCount ;
#             fieldData[fi] += opCount ;
#          }
#          ++pmsg ;
#       }
#       if (colMax) {
#          /* contiguous memory */
#          srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
#          MPI_Wait(&domain->recvRequest[pmsg], &status) ;
#          for (Index_t fi=0 ; fi<xferFields; ++fi) {
#             Real_t *destAddr = fieldData[fi] ;
#             for (Index_t i=0; i<opCount; ++i) {
#                destAddr[i] = srcAddr[i] ;
#             }
#             srcAddr += opCount ;
#          }
#          ++pmsg ;
end

function commSyncPosVel(domain::Domain)
    if domain.comm === nothing
        return
    end
    error("not implemented")
end
