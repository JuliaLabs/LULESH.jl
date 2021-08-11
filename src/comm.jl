   # assume communication to 6 neighbors by default
using MPI
function  get_neighbors(domain::Domain)
   rowMin = !(domain.m_rowLoc == 0)
   rowMax = !(domain.m_rowLoc == domain.m_tp - 1)
   colMin = !(domain.m_colLoc == 0)
   colMax = !(domain.m_colLoc == domain.m_tp - 1)
   planeMin = !(domain.m_planeLoc == 0)
   planeMax = !(domain.m_planeLoc == domain.m_tp - 1)

   return rowMin, rowMax, colMin, colMax, planeMin, planeMax
end

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
   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

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

   # TODO:
   # Check copyto! offset semantics

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

   if rowMin || rowMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dx * dz

      if rowMin
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

      if rowMax
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
   if colMin || colMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dy * dz

      if colMin
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
      if colMax
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
   return nothing
end

function commSyncPosVel(domain::Domain)
   if domain.comm === nothing
      return
   end
   error("not implemented")
end
