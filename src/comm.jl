   # assume communication to 6 neighbors by default
using MPI
using Enzyme
function  get_neighbors(domain::Domain)
   rowMin = !(domain.m_rowLoc == 0)
   rowMax = !(domain.m_rowLoc == domain.m_tp - 1)
   colMin = !(domain.m_colLoc == 0)
   colMax = !(domain.m_colLoc == domain.m_tp - 1)
   planeMin = !(domain.m_planeLoc == 0)
   planeMax = !(domain.m_planeLoc == domain.m_tp - 1)

   return rowMin, rowMax, colMin, colMax, planeMin, planeMax
end

copyto_zero!(dest, doffs, src, soffs, nelems) = copyto!(dest, doffs+1, src, soffs+1, nelems)


function commRecv(domain::Domain, msgType, xferFields, dx, dy, dz, doRecv, planeOnly)
   comm = domain.comm
   comm = comm::MPI.Comm

   # post receive buffers for all incoming messages
   maxPlaneComm = xferFields * domain.maxPlaneSize
   maxEdgeComm  = xferFields * domain.maxEdgeSize
   pmsg = 0 # plane comm msg
   emsg = 0 # edge comm msg
   cmsg = 0 # corner comm msg

   baseType = MPI.Datatype(Float64) # TODO support Float32

   # assume communication to 6 neighbors by default
   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

   myRank = MPI.Comm_rank(comm)

   # post receives
   if planeMax
      # contiguous memory
      fromProc = myRank + domain.m_tp^2
      data = MPI.Buffer(view(domain.commDataRecv, 1:2))
      req = MPI.Irecv!(data, fromProc, msgType, MPI.COMM_WORLD)
      domain.recvRequest[1] = req
   end

end

function commSend(domain::Domain, msgType, fields,
                  dx, dy, dz, doSend, planeOnly)

   comm = domain.comm

   xferFields = length(fields)

   # post recieve buffers for all incoming messages
   maxPlaneComm = xferFields * domain.maxPlaneSize
   maxEdgeComm  = xferFields * domain.maxEdgeSize
   pmsg = 0 # plane comm msg
   emsg = 0 # edge comm msg
   cmsg = 0 # corner comm msg


   # MPI_Status status[26] ;

   # assume communication to 6 neighbors by default
   rowMin,rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

end


function commMonoQ(domain::Domain)
   comm = domain.comm
   if comm === nothing
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
   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

   numElem = domain.numElem

   # point into ghost data area
   fieldOffsets = [numElem, numElem, numElem]

   myRank = MPI.Comm_rank(comm)

   if planeMin | planeMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
      opCount = dx * dy

      if planeMin
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])

         offset = pmsg * maxPlaneComm
         for (fi, field) in enumerate(fields)
            fOffset = fieldOffsets[fi]
            copyto_zero!(field, fOffset, domain.commDataRecv, offset, opCount)

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
            copyto_zero!(field, fOffset, domain.commDataRecv, offset, opCount)

            offset += opCount
            fieldOffsets[fi] += opCount
         end
         pmsg += 1
      end
   end

   if rowMin | rowMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dx * dz

      if rowMin
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for (fi, field) in enumerate(fields)
            fOffset = fieldOffsets[fi]
            copyto_zero!(field, fOffset, domain.commDataRecv, offset, opCount)

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
            copyto_zero!(field, fOffset, domain.commDataRecv, offset, opCount)
            offset += opCount
            fieldOffsets[fi] += opCount
         end
         pmsg += 1
      end
   end
   if colMin | colMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dy * dz

      if colMin
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for (fi, field) in enumerate(fields)
            fOffset = fieldOffsets[fi]
            copyto_zero!(field, fOffset, domain.commDataRecv, offset, opCount)

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
            copyto_zero!(field, fOffset, domain.commDataRecv, offset, opCount)

            offset += opCount
            fieldOffsets[fi] += opCount
         end
         pmsg += 1
      end
   end
   printAllFields(domain, "$(@__FILE__):$(@__LINE__)")
   return nothing
end

function commSyncPosVel(domain::Domain)
   


   return nothing
end
