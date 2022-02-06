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
   if comm === nothing
       return
   end
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

   fill!(domain.recvRequest, MPI.Request())

   myRank = MPI.Comm_rank(comm)

   # post receives
   function irecv!(fromProc, offset, recvCount)
      idx = offset + 1
      data = MPI.Buffer(view(domain.commDataRecv, idx:(idx+recvCount-1)))
      return MPI.Irecv!(data, fromProc, msgType, domain.comm::MPI.Comm)
   end

   # receive data from neighboring domain faces
   if planeMin && doRecv
      # contiguous memory
      fromProc = myRank - domain.m_tp^2
      recvCount = dx * dy * xferFields
      req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
      domain.recvRequest[pmsg+1] = req
      pmsg += 1
   end

   if planeMax
      # contiguous memory
      fromProc = myRank + domain.m_tp^2
      recvCount = dx * dy * xferFields
      req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
      domain.recvRequest[pmsg+1] = req
      pmsg += 1
   end

   if rowMin && doRecv
      # semi-contiguous memory
      fromProc = myRank - domain.m_tp
      recvCount = dx * dz * xferFields
      req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
      domain.recvRequest[pmsg+1] = req
      pmsg += 1
   end

   if rowMax
      # semi-contiguous memory
      fromProc = myRank + domain.m_tp
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
         fromProc = myRank - domain.m_tp - 1
         recvCount = dz * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && planeMin && doRecv
         fromProc = myRank - domain.m_tp^2 - domain.m_tp
         recvCount = dx * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMin && planeMin && doRecv
         fromProc = myRank - domain.m_tp^2 - 1
         recvCount = dy * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && colMax
         fromProc = myRank + domain.m_tp + 1
         recvCount = dz * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && planeMax
         fromProc = myRank + domain.m_tp^2 + domain.m_tp
         recvCount = dx * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMax && planeMax
         fromProc = myRank + domain.m_tp^2 + 1
         recvCount = dy * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && colMin
         fromProc = myRank + domain.m_tp - 1
         recvCount = dz * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && planeMax
         fromProc = myRank + domain.m_tp^2 - domain.m_tp
         recvCount = dx * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMin && planeMax
         fromProc = myRank + domain.m_tp^2 - 1
         recvCount = dy * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && colMax && doRecv
         fromProc = myRank - domain.m_tp + 1
         recvCount = dz * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && planeMin && doRecv
         fromProc = myRank - domain.m_tp^2 + domain.m_tp
         recvCount = dx * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMax && planeMin && doRecv
         fromProc = myRank - domain.m_tp^2 + 1
         recvCount = dy * xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      # receive data from domains connected only by a corner
      if rowMin && colMin && planeMin && doRecv
         # corner at domain logical coord (0, 0, 0)
         fromProc = myRank - domain.m_tp^2 - domain.m_tp - 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMin && colMin && planeMax
         # corner at domain logical coord (0, 0, 1)
         fromProc = myRank + domain.m_tp^2 - domain.m_tp - 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMin && colMax && planeMin && doRecv
         # corner at domain logical coord (1, 0, 0)
         fromProc = myRank - domain.m_tp^2 - domain.m_tp + 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMin && colMax && planeMax
         # corner at domain logical coord (1, 0, 1)
         fromProc = myRank + domain.m_tp^2 - domain.m_tp + 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMin && planeMin && doRecv
         # corner at domain logical coord (0, 1, 0)
         fromProc = myRank - domain.m_tp^2 + domain.m_tp - 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMin && planeMax
         # corner at domain logical coord (0, 1, 1)
         fromProc = myRank + domain.m_tp^2 + domain.m_tp - 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMax && planeMin && doRecv
         # corner at domain logical coord (1, 1, 0)
         fromProc = myRank - domain.m_tp^2 + domain.m_tp + 1
         recvCount = xferFields
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         req = irecv!(fromProc, offset, recvCount)
         domain.recvRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMax && planeMax
         # corner at domain logical coord (1, 1, 1)
         fromProc = myRank + domain.m_tp^2 + domain.m_tp + 1
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

   comm = domain.comm
   if comm === nothing
      return
   end

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

   fill!(domain.sendRequest, MPI.Request())
   myRank = MPI.Comm_rank(comm)

   # post sends
   if planeMin | planeMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      sendCount = dx * dy

      if planeMin
         # contiguous memory
         srcOffset = 0
         offset = pmsg * maxPlaneComm
         for field in fields
            copyto_zero!(domain.commDataSend, offset, field, srcOffset, sendCount)
            offset += sendCount
         end
         idx = pmsg * maxPlaneComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount-1))))

         otherRank = myRank - domain.m_tp^2
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end

      if planeMax && doSend
         # contiguous memory
         srcOffset = dx*dy*(dz - 1)
         offset = pmsg * maxPlaneComm
         for field in fields
            copyto_zero!(domain.commDataSend, offset, field, srcOffset, sendCount)
            offset += sendCount
         end
         idx = pmsg * maxPlaneComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount -1))))

         otherRank = myRank + domain.m_tp^2
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end
   end

   if rowMin | rowMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      sendCount = dx * dz

      if rowMin
         # contiguous memory
         srcOffset = 0
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               copyto_zero!(domain.commDataSend, offset+i*dx, field, srcOffset+i*dx*dy, dx)
            end
            offset += sendCount
         end
         idx = pmsg * maxPlaneComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount-1))))
         otherRank = myRank - domain.m_tp
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end
      if rowMax && doSend
         # contiguous memory
         srcOffset = dx*(dy - 1)
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               copyto_zero!(domain.commDataSend, offset+i*dx , field, srcOffset+i*dx*dy, dx)
            end
            offset += sendCount
         end
         idx = pmsg * maxPlaneComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount -1))))
         otherRank = myRank + domain.m_tp
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end
   end

   if colMin | colMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      sendCount = dy * dz

      if colMin
         # contiguous memory
         srcOffset = 0
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dy-1)
                  domain.commDataSend[offset+i*dy+j + 1] = field[srcOffset+i*dx*dy+j*dx + 1]
               end
            end
            offset += sendCount
         end
         idx = pmsg * maxPlaneComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount -1))))
         otherRank = myRank - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end

      if colMax && doSend
         # contiguous memory
         srcOffset = dx - 1
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dy-1)
                  domain.commDataSend[offset+i*dy+j + 1] = field[srcOffset+j*dx + 1]
               end
            end
            offset += sendCount
         end
         idx = pmsg * maxPlaneComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount -1))))
         otherRank = myRank + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end
   end
   if !planeOnly
      if rowMin && colMin
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         for field in fields
            for i in 0:(dz-1)
               domain.commDataSend[offset+i + 1] = field[i*dx*dy + 1]
            end
            offset += dz
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dz -1))))
         otherRank = myRank - domain.m_tp - 1

         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && planeMin
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         for field in fields
            for i in 0:(dx-1)
               domain.commDataSend[offset+i + 1] = field[i + 1]
            end
            offset += dx
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dx -1))))
         otherRank = myRank - domain.m_tp^2 - domain.m_tp
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMin && planeMin
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         for field in fields
            for i in 0:(dy-1)
               domain.commDataSend[offset+i + 1] = field[i*dx + 1]
            end
            offset += dy
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dy -1))))
         otherRank = myRank - domain.m_tp^2 - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && colMax && doSend
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*dy - 1
         for field in fields
            for i in 0:(dz-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx*dy + 1]
            end
            offset += dz
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dz -1))))
         otherRank = myRank + domain.m_tp + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && planeMax && doSend
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*(dy-1) + dx*dy*(dz-1)
         for field in fields
            for i in 0:(dx-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i + 1]
            end
            offset += dx
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dx -1))))
         toRank = myRank + domain.m_tp^2 + domain.m_tp
         req = MPI.Isend(src, toRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMax && planeMax && doSend
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*dy*(dz-1) + dx - 1
         for field in fields
            for i in 0:(dy-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx + 1]
            end
            offset += dy
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dy -1))))
         otherRank = myRank + domain.m_tp^2 + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && colMin && doSend
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*(dy-1)
         for field in fields
            for i in 0:(dz-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx*dy + 1]
            end
            offset += dz
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dz -1))))
         otherRank = myRank + domain.m_tp - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && planeMax && doSend
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*dy*(dz-1)
         for field in fields
            for i in 0:(dx-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i + 1]
            end
            offset += dx
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dx -1))))
         otherRank = myRank + domain.m_tp^2 - domain.m_tp
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMin && planeMax && doSend
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*dy*(dz-1)
         for field in fields
            for i in 0:(dy-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx + 1]
            end
            offset += dy
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dy -1))))
         otherRank = myRank + domain.m_tp^2 - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && colMax
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx - 1
         for field in fields
            for i in 0:(dz-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx*dy + 1]
            end
            offset += dz
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dz -1))))
         otherRank = myRank - domain.m_tp + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMax && planeMin
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx*(dy - 1)
         for field in fields
            for i in 0:(dx-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i + 1]
            end
            offset += dx
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dx -1))))
         otherRank = myRank - domain.m_tp^2 + domain.m_tp
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if colMax && planeMin
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
         srcOffset = dx - 1
         for field in fields
            for i in 0:(dy-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx + 1]
            end
            offset += dy
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * dy -1))))
         otherRank = myRank - domain.m_tp^2 + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end

      if rowMin && colMin && planeMin
         # corner at domain logical coord (0, 0, 0) */
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm +  cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank - domain.m_tp^2 - domain.m_tp - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMin && colMin && planeMax && doSend
         # corner at domain logical coord (0, 0, 1)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[dx*dy*(dz-1) + 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank + domain.m_tp^2 - domain.m_tp - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMin && colMax && planeMin
         # corner at domain logical coord (1, 0, 0)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[(dx - 1) + 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank - domain.m_tp^2 - domain.m_tp + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMin && colMax && planeMax && doSend
         # corner at domain logical coord (1, 0, 1)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[dx*dy*(dz - 1) + (dx - 1) + 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank + domain.m_tp^2 - domain.m_tp + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMin && planeMin
         # corner at domain logical coord (0, 1, 0)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[dx*(dy - 1) + 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank - domain.m_tp^2 + domain.m_tp - 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMin && planeMax && doSend
         # corner at domain logical coord (0, 1, 1)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[dx*dy*(dz - 1) + dx*(dy - 1)+ 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank + domain.m_tp^2 + domain.m_tp - 1
         req = MPI.Isend(MPI.Buffer(src), otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMax && planeMin
         # corner at domain logical coord (1, 1, 0)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[dx*dy - 1 + 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank - domain.m_tp^2 + domain.m_tp + 1
         req = MPI.Isend(MPI.Buffer(src), otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end

      if rowMax && colMax && planeMax && doSend
         # corner at domain logical coord (1, 1, 1)
         offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
         for (fi, field) in enumerate(fields)
            # fi is one based already
            domain.commDataSend[offset+fi] = field[dx*dy*dz - 1 + 1]
         end
         idx = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL + 1
         src = MPI.Buffer(view(domain.commDataSend, idx:(idx+xferFields -1)))
         otherRank = myRank + domain.m_tp^2 + domain.m_tp + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
         domain.sendRequest[pmsg+emsg+cmsg+1] = req
         cmsg += 1
      end
   end

   for i in 1:(pmsg+emsg+cmsg)
      MPI.Wait!(domain.sendRequest[i])
   end
   # MPI.Waitall!(domain.sendRequest)
end

function commSBN(domain::Domain, fields)
   comm = domain.comm
   if comm === nothing
      return
   end

   xferFields = length(fields)

   maxPlaneComm = xferFields * domain.maxPlaneSize
   maxEdgeComm  = xferFields * domain.maxEdgeSize
   pmsg = 0 # plane comm msg
   emsg = 0 # edge comm msg
   cmsg = 0 # corner comm msg

   dx = domain.sizeX + 1
   dy = domain.sizeY + 1
   dz = domain.sizeZ + 1

   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

   myRank = MPI.Comm_rank(comm)

   if planeMin | planeMax
      opCount = dx * dy ;

      if planeMin
         srcaddr = pmsg * maxPlaneComm;
         status = MPI.Wait!(domain.recvRequest[pmsg + 1]);
         for field in fields
            for i in 1:opCount
               field[i] += domain.commDataRecv[srcaddr+i]
            end
            srcaddr += opCount
         end
         pmsg += 1
      end
      if planeMax
         srcaddr = pmsg * maxPlaneComm;
         status = MPI.Wait!(domain.recvRequest[pmsg + 1]);
         for field in fields
            for i in 1:opCount
               field[dx*dy*(dz-1) + i] += domain.commDataRecv[srcaddr+i]
            end
            srcaddr += opCount
         end
         pmsg += 1
      end
   end

   if rowMin | rowMax
      opCount = dx * dz

      if rowMin
         srcaddr = pmsg * maxPlaneComm;
         status = MPI.Wait!(domain.recvRequest[pmsg + 1]);
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dx-1)
                  field[i*dx*dy + j + 1] += domain.commDataRecv[srcaddr + i*dx + j + 1]
               end
            end
            srcaddr += opCount
         end
         pmsg += 1
      end

      if rowMax
         srcaddr = pmsg * maxPlaneComm;
         status = MPI.Wait!(domain.recvRequest[pmsg + 1]);
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dx-1)
                  field[dx*(dy-1) + i*dx*dy + j + 1] += domain.commDataRecv[srcaddr + i*dx + j + 1]
               end
            end
            srcaddr += opCount
         end
         pmsg += 1
      end
   end

   if colMin | colMax
      opCount = dy * dz ;

      if colMin
         srcaddr = pmsg * maxPlaneComm;
         status = MPI.Wait!(domain.recvRequest[pmsg + 1]);
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dy-1)
                  field[i*dx*dy + j*dx + 1] += domain.commDataRecv[srcaddr + i*dy + j + 1]
               end
            end
            srcaddr += opCount
         end
         pmsg += 1
      end

      if colMax
         srcaddr = pmsg * maxPlaneComm;
         status = MPI.Wait!(domain.recvRequest[pmsg + 1]);
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dy-1)
                  field[dx - 1 + i*dx*dy + j*dx + 1] += domain.commDataRecv[srcaddr + i*dy + j + 1]
               end
            end
            srcaddr += opCount
         end
         pmsg += 1
      end
   end

   if rowMin && colMin
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dz-1)
            # (domain.*dest)(i*dx*dy) += srcaddr[i] ;
            field[i*dx*dy + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dz
      end
      emsg += 1
   end

   if rowMin && planeMin
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dx-1)
            # (domain.*dest)(i) += srcaddr[i] ;
            field[i + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dx
      end
      emsg += 1
   end

   if colMin && planeMin
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dy-1)
            # (domain.*dest)(i*dx) += srcaddr[i] ;
            field[i*dx + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dy
      end
      emsg += 1
   end

   if rowMax && colMax
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dz-1)
            # (domain.*dest)(dx*dy - 1 + i*dx*dy) += srcaddr[i] ;
            field[dx*dy - 1 + i*dx*dy + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dz
      end
      emsg += 1
   end

   if rowMax && planeMax
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dx-1)
            # (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) += srcaddr[i] ;
            field[dx*(dy-1) + dx*dy*(dz-1) + i + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dx
      end
      emsg += 1
   end

   if colMax && planeMax
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dy-1)
            # (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) += srcaddr[i] ;
            field[dx*dy*(dz-1) + dx - 1 + i*dx + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dy
      end
      emsg += 1
   end

   if rowMax && colMin
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dz-1)
            # (domain.*dest)(dx*(dy-1) + i*dx*dy) += srcaddr[i] ;
            field[dx*(dy-1) + i*dx*dy + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dz
      end
      emsg += 1
   end

   if rowMin && planeMax
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dx-1)
            # (domain.*dest)(dx*dy*(dz-1) + i) += srcaddr[i] ;
            field[dx*dy*(dz-1) + i + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dx
      end
      emsg += 1
   end

   if colMin && planeMax
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dy-1)
            # (domain.*dest)(dx*dy*(dz-1) + i*dx) += srcaddr[i] ;
            field[dx*dy*(dz-1) + i*dx + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dy
      end
      emsg += 1
   end

   if rowMin && colMax
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dz-1)
            # (domain.*dest)(dx - 1 + i*dx*dy) += srcaddr[i] ;
            field[dx - 1 + i*dx*dy + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dz
      end
      emsg += 1
   end

   if rowMax && planeMin
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dx-1)
            # (domain.*dest)(dx*(dy - 1) + i) += srcaddr[i] ;
            field[dx*(dy-1) + i + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dx
      end
      emsg += 1
   end

   if colMax && planeMin
      srcaddr = pmsg * maxPlaneComm + emsg * maxEdgeComm
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg + 1])
      for field in fields
         for i in 0:(dy-1)
            # (domain.*dest)(dx - 1 + i*dx) += srcaddr[i] ;
            field[dx - 1 + i*dx + 1] += domain.commDataRecv[srcaddr + i + 1]
         end
         srcaddr += dy
      end
      emsg += 1
   end

   if rowMin && colMin && planeMin
      idx = 0
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMin && colMin && planeMax
      idx = dx*dy*(dz - 1)
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMin && colMax && planeMin
      idx = dx - 1
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMin && colMax && planeMax
      idx = dx*dy*(dz - 1) + (dx - 1)
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMax && colMin && planeMin
      idx = dx*(dy - 1)
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMax && colMin && planeMax
      idx = dx*dy*(dz - 1) + dx*(dy - 1)
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMax && colMax && planeMin
      idx = dx*dy - 1
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end

   if rowMax && colMax && planeMax
      idx = dx*dy*dz - 1
      comBuf = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      status = MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg + 1])
      for (fi, field) in enumerate(fields)
         field[idx+1] += domain.commDataRecv[comBuf + fi] # No 1 here as already taken care of by fi
      end
      cmsg += 1
   end
   printAllFields(domain, "$(@__FILE__):$(@__LINE__)")
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
   comm = domain.comm
   if comm === nothing
      return
   end
   doRecv = false
   xferFields = 6 ; # x, y, z, xd, yd, zd
   fields = (domain.x, domain.y, domain.z, domain.xd, domain.yd, domain.zd)
   maxPlaneComm = xferFields * domain.maxPlaneSize
   maxEdgeComm  = xferFields * domain.maxEdgeSize
   pmsg = 0 # plane comm msg
   emsg = 0 # edge comm msg
   cmsg = 0 # corner comm msg
   dx = domain.sizeX + 1
   dy = domain.sizeY + 1
   dz = domain.sizeZ + 1

   # assume communication to 6 neighbors by default
   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

   myRank = MPI.Comm_rank(comm)

   if planeMin | planeMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dx * dy

      if planeMin && doRecv
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for field in fields
            copyto_zero!(field, 0, domain.commDataRecv, offset, opCount)
            offset += opCount
         end
         pmsg += 1
      end

      if planeMax
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])

         destOffset = dx*dy*(dz - 1)
         offset = pmsg * maxPlaneComm
         for field in fields
            copyto_zero!(field, destOffset, domain.commDataRecv, offset, opCount)
            offset += opCount
         end
         pmsg += 1
      end

   end
   if rowMin | rowMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dx * dz

      if rowMin && doRecv
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               destOffset = i*dx*dy
               copyto_zero!(field, destOffset, domain.commDataRecv, offset + i*dx, dx)
            end
            offset += opCount
         end
         pmsg += 1
      end
      if rowMax
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               destOffset = dx*(dy - 1) + i*dx*dy
               copyto_zero!(field, destOffset, domain.commDataRecv, offset + i*dx, dx)
            end
            offset += opCount
         end
         pmsg += 1
      end
   end
   if colMin | colMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dy * dz

      if colMin && doRecv
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dy-1)
                  field[i*dx*dy+j*dx + 1] = domain.commDataRecv[offset + i*dy + j + 1]
               end
            end
            offset += opCount
         end
         pmsg += 1
      end
      if colMax
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         destOffset = dx-1
         for field in fields
            for i in 0:(dz-1)
               for j in 0:(dy-1)
                  field[destOffset + i*dx*dy + j*dx + 1] = domain.commDataRecv[offset + i*dy + j + 1]
               end
            end
            offset += opCount
         end
         pmsg += 1
      end
   end

   if rowMin && colMin && doRecv
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      for field in fields
         for i in 0:(dz-1)
            field[i*dx*dy + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dz
      end
      emsg += 1
   end

   if rowMin && planeMin && doRecv
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      for field in fields
         copyto_zero!(field, 0, domain.commDataRecv, offset, dx)
         offset += dx
      end
      emsg += 1
   end

   if colMin && planeMin && doRecv
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      for field in fields
         for i in 0:(dy-1)
            field[i*dx + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dy
      end
      emsg += 1
   end

   if rowMax && colMax
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*dy - 1
      for field in fields
         for i in 0:(dz-1)
            field[destOffset + i*dx*dy + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dz
      end
      emsg += 1
   end

   if rowMax && planeMax
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*(dy-1) + dx*dy*(dz-1)
      for field in fields
         copyto_zero!(field, destOffset, domain.commDataRecv, offset, dx)
         offset += dx
      end
      emsg += 1
   end

   if colMax && planeMax
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*dy*(dz-1) + dx - 1
      for field in fields
         for i in 0:(dy-1)
            field[destOffset + i*dx + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dy
      end
      emsg += 1
   end

   if rowMax && colMin
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*(dy-1)
      for field in fields
         for i in 0:(dz-1)
            field[destOffset + i*dx*dy + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dz
      end
      emsg += 1
   end

   if rowMin && planeMax
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*dy*(dz-1)
      for field in fields
         copyto_zero!(field, destOffset, domain.commDataRecv, offset, dx)
         offset += dx
      end
      emsg += 1
   end

   if colMin && planeMax
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*dy*(dz-1)
      for field in fields
         for i in 0:(dy-1)
            field[destOffset + i*dx + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dy
      end
      emsg += 1
   end

   if rowMin && colMax && doRecv
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx-1
      for field in fields
         for i in 0:(dz-1)
            field[destOffset + i*dx*dy + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dz
      end
      emsg += 1
   end

   if rowMax && planeMin && doRecv
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx*(dy - 1)
      for field in fields
         copyto_zero!(field, destOffset, domain.commDataRecv, offset, dx)
         offset += dx
      end
      emsg += 1
   end

   if colMax && planeMin && doRecv
      MPI.Wait!(domain.recvRequest[pmsg+emsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm
      destOffset = dx - 1
      for field in fields
         for i in 0:(dy-1)
            field[destOffset + i*dx + 1] = domain.commDataRecv[offset + i + 1]
         end
         offset += dy
      end
      emsg += 1
   end

   if rowMin && colMin && planeMin && doRecv
      # corner at domain logical coord (0, 0, 0)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMin && colMin && planeMax
      # corner at domain logical coord (0, 0, 1)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx*dy*(dz - 1)
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMin && colMax && planeMin && doRecv
      # corner at domain logical coord (1, 0, 0)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx - 1
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMin && colMax && planeMax
      # corner at domain logical coord (1, 0, 1)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx*dy*(dz - 1) + (dx - 1)
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMax && colMin && planeMin && doRecv
      # corner at domain logical coord (0, 1, 0)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx*(dy - 1)
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMax && colMin && planeMax
      # corner at domain logical coord (0, 1, 1)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx*dy*(dz - 1) + dx*(dy - 1)
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMax && colMax && planeMin && doRecv
      # corner at domain logical coord (1, 1, 0)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx*dy - 1
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   if rowMax && colMax && planeMax
      # corner at domain logical coord (1, 1, 1)
      MPI.Wait!(domain.recvRequest[pmsg+emsg+cmsg+1])
      offset = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg * CACHE_COHERENCE_PAD_REAL
      idx = dx*dy*dz - 1
      for (fi, field) in enumerate(fields)
         # fi is one based
         field[idx + 1] = domain.commDataRecv[offset + fi]
      end
      cmsg += 1
   end

   return nothing
end
