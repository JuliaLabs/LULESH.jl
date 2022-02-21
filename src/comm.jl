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
   function irecv!(fromProc, offset, recvCount)
      idx = offset + 1
      data = MPI.Buffer(view(domain.commDataRecv, 1:2))
      return MPI.Irecv!(data, fromProc, msgType, domain.comm::MPI.Comm)
   end

   if !planeOnly

      if rowMax && colMax
         fromProc = myRank + domain.m_tp + 1
         recvCount = dz * xferFields
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

   myRank = MPI.Comm_rank(comm)
      
      if rowMin && colMin
         src = MPI.Buffer(view(domain.commDataSend, 1:2))
         otherRank = myRank - domain.m_tp - 1

         req = MPI.Isend(src, otherRank, msgType, comm)
      	 MPI.Wait!(req)
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
         src = MPI.Buffer(view(domain.commDataSend, 1:2))
         otherRank = myRank - domain.m_tp + 1
         req = MPI.Isend(src, otherRank, msgType, comm)
      	 MPI.Wait!(req)
         domain.sendRequest[pmsg+emsg+1] = req
         emsg += 1
      end
end

function commSBN(domain::Domain, fields)
   comm = domain.comm

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
