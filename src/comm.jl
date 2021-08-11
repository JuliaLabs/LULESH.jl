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

   # post recieve buffers for all incoming messages
   maxPlaneComm = xferFields * domain.maxPlaneSize
   maxEdgeComm  = xferFields * domain.maxEdgeSize
   pmsg = 0 # plane comm msg
   emsg = 0 # edge comm msg
   cmsg = 0 # corner comm msg

   # MPI_Status status[26] ;

   # assume communication to 6 neighbors by default
   rowMin,rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

   packable = true
   # TODO: What is this checking for?
   # Probably that fieldData is all the same
   for i = 1:(xferFields - 2)
       if fieldData[i+1] - fieldData[i] !=
          fieldData[i+2] - fieldData[i+1]
           packable = false
           break
       end
   end

   # XXX:
   # for (Index_t i=0; i<26; ++i) {
   #    domain->sendRequest[i] = MPI_REQUEST_NULL ;
   # }

   myRank = MPI.Comm_rank(domain.comm)

   # post sends
   if planeMin || planeMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      # static MPI_Datatype msgTypePlane ;
      # static bool packPlane ;
      sendCount = dx * dy

      # if (msgTypePlane == 0) {
      #    /* Create an MPI_struct for field data */
      #    if (ALLOW_UNPACKED_PLANE && packable) {

      #       MPI_Type_vector(xferFields, sendCount,
      #                       (fieldData[1] - fieldData[0]),
      #                       baseType, &msgTypePlane) ;
      #       MPI_Type_commit(&msgTypePlane) ;
      #       packPlane = false ;
      #    }
      #    else {
      #       msgTypePlane = baseType ;
      #       packPlane = true ;
      #    }
      # }
      @assert !ALLOW_UNPACKED_PLANE
      packPlane = true

      if planeMin
         # contiguous memory
         if packPlane
            srcOffset = 0
            offset = pmsg * maxPlaneComm
            for field in fields
               copyto!(domain.commDataSend, offset, field, srcOffset, sendCount)
               offset += sendCount
            end
            idx = pmsg * maxPlaneComm + 1
            src = view(domain.commDataSend, idx:(idx+(xferFields * sendCount)))
         else
            error("unimplemented")
            # destAddr = fieldData[0] ;
         end

         otherRank = myRank - domain.tp^2
         req = MPI.Isend(src, otherRank, msgType + myRank, domain.comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end

      if planeMax && doSend
         # contiguous memory
         if packPlane
            srcOffset = dx*dy*(dz - 1)
            offset = pmsg * maxPlaneComm
            for field in fields
               copyto!(domain.commDataSend, offset, field, srcOffset, sendCount)
               offset += sendCount
            end
            idx = pmsg * maxPlaneComm + 1
            src = view(domain.commDataSend, idx:(idx+(xferFields * sendCount)))
         else
            error("unimplemented")
         end

         otherRank = myRank + domain.tp^2
         req = MPI.Isend(src, otherRank, msgType + myRank, domain.comm)
         domain.sendRequest[pmsg+1] = req
         pmsg += 1
      end


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
   doRecv = false
   xferFields = 6 ; # x, y, z, xd, yd, zd
   fields = (x, y, z, xd, yd, zd)
   maxPlaneComm = xferFields * domain->maxPlaneSize ;
   maxEdgeComm  = xferFields * domain->maxEdgeSize ;
   pmsg = 0 # plane comm msg
   emsg = 0 # edge comm msg
   cmsg = 0 # corner comm msg
   dx = domain.sizeX + 1
   dy = domain.sizeY + 1
   dz = domain.sizeZ + 1

   # assume communication to 6 neighbors by default
   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

   myRank = MPI.Comm_rank(domain.comm)

   if planeMin | planeMax
      # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
      opCount = dx * dy

      if planeMin && doRecv
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])
         offset = pmsg * maxPlaneComm
         for field in fields
            copyto!(field, 0, domain.commDataRecv, offset, opCount)
            offset += opCount
         end
         pmsg += 1
      end

      if planeMax
         # contiguous memory
         MPI.Wait!(domain.recvRequest[pmsg+1])

         dest_offset = dx*dy*(dz - 1)
         offset = pmsg * maxPlaneComm
         for field in fields
            copyto!(field, dest_offset, domain.commDataRecv, offset, opCount)
            offset += opCount
         end
         pmsg += 1
      end

   if (rowMin | rowMax) {
      /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
      Index_t opCount = dx * dz ;

      if (rowMin && doRecv) {
         /* contiguous memory */
         srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
         MPI_Wait(&domain->recvRequest[pmsg], &status) ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            for (Index_t i=0; i<dz; ++i) {
               Real_t *destAddr = &fieldData[fi][i*dx*dy] ;
               for (Index_t j=0; j<dx; ++j) {
                  destAddr[j] = srcAddr[i*dx + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (rowMax) {
         /* contiguous memory */
         Index_t offset = dx*(dy - 1) ;
         srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
         MPI_Wait(&domain->recvRequest[pmsg], &status) ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            for (Index_t i=0; i<dz; ++i) {
               Real_t *destAddr = &fieldData[fi][offset + i*dx*dy] ;
               for (Index_t j=0; j<dx; ++j) {
                  destAddr[j] = srcAddr[i*dx + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }
   if (colMin | colMax) {
      /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
      Index_t opCount = dy * dz ;

      if (colMin && doRecv) {
         /* contiguous memory */
         srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
         MPI_Wait(&domain->recvRequest[pmsg], &status) ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            for (Index_t i=0; i<dz; ++i) {
               Real_t *destAddr = &fieldData[fi][i*dx*dy] ;
               for (Index_t j=0; j<dy; ++j) {
                  destAddr[j*dx] = srcAddr[i*dy + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
      if (colMax) {
         /* contiguous memory */
         Index_t offset = dx - 1 ;
         srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
         MPI_Wait(&domain->recvRequest[pmsg], &status) ;
         for (Index_t fi=0 ; fi<xferFields; ++fi) {
            for (Index_t i=0; i<dz; ++i) {
               Real_t *destAddr = &fieldData[fi][offset + i*dx*dy] ;
               for (Index_t j=0; j<dy; ++j) {
                  destAddr[j*dx] = srcAddr[i*dy + j] ;
               }
            }
            srcAddr += opCount ;
         }
         ++pmsg ;
      }
   }

   if (rowMin && colMin && doRecv) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = fieldData[fi] ;
         for (Index_t i=0; i<dz; ++i) {
            destAddr[i*dx*dy] = srcAddr[i] ;
         }
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMin && planeMin && doRecv) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = fieldData[fi] ;
         for (Index_t i=0; i<dx; ++i) {
            destAddr[i] = srcAddr[i] ;
         }
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMin && planeMin && doRecv) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = fieldData[fi] ;
         for (Index_t i=0; i<dy; ++i) {
            destAddr[i*dx] = srcAddr[i] ;
         }
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMax && colMax) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*dy - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dz; ++i) {
            destAddr[i*dx*dy] = srcAddr[i] ;
         }
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMax && planeMax) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*(dy-1) + dx*dy*(dz-1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dx; ++i) {
            destAddr[i] = srcAddr[i] ;
         }
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMax && planeMax) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*dy*(dz-1) + dx - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dy; ++i) {
            destAddr[i*dx] = srcAddr[i] ;
         }
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMax && colMin) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*(dy-1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dz; ++i) {
            destAddr[i*dx*dy] = srcAddr[i] ;
         }
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMin && planeMax) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*dy*(dz-1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dx; ++i) {
            destAddr[i] = srcAddr[i] ;
         }
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMin && planeMax) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*dy*(dz-1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dy; ++i) {
            destAddr[i*dx] = srcAddr[i] ;
         }
         srcAddr += dy ;
      }
      ++emsg ;
   }

   if (rowMin && colMax && doRecv) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dz; ++i) {
            destAddr[i*dx*dy] = srcAddr[i] ;
         }
         srcAddr += dz ;
      }
      ++emsg ;
   }

   if (rowMax && planeMin && doRecv) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx*(dy - 1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dx; ++i) {
            destAddr[i] = srcAddr[i] ;
         }
         srcAddr += dx ;
      }
      ++emsg ;
   }

   if (colMax && planeMin && doRecv) {
      srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                       emsg * maxEdgeComm] ;
      Index_t offset = dx - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
      for (Index_t fi=0 ; fi<xferFields; ++fi) {
         Real_t *destAddr = &fieldData[fi][offset] ;
         for (Index_t i=0; i<dy; ++i) {
            destAddr[i*dx] = srcAddr[i] ;
         }
         srcAddr += dy ;
      }
      ++emsg ;
   }


   if (rowMin && colMin && planeMin && doRecv) {
      /* corner at domain logical coord (0, 0, 0) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][0] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMin && colMin && planeMax) {
      /* corner at domain logical coord (0, 0, 1) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMin && colMax && planeMin && doRecv) {
      /* corner at domain logical coord (1, 0, 0) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMin && colMax && planeMax) {
      /* corner at domain logical coord (1, 0, 1) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMax && colMin && planeMin && doRecv) {
      /* corner at domain logical coord (0, 1, 0) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*(dy - 1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMax && colMin && planeMax) {
      /* corner at domain logical coord (0, 1, 1) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMax && colMax && planeMin && doRecv) {
      /* corner at domain logical coord (1, 1, 0) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }
      ++cmsg ;
   }
   if (rowMax && colMax && planeMax) {
      /* corner at domain logical coord (1, 1, 1) */
      Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm +
                                      cmsg * CACHE_COHERENCE_PAD_REAL] ;
      Index_t idx = dx*dy*dz - 1 ;
      MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
      for (Index_t fi=0; fi<xferFields; ++fi) {
         fieldData[fi][idx] = comBuf[fi] ;
      }void CommSyncPosVel(Domain *domain) {

         if (domain->numProcs == 1) return ;

         int myRank ;
         bool doRecv = false ;
         Index_t xferFields = 6 ; /* x, y, z, xd, yd, zd */
         Real_t *fieldData[6] ;
         Index_t maxPlaneComm = xferFields * domain->maxPlaneSize ;
         Index_t maxEdgeComm  = xferFields * domain->maxEdgeSize ;
         Index_t pmsg = 0 ; /* plane comm msg */
         Index_t emsg = 0 ; /* edge comm msg */
         Index_t cmsg = 0 ; /* corner comm msg */
         Index_t dx = domain->sizeX + 1 ;
         Index_t dy = domain->sizeY + 1 ;
         Index_t dz = domain->sizeZ + 1 ;
         MPI_Status status ;
         Real_t *srcAddr ;
         bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
         /* assume communication to 6 neighbors by default */
         rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
         if (domain->rowLoc == 0) {
            rowMin = false ;
         }
         if (domain->rowLoc == (domain->tp-1)) {
            rowMax = false ;
         }
         if (domain->colLoc == 0) {
            colMin = false ;
         }
         if (domain->colLoc == (domain->tp-1)) {
            colMax = false ;
         }
         if (domain->planeLoc == 0) {
            planeMin = false ;
         }
         if (domain->planeLoc == (domain->tp-1)) {
            planeMax = false ;
         }

         fieldData[0] = domain->x ;
         fieldData[1] = domain->y ;
         fieldData[2] = domain->z ;
         fieldData[3] = domain->xd ;
         fieldData[4] = domain->yd ;
         fieldData[5] = domain->zd ;

         MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

         if (planeMin | planeMax) {
            /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
            Index_t opCount = dx * dy ;

            if (planeMin && doRecv) {
               /* contiguous memory */
               srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
               MPI_Wait(&domain->recvRequest[pmsg], &status) ;
               for (Index_t fi=0 ; fi<xferFields; ++fi) {
                  Real_t *destAddr = fieldData[fi] ;
                  for (Index_t i=0; i<opCount; ++i) {
                     destAddr[i] = srcAddr[i] ;
                  }
                  srcAddr += opCount ;
               }
               ++pmsg ;
            }
            if (planeMax) {
               /* contiguous memory */
               Index_t offset = dx*dy*(dz - 1) ;
               srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
               MPI_Wait(&domain->recvRequest[pmsg], &status) ;
               for (Index_t fi=0 ; fi<xferFields; ++fi) {
                  Real_t *destAddr = &fieldData[fi][offset] ;
                  for (Index_t i=0; i<opCount; ++i) {
                     destAddr[i] = srcAddr[i] ;
                  }
                  srcAddr += opCount ;
               }
               ++pmsg ;
            }
         }

         if (rowMin | rowMax) {
            /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
            Index_t opCount = dx * dz ;

            if (rowMin && doRecv) {
               /* contiguous memory */
               srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
               MPI_Wait(&domain->recvRequest[pmsg], &status) ;
               for (Index_t fi=0 ; fi<xferFields; ++fi) {
                  for (Index_t i=0; i<dz; ++i) {
                     Real_t *destAddr = &fieldData[fi][i*dx*dy] ;
                     for (Index_t j=0; j<dx; ++j) {
                        destAddr[j] = srcAddr[i*dx + j] ;
                     }
                  }
                  srcAddr += opCount ;
               }
               ++pmsg ;
            }
            if (rowMax) {
               /* contiguous memory */
               Index_t offset = dx*(dy - 1) ;
               srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
               MPI_Wait(&domain->recvRequest[pmsg], &status) ;
               for (Index_t fi=0 ; fi<xferFields; ++fi) {
                  for (Index_t i=0; i<dz; ++i) {
                     Real_t *destAddr = &fieldData[fi][offset + i*dx*dy] ;
                     for (Index_t j=0; j<dx; ++j) {
                        destAddr[j] = srcAddr[i*dx + j] ;
                     }
                  }
                  srcAddr += opCount ;
               }
               ++pmsg ;
            }
         }
         if (colMin | colMax) {
            /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
            Index_t opCount = dy * dz ;

            if (colMin && doRecv) {
               /* contiguous memory */
               srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
               MPI_Wait(&domain->recvRequest[pmsg], &status) ;
               for (Index_t fi=0 ; fi<xferFields; ++fi) {
                  for (Index_t i=0; i<dz; ++i) {
                     Real_t *destAddr = &fieldData[fi][i*dx*dy] ;
                     for (Index_t j=0; j<dy; ++j) {
                        destAddr[j*dx] = srcAddr[i*dy + j] ;
                     }
                  }
                  srcAddr += opCount ;
               }
               ++pmsg ;
            }
            if (colMax) {
               /* contiguous memory */
               Index_t offset = dx - 1 ;
               srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm] ;
               MPI_Wait(&domain->recvRequest[pmsg], &status) ;
               for (Index_t fi=0 ; fi<xferFields; ++fi) {
                  for (Index_t i=0; i<dz; ++i) {
                     Real_t *destAddr = &fieldData[fi][offset + i*dx*dy] ;
                     for (Index_t j=0; j<dy; ++j) {
                        destAddr[j*dx] = srcAddr[i*dy + j] ;
                     }
                  }
                  srcAddr += opCount ;
               }
               ++pmsg ;
            }
         }

         if (rowMin && colMin && doRecv) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = fieldData[fi] ;
               for (Index_t i=0; i<dz; ++i) {
                  destAddr[i*dx*dy] = srcAddr[i] ;
               }
               srcAddr += dz ;
            }
            ++emsg ;
         }

         if (rowMin && planeMin && doRecv) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = fieldData[fi] ;
               for (Index_t i=0; i<dx; ++i) {
                  destAddr[i] = srcAddr[i] ;
               }
               srcAddr += dx ;
            }
            ++emsg ;
         }

         if (colMin && planeMin && doRecv) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = fieldData[fi] ;
               for (Index_t i=0; i<dy; ++i) {
                  destAddr[i*dx] = srcAddr[i] ;
               }
               srcAddr += dy ;
            }
            ++emsg ;
         }

         if (rowMax && colMax) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*dy - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dz; ++i) {
                  destAddr[i*dx*dy] = srcAddr[i] ;
               }
               srcAddr += dz ;
            }
            ++emsg ;
         }

         if (rowMax && planeMax) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*(dy-1) + dx*dy*(dz-1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dx; ++i) {
                  destAddr[i] = srcAddr[i] ;
               }
               srcAddr += dx ;
            }
            ++emsg ;
         }

         if (colMax && planeMax) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*dy*(dz-1) + dx - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dy; ++i) {
                  destAddr[i*dx] = srcAddr[i] ;
               }
               srcAddr += dy ;
            }
            ++emsg ;
         }

         if (rowMax && colMin) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*(dy-1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dz; ++i) {
                  destAddr[i*dx*dy] = srcAddr[i] ;
               }
               srcAddr += dz ;
            }
            ++emsg ;
         }

         if (rowMin && planeMax) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*dy*(dz-1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dx; ++i) {
                  destAddr[i] = srcAddr[i] ;
               }
               srcAddr += dx ;
            }
            ++emsg ;
         }

         if (colMin && planeMax) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*dy*(dz-1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dy; ++i) {
                  destAddr[i*dx] = srcAddr[i] ;
               }
               srcAddr += dy ;
            }
            ++emsg ;
         }

         if (rowMin && colMax && doRecv) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dz; ++i) {
                  destAddr[i*dx*dy] = srcAddr[i] ;
               }
               srcAddr += dz ;
            }
            ++emsg ;
         }

         if (rowMax && planeMin && doRecv) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx*(dy - 1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dx; ++i) {
                  destAddr[i] = srcAddr[i] ;
               }
               srcAddr += dx ;
            }
            ++emsg ;
         }

         if (colMax && planeMin && doRecv) {
            srcAddr = &domain->commDataRecv[pmsg * maxPlaneComm +
                                             emsg * maxEdgeComm] ;
            Index_t offset = dx - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg], &status) ;
            for (Index_t fi=0 ; fi<xferFields; ++fi) {
               Real_t *destAddr = &fieldData[fi][offset] ;
               for (Index_t i=0; i<dy; ++i) {
                  destAddr[i*dx] = srcAddr[i] ;
               }
               srcAddr += dy ;
            }
            ++emsg ;
         }


         if (rowMin && colMin && planeMin && doRecv) {
            /* corner at domain logical coord (0, 0, 0) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][0] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMin && colMin && planeMax) {
            /* corner at domain logical coord (0, 0, 1) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx*dy*(dz - 1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMin && colMax && planeMin && doRecv) {
            /* corner at domain logical coord (1, 0, 0) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMin && colMax && planeMax) {
            /* corner at domain logical coord (1, 0, 1) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMax && colMin && planeMin && doRecv) {
            /* corner at domain logical coord (0, 1, 0) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx*(dy - 1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMax && colMin && planeMax) {
            /* corner at domain logical coord (0, 1, 1) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMax && colMax && planeMin && doRecv) {
            /* corner at domain logical coord (1, 1, 0) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx*dy - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
         if (rowMax && colMax && planeMax) {
            /* corner at domain logical coord (1, 1, 1) */
            Real_t *comBuf = &domain->commDataRecv[pmsg * maxPlaneComm +
                                                   emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL] ;
            Index_t idx = dx*dy*dz - 1 ;
            MPI_Wait(&domain->recvRequest[pmsg+emsg+cmsg], &status) ;
            for (Index_t fi=0; fi<xferFields; ++fi) {
               fieldData[fi][idx] = comBuf[fi] ;
            }
            ++cmsg ;
         }
      }
      ++cmsg ;
   }
}
end
