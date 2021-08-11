   # assume communication to 6 neighbors by default
using MPI
function  get_neighbors(domain::Domain)
   rowMin = rowMax = colMin = colMax = planeMin = planeMax = true

   domain.m_rowLoc == 0 ? rowMin = false : rowMin = true
   domain.m_rowLoc == domain.m_tp - 1 ? rowMax = false : rowMax = true
   domain.m_colLoc == 0 ? colMin = false : colMin = true
   domain.m_colLoc == domain.m_tp - 1 ? colMax = false : colMax = true
   domain.m_planeLoc == 0 ? planeMin = false : planeMin = true
   domain.m_planeLoc == domain.m_tp - 1 ? planeMax = false : planeMax = true

   return rowMin, rowMax, colMin, colMax, planeMin, planeMax
end

# doRecv flag only works with regular block structure
function commRecv(
      domain::Domain, msgType::Int, xferFields::IndexT,
      dx::IndexT, dy::IndexT, dz::IndexT, doRecv::Bool, planeOnly::Bool, comm::MPI.Comm
   )
   maxPlaneComm = xferFields * domain.maxPlaneSize
   maxEdgeComm  = xferFields * domain.maxEdgeSize
   rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)
   pmsg = emsg = cmsg = 1
   # for r in domain.recvRequest
   #    r = MPI.REQUEST_NULL
   # end
   myRank = getMyRank(comm)

   # post receives

   # receive data from neighboring domain faces
   if (planeMin && doRecv)
      # contiguous memory
      fromRank = myRank - domain.m_tp*domain.m_tp
      recvCount = dx * dy * xferFields
      buf = @view domain.commDataRecv[pmsg * maxPlaneComm:pmsg * maxPlaneComm + recvCount]
      domain.recvRequest[pmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
      pmsg += 1
   end
   if (planeMax)
      # contiguous memory
       fromRank = myRank + domain.m_tp*domain.m_tp
       recvCount = dx * dy * xferFields
      buf = @view domain.commDataRecv[pmsg * maxPlaneComm:pmsg * maxPlaneComm + recvCount]
      domain.recvRequest[pmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
      pmsg += 1
   end
   if (rowMin && doRecv)
      # semi-contiguous memory
       fromRank = myRank - domain.m_tp
       recvCount = dx * dz * xferFields
      buf = @view domain.commDataRecv[pmsg * maxPlaneComm:pmsg * maxPlaneComm + recvCount]
      domain.recvRequest[pmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
      pmsg += 1
   end
   if (rowMax)
      # semi-contiguous memory
       fromRank = myRank + domain.m_tp
       recvCount = dx * dz * xferFields
      buf = @view domain.commDataRecv[pmsg * maxPlaneComm:pmsg * maxPlaneComm + recvCount]
      domain.recvRequest[pmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
      pmsg += 1
   end
   if (colMin && doRecv)
      # scattered memory
       fromRank = myRank - 1
       recvCount = dy * dz * xferFields
      buf = @view domain.commDataRecv[pmsg * maxPlaneComm:pmsg * maxPlaneComm + recvCount]
      domain.recvRequest[pmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
      pmsg += 1
   end
   if (colMax)
      # scattered memory
       fromRank = myRank + 1
       recvCount = dy * dz * xferFields
      buf = @view domain.commDataRecv[pmsg * maxPlaneComm:pmsg * maxPlaneComm + recvCount]
      domain.recvRequest[pmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
      pmsg += 1
   end

   if (!planeOnly)
      # receive data from domains connected only by an edge
      if (rowMin && colMin && doRecv)
          fromRank = myRank - domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dz *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMin && planeMin && doRecv)
          fromRank = myRank - domain.m_tp*domain.m_tp - domain.m_tp
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dx *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (colMin && planeMin && doRecv)
          fromRank = myRank - domain.m_tp*domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dy *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMax && colMax)
          fromRank = myRank + domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dz *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMax && planeMax)
          fromRank = myRank + domain.m_tp*domain.m_tp + domain.m_tp
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dx *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (colMax && planeMax)
          fromRank = myRank + domain.m_tp*domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dy *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMax && colMin)
          fromRank = myRank + domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dz *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMin && planeMax)
          fromRank = myRank + domain.m_tp*domain.m_tp - domain.m_tp
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dx *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (colMin && planeMax)
          fromRank = myRank + domain.m_tp*domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dy *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMin && colMax && doRecv)
          fromRank = myRank - domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dz *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (rowMax && planeMin && doRecv)
          fromRank = myRank - domain.m_tp*domain.m_tp + domain.m_tp
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dx *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      if (colMax && planeMin && doRecv)
          fromRank = myRank - domain.m_tp*domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm
      buf = @view domain.commDataRecv[start:start + dy *xferFields]
      domain.recvRequest[pmsg+emsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         emsg += 1
      end

      # receive data from domains connected only by a corner
      if (rowMin && colMin && planeMin && doRecv)
         # corner at domain logical coord (0, 0, 0)
          fromRank = myRank - domain.m_tp*domain.m_tp - domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMin && colMin && planeMax)
         # corner at domain logical coord (0, 0, 1)
          fromRank = myRank + domain.m_tp*domain.m_tp - domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMin && colMax && planeMin && doRecv)
         # corner at domain logical coord (1, 0, 0)
          fromRank = myRank - domain.m_tp*domain.m_tp - domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMin && colMax && planeMax)
         # corner at domain logical coord (1, 0, 1)
          fromRank = myRank + domain.m_tp*domain.m_tp - domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMax && colMin && planeMin && doRecv)
         # corner at domain logical coord (0, 1, 0)
          fromRank = myRank - domain.m_tp*domain.m_tp + domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMax && colMin && planeMax)
         # corner at domain logical coord (0, 1, 1)
          fromRank = myRank + domain.m_tp*domain.m_tp + domain.m_tp - 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMax && colMax && planeMin && doRecv)
         # corner at domain logical coord (1, 1, 0)
          fromRank = myRank - domain.m_tp*domain.m_tp + domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
      if (rowMax && colMax && planeMax)
         # corner at domain logical coord (1, 1, 1)
         fromRank = myRank + domain.m_tp*domain.m_tp + domain.m_tp + 1
        start = pmsg * maxPlaneComm + emsg * maxEdgeComm + cmsg
      buf = @view domain.commDataRecv[start:start + xferFields]
      domain.recvRequest[pmsg+emsg+cmsg] = MPI.Irecv!(buf, fromRank, msgType, comm)
         cmsg += 1
      end
   end
end

# function CommSend(
# 		domain::Domain, msgType::Int,
# 		dx::IndexT, dy::IndexT, dz::IndexT, doSend::Bool, planeOnly::Bool, comm::MPI.Comm
# 	)
# 	maxPlaneComm = xferFields * domain.maxPlaneSize()
# 	maxEdgeComm  = xferFields * domain.maxEdgeSize()
# 	rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)
# 	pmsg = emsg = cmsg = 1
# 	domain.sendRequest .= MPI.MPI_REQUEST_NULL
# 	myRank = getMyRank(comm)

# #    Index_t xferFields = sizeof...(fields);
# #    Domain_member fieldData[xferFields];
# #    {size_t i=0;
# #        ((fieldData[i++] = fields),...);
# #     }
