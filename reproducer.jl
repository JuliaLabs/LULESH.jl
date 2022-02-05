using MPI
using Enzyme
MPI.Init()

# function Isend(buf::Buffer, dest::Integer, tag::Integer, comm::Comm)
#     req = Request()
#     # int MPI_Isend(const void* buf, int count, MPI_Datatype datatype, int dest,
#     #               int tag, MPI_Comm comm, MPI_Request *request)
#     @mpichk ccall((:MPI_Isend, MPI.libmpi), Cint,
#           (MPIPtr, Cint, MPI_Datatype, Cint, Cint, MPI_Comm, Ptr{MPI_Request}),
#                   buf.data, buf.count, buf.datatype, dest, tag, comm, req)
#     req.buffer = buf
#     finalizer(free, req)
#     return req
# end

import MPI: libmpi, MPIPtr, MPI_Datatype, MPI_Request,
            MPI_Comm, @mpichk, Buffer, Request, free, Comm

function billysIrecv!(buf::Buffer)
    req = Request()
    # int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source,
    #               int tag, MPI_Comm comm, MPI_Request *request)
    ccall((:MPI_Irecv, libmpi), Cint,
                  (MPIPtr, Cint, MPI_Datatype, Cint, Cint, MPI_Comm, Ptr{MPI_Request}),
                  buf.data, 2, buf.datatype, 0, 0, MPI.COMM_WORLD, req)
    req.buffer = buf
    finalizer(req) do r
    end
    return nothing
end

function mycalcForceForNodes(data, myRank)
    data = MPI.Buffer(view(data, 1:2))
    billysIrecv!(data) #domain.comm::MPI.Comm)
   return nothing
end

Enzyme.API.printall!(true)
Enzyme.API.typeWarning!(false)
Enzyme.API.maxtypeoffset!(32)
Enzyme.API.inlineall!(true)
myRank = MPI.Comm_rank(MPI.COMM_WORLD)
x = [1.0, 2.0]
dx = [3.0, 5.0]
Enzyme.autodiff(mycalcForceForNodes, Duplicated(x, dx), myRank)

# using MPI
# using LULESH
# using Enzyme

# MPI.Init()

# function  get_neighbors(domain::Domain)
#    rowMin = !(domain.m_rowLoc == 0)
#    rowMax = !(domain.m_rowLoc == domain.m_tp - 1)
#    colMin = !(domain.m_colLoc == 0)
#    colMax = !(domain.m_colLoc == domain.m_tp - 1)
#    planeMin = !(domain.m_planeLoc == 0)
#    planeMax = !(domain.m_planeLoc == domain.m_tp - 1)

#    return rowMin, rowMax, colMin, colMax, planeMin, planeMax
# end

# copyto_zero!(dest, doffs, src, soffs, nelems) = copyto!(dest, doffs+1, src, soffs+1, nelems)

# function commRecv(domain::Domain, msgType, xferFields, dx, dy, dz, doRecv, planeOnly)
#    comm = domain.comm
#    if comm === nothing
#        return
#    end
#    comm = comm::MPI.Comm

#    # post receive buffers for all incoming messages
#    maxPlaneComm = xferFields * domain.maxPlaneSize
#    maxEdgeComm  = xferFields * domain.maxEdgeSize
#    pmsg = 0 # plane comm msg
#    emsg = 0 # edge comm msg
#    cmsg = 0 # corner comm msg

#    baseType = MPI.Datatype(Float64) # TODO support Float32

#    # assume communication to 6 neighbors by default
#    rowMin, rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

#    fill!(domain.recvRequest, MPI.Request())

#    myRank = MPI.Comm_rank(comm)

#    # post receives
#    function irecv!(fromProc, offset, recvCount)
#       idx = offset + 1
#       data = MPI.Buffer(view(domain.commDataRecv, idx:(idx+recvCount-1)))
#       return MPI.Irecv!(data, fromProc, msgType, domain.comm::MPI.Comm)
#    end

#    # receive data from neighboring domain faces
#    if planeMin && doRecv
#       # contiguous memory
#       fromProc = myRank - domain.m_tp^2
#       recvCount = dx * dy * xferFields
#       req = irecv!(fromProc, pmsg*maxPlaneComm, recvCount)
#       domain.recvRequest[pmsg+1] = req
#       pmsg += 1
#    end
# end

# function commSend(domain::Domain, msgType, fields,
#                   dx, dy, dz, doSend, planeOnly)

#    comm = domain.comm
#    if comm === nothing
#       return
#    end

#    xferFields = length(fields)

#    # post recieve buffers for all incoming messages
#    maxPlaneComm = xferFields * domain.maxPlaneSize
#    maxEdgeComm  = xferFields * domain.maxEdgeSize
#    pmsg = 0 # plane comm msg
#    emsg = 0 # edge comm msg
#    cmsg = 0 # corner comm msg


#    # MPI_Status status[26] ;

#    # assume communication to 6 neighbors by default
#    rowMin,rowMax, colMin, colMax, planeMin, planeMax = get_neighbors(domain)

#    fill!(domain.sendRequest, MPI.Request())
#    myRank = MPI.Comm_rank(comm)

#    # post sends
#    if planeMin | planeMax
#       # ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
#       sendCount = dx * dy

#       if planeMax && doSend
#          # contiguous memory
#          srcOffset = dx*dy*(dz - 1)
#          offset = pmsg * maxPlaneComm
#          for field in fields
#             copyto_zero!(domain.commDataSend, offset, field, srcOffset, sendCount)
#             offset += sendCount
#          end
#          idx = pmsg * maxPlaneComm + 1
#          src = MPI.Buffer(view(domain.commDataSend, idx:(idx+(xferFields * sendCount -1))))

#          otherRank = myRank + domain.m_tp^2
#          req = MPI.Isend(src, otherRank, msgType, comm)
#          domain.sendRequest[pmsg+1] = req
#          pmsg += 1
#       end
#    end

#    for i in 1:(pmsg+emsg+cmsg)
#     MPI.Wait!(domain.sendRequest[i])
#    end
# end

# function calcForceForNodes(domain::Domain)
#     commRecv(domain, MSG_COMM_SBN, 3,
#              domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
#              true, false)

#     fields = (domain.fx, domain.fy, domain.fz)
#     commSend(domain, MSG_COMM_SBN, fields,
#              domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
#              true, false)
#     # commSBN(domain, fields)
# end

# prop = LuleshProblem(1, true, 45, 1, 1, 1, Vector, Float64, MPI.COMM_WORLD)
# domain = Domain(prop)
# shadowDomain = Domain(prop)
# # LULESH.calcForceForNodes(domain)
# Enzyme.autodiff(LULESH.calcForceForNodes, Duplicated(domain, shadowDomain))
