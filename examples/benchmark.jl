using LLVM
LLVM.clopts("-memdep-block-scan-limit=70000")
LLVM.clopts("-dse-memoryssa-walklimit=10000")
LLVM.clopts("-attributor-max-iterations=128")
LLVM.clopts("-capture-tracking-max-uses-to-explore=256")

using MPI
using Enzyme

Enzyme.API.printperf!(true)
Enzyme.API.printall!(true)
Enzyme.API.instname!(true)

Enzyme.API.inlineall!(true)
Enzyme.API.maxtypeoffset!(1024)
isdefined(Enzyme.API, :strictAliasing!) && Enzyme.API.strictAliasing!(true)
isdefined(Enzyme.API, :typeWarning!) &&  Enzyme.API.typeWarning!(false)
Enzyme.API.looseTypeAnalysis!(true)

mutable struct Data
   commDataSend::Vector{Float64}
end

function free(buf)
  return nothing
end

function Isend(ar, buf, count, datatype, comm)
    req = MPI.Request()
    ccall((:MPI_Isend, MPI.libmpi), Cint,
          (MPI.MPIPtr, Cint, MPI.MPI_Datatype, Cint, Cint, MPI.MPI_Comm, Ptr{MPI.MPI_Request}),
                  buf.data, buf.count, datatype, 0, 0, comm, req)
    finalizer(free, req)
    return req
end

function fooSend(domain, fields,
                  dx, dy, dz, comm)
   	maxEdgeComm  = 6 * 32

	 offset = maxEdgeComm
         srcOffset = dx - 1
         for field in fields
            for i in 0:(dz-1)
               domain.commDataSend[offset+i + 1] = field[srcOffset+i*dx*dy + 1]
            end
            offset += dz
         end
	ar = domain.commDataSend
	datatype = MPI.Datatype(Float64)
         buf = MPI.Buffer(ar)
         req = Isend(ar, buf, buf.count, datatype, comm)
	 
    st = Ref{MPI.Status}(MPI.Status(0, 0, 0, 0, 0, 0))
    ccall((:MPI_Recv, MPI.libmpi), Cint,
                  (MPI.MPIPtr, Cint, MPI.MPI_Datatype, Cint, Cint, MPI.MPI_Comm, Ptr{MPI.Status}),
                   ar, buf.count, datatype, 0, 0, comm, st)
    
    stat_ref = Ref{MPI.Status}(MPI.Status(0, 0, 0, 0, 0, 0))
    ccall((:MPI_Wait, MPI.libmpi), Cint,
                  (Ptr{MPI.MPI_Request}, Ptr{MPI.Status}),
                  req, stat_ref)
    return nothing 
end
function foo(domain, domx, dx, dy, dz)

   # assume communication to 6 neighbors by default
   comm = MPI.COMM_WORLD
        
      fields = (domx, domx, domx, domx, domx, domx)
      fooSend(domain, fields,
                 dx, dy, dz,
		 comm)

    return nothing
end

function main(enzyme)
        comm = MPI.COMM_WORLD
     
	domain = Data(Vector{Float64}(undef, 30*30*30))
        shadowDomain = Data(Vector{Float64}(undef, 30*30*30))

   dx = 30 + 1
   dy = 30 + 1
   dz = 30 + 1
	domx = Vector{Float64}(undef, 29791)
	sdomx = Vector{Float64}(undef, 29791)

	if enzyme
            Enzyme.autodiff(foo, Duplicated(domain, shadowDomain), Duplicated(domx, sdomx), dx, dy, dz)
        else
            foo(domain, domx, dx, dy, dz)
        end
end

if !isinteractive()
    !MPI.Initialized() && MPI.Init()
    main(false)
    @show "ran primal"
    flush(stdout)
    main(true)
    MPI.Finalize()
end
