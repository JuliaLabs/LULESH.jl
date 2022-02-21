using LLVM
LLVM.clopts("-memdep-block-scan-limit=70000")
LLVM.clopts("-dse-memoryssa-walklimit=10000")
LLVM.clopts("-attributor-max-iterations=128")
LLVM.clopts("-capture-tracking-max-uses-to-explore=256")

using LULESH
using MPI
using Enzyme

Enzyme.API.printperf!(true)
Enzyme.API.printall!(true)
Enzyme.API.instname!(true)

Enzyme.API.inlineall!(true)
# Size of Domain is around 1024
Enzyme.API.maxtypeoffset!(1024)
# Enzyme.API.maxtypeoffset!(32)
isdefined(Enzyme.API, :strictAliasing!) && Enzyme.API.strictAliasing!(true)
isdefined(Enzyme.API, :typeWarning!) &&  Enzyme.API.typeWarning!(false)
Enzyme.API.looseTypeAnalysis!(true)

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
         src = MPI.Buffer(domain.commDataSend)
         req = MPI.Isend(src, 0, 0, comm)
	 
         MPI.Recv!(src, 0, 0, comm)
         
	MPI.Wait!(req)
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

function main(nx, structured, num_iters, mpi, enzyme)
    # TODO: change default nr to 11
    nr = 1
    balance = 1
    cost = 1
    floattype = Float64
    devicetype = Vector

        !MPI.Initialized() && MPI.Init()
        comm = MPI.COMM_WORLD

    prob = LuleshProblem(num_iters, structured, nx, nr, balance, cost, devicetype, floattype, comm)

    # Set up the mesh and decompose. Assumes regular cubes for now
    # TODO: modify this constructor to account for new fields
    # TODO: setup communication buffers

    domain = Domain(prob)
        shadowDomain = Domain(prob)
     
	commDataSend = Vector{Float64}(undef, domain.sizeX * domain.sizeY * domain.sizeZ)
        scommDataSend = Vector{Float64}(undef, domain.sizeX * domain.sizeY * domain.sizeZ)

   dx = domain.sizeX + 1
   dy = domain.sizeY + 1
   dz = domain.sizeZ + 1
	domx = Vector{Float64}(undef, 29791)
	sdomx = Vector{Float64}(undef, 29791)

	if enzyme
            Enzyme.autodiff(foo, Duplicated(domain, shadowDomain), Duplicated(domx, sdomx), dx, dy, dz)
        else
            foo(domain, domx, dx, dy, dz)
        end
        MPI.Finalize()
end

if !isinteractive()
    args = LULESH.parse_cmd()

    # assume cube subdomain geometry for now (nx)
    nx = args["N"]
    structured = args["s"]
    num_iters = args["num_iters"]
    mpi = args["mpi"]
    enzyme = args["enzyme"]
    main(nx, structured, num_iters, mpi, enzyme)
end
