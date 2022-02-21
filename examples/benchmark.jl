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
Enzyme.API.maxtypeoffset!(32)
isdefined(Enzyme.API, :strictAliasing!) && Enzyme.API.strictAliasing!(true)
isdefined(Enzyme.API, :typeWarning!) &&  Enzyme.API.typeWarning!(false)
Enzyme.API.looseTypeAnalysis!(true)

function main(nx, structured, num_iters, mpi, enzyme)
    # TODO: change default nr to 11
    nr = 1
    balance = 1
    cost = 1
    floattype = Float64
    devicetype = Vector

        !MPI.Initialized() && MPI.Init()
        comm = MPI.COMM_WORLD


      comm = MPI.COMM_WORLD
      myRank = MPI.Comm_rank(comm)
    domain = [1.0, 2.0, 3.0]
        shadowDomain = [4.0, 5.0, 6.0]
        if enzyme
            Enzyme.autodiff(lagrangeLeapFrog, Duplicated(domain, shadowDomain), myRank)
        else
            lagrangeLeapFrog(domain, myRank)
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
