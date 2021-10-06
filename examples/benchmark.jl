using LULESH
using MPI
# using Enzyme

# Enzyme.API.printperf!(true)
# Enzyme.API.printall!(true)
# Enzyme.API.instname!(true)

# Enzyme.API.inlineall!(true)
# Enzyme.API.maxtypeoffset!(32)
# ccall((:EnzymeSetCLInteger, Enzyme.API.libEnzyme), Cvoid, (Ptr{Cvoid}, Int64), cglobal((:MaxTypeOffset, Enzyme.API.libEnzyme)), 32)

function main(nx, structured, num_iters, mpi)
    # TODO: change default nr to 11
    nr = 1
    balance = 1
    cost = 1
    floattype = Float64
    devicetype = Vector

    if mpi
        !MPI.Initialized() && MPI.Init()
        comm = MPI.COMM_WORLD
    else
        comm = nothing
    end

    prob = LuleshProblem(num_iters, structured, nx, nr, balance, cost, devicetype, floattype, comm)

    # Set up the mesh and decompose. Assumes regular cubes for now
    # TODO: modify this constructor to account for new fields
    # TODO: setup communication buffers

    domain = Domain(prob)
    shadowDomain = Domain(prob)

    # getnodalMass = nodalMass(domain)

    # Initial domain boundary communication
    # commRecv(domain, MSG_COMM_SBN, 1,
    #             domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
    #             true, false, prob.comm)
    #    CommSend<&Domain::nodalMass>(*domain, MSG_COMM_SBN,
    #             domain->sizeX() + 1, domain->sizeY() + 1, domain->sizeZ() +  1,
    #             true, false)
    #    CommSBN<&Domain::nodalMass>(*domain)

    # End initialization
    if mpi
        MPI.Barrier(prob.comm::MPI.Comm)
    end

    if getMyRank(prob.comm) == 0
        if (structured)
            @info "Running" until=domain.stoptime domain=(nx,nx,nx)
        else
            @info "Running" until=domain.stoptime domain=domain.numElem
            @warn "Unstructured setup not supported"
        end
    end

    # timestep to solution
    start = getWtime(prob.comm)
    while domain.time < domain.stoptime
        # this has been moved after computation of volume forces to hide launch latencies
        timeIncrement!(domain)
        lagrangeLeapFrog(domain)

        # checkErrors(domain, its, myRank)
        if getMyRank(prob.comm) == 0
            @info "Completed" cycle=domain.cycle time=domain.time dt=domain.deltatime
        end
        if domain.cycle >= num_iters
            break
        end
    end

    # Use reduced max elapsed time
    elapsed_time = getWtime(prob.comm) - start
    elapsed_timeG = comm_max(elapsed_time, prob.comm)

    if getMyRank(prob.comm) == 0
        verifyAndWriteFinalOutput(elapsed_timeG, domain, nx, getNumRanks(prob.comm))
    end

    if mpi
        MPI.Finalize()
    end
end

if !isinteractive()
    args = LULESH.parse_cmd()

    # assume cube subdomain geometry for now (nx)
    nx = args["N"]
    structured = args["s"]
    num_iters = args["num_iters"]
    mpi = args["mpi"]
    main(nx, structured, num_iters, mpi)
end
