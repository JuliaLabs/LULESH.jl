using LULESH
using MPI
using Printf
using Enzyme

function main(nx, structured, num_iters, mpi, cuda)
    # TODO: change default nr to 11
    nr = 1
    balance = 1
    cost = 1
    floattype = Float64

    if cuda
        # devicetype = CUDA.CuArray
        error("CUDA not yet supported")
    else
        devicetype = Vector
    end

    if mpi
        !MPI.Initialized() && MPI.Init()
        comm = MPI.COMM_WORLD
    else
        comm = nothing
    end

    prob = LuleshProblem(num_iters, structured, nx, nr, balance, cost, devicetype, floattype, comm)

    if comm !== nothing && cuda
        local_comm = MPI.Comm_split_type(comm, MPI.MPI_COMM_TYPE_SHARED, MPI.Comm_rank(comm))
        @assert MPI.Comm_size(local_comm) <= length(CUDA.devices())
        CUDA.device!(MPI.Comm_rank(local_comm))
    end

    # Set up the mesh and decompose. Assumes regular cubes for now
    # TODO: modify this constructor to account for new fields
    # TODO: setup communication buffers

    domain = Domain(prob)
    shadowDomain = Domain(prob)

    getnodalMass = nodalMass(domain)

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
        MPI.Barrier(prob.comm)
    end

    if getMyRank(prob.comm) == 0
        if (structured)
            @printf("Running until t=%f, Problem size=%dx%dx%d\n", domain.stoptime, nx, nx, nx)
        else
            @printf("Running until t=%f, Problem size=%d \n", domain.stoptime, domain.numElem)
            @warn "Unstructured setup not supported"
        end
    end

    if cuda
        CUDA.Profile.start()
    end

    # timestep to solution
    start = getWtime(prob.comm)
    while domain.time < domain.stoptime
        # this has been moved after computation of volume forces to hide launch latencies
        timeIncrement!(domain)
        lagrangeLeapFrog(domain)
        # checkErrors(domain, its, myRank)
        if getMyRank(prob.comm) == 0
            @printf("cycle = %d, time = %e, dt=%e\n", domain.cycle, domain.time, domain.deltatime)
        end
        if domain.cycle >= num_iters
            break
        end
    end

    # make sure GPU finished its work
    if cuda
        CUDA.synchronize()
    end

    # Use reduced max elapsed time
    elapsed_time = getWtime(prob.comm) - start
    elapsed_timeG = comm_max(elapsed_time, prob.comm)

    if cuda
        CUDA.Profile.stop()
    end

    #   if (myRank == 0)
    #     VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, its, nx, numRanks, structured);

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
    main(nx, structured, num_iters, mpi, false)
end
