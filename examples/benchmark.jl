using LULESH
using MPI
using Enzyme
using LLVM
using Printf

LLVM.clopts("-memdep-block-scan-limit=70000")
LLVM.clopts("-dse-memoryssa-walklimit=10000")
LLVM.clopts("-attributor-max-iterations=128")
LLVM.clopts("-capture-tracking-max-uses-to-explore=256")

# Enzyme.API.printperf!(true)
# Enzyme.API.printall!(true)
# Enzyme.API.instname!(true)

Enzyme.API.inlineall!(true)
# Size of Domain is around 1024
Enzyme.API.maxtypeoffset!(32)
isdefined(Enzyme.API, :strictAliasing!) && Enzyme.API.strictAliasing!(true)
isdefined(Enzyme.API, :typeWarning!) &&  Enzyme.API.typeWarning!(false)
Enzyme.API.looseTypeAnalysis!(true)

function main(nx, structured, num_iters, mpi, enzyme, verification)
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

    if getMyRank(comm) == 0
        @info "Constructing LuleshProblem" num_iters structured nx nr balance cost ranks=getNumRanks(comm)
    end

    prob = LuleshProblem(num_iters, structured, nx, nr, balance, cost, devicetype, floattype, comm)

    # Set up the mesh and decompose. Assumes regular cubes for now
    # TODO: modify this constructor to account for new fields
    # TODO: setup communication buffers

    domain = Domain(prob)
    if enzyme
        shadowDomain = Domain(prob)
        if verification
            @time Enzyme.autodiff(doubleFrog, Duplicated(domain, shadowDomain))
        else
            @time Enzyme.autodiff(lagrangeLeapFrog, Duplicated(domain, shadowDomain))
        end
        domain = Domain(prob)
    elseif verification
        # not enzyme and verification means we need perturbed domain for finite differences
        perturbedDomain = Domain(prob)
    end

    # getnodalMass = nodalMass(domain)

    # Initial domain boundary communication
    commRecv(domain, MSG_COMM_SBN, 1,
                domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
                true, false)
    fields = (domain.nodalMass,)
    commSend(domain, MSG_COMM_SBN, fields,
             domain.sizeX + 1, domain.sizeY + 1, domain.sizeZ + 1,
             true, false)
    commSBN(domain, fields)

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

    eps = 1e-6
    if verification
        num_iters = 1
        if enzyme
            println("seeding gradients")
            for f in fieldnames(Domain)
                prop = getproperty(shadowDomain, f)
                if(typeof(prop) == Float64)
                    setproperty!(shadowDomain,f,0.0)
                    @printf("set %s to %f\n",f,sum(getproperty(shadowDomain,f)))
                elseif(typeof(prop) == Vector{Float64})
                    fill!(prop, 0.0)
                    @printf("filled %s with %f\n",f,sum(getproperty(shadowDomain,f)))
                end
            end
            fill!(shadowDomain.e, 1.0)
        else
            # finite difference initialization
            println("seeding FD")
            domain.e[1] += eps
            perturbedDomain.e[1] -= eps
        end
    end

    # timestep to solution
    start = getWtime(prob.comm)
    while domain.time < domain.stoptime
        # this has been moved after computation of volume forces to hide launch latencies
        timeIncrement!(domain)
        if enzyme
            if verification
                Enzyme.autodiff(doubleFrog, Duplicated(domain, shadowDomain))
            else
                Enzyme.autodiff(lagrangeLeapFrog, Duplicated(domain, shadowDomain))
            end
        else
            if verification
                doubleFrog(domain)
                timeIncrement!(perturbedDomain)
                doubleFrog(perturbedDomain)
            else
                lagrangeLeapFrog(domain)
            end
        end

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

    if verification
        if enzyme
            println("harvesting gradients")
            @printf("gradient checksum: %f\n", shadowDomain.e[1])
        else
            @printf("FD checksum: %f\n",sum(domain.e - perturbedDomain.e) / (2*eps))
        end
    end

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
    enzyme = args["enzyme"]
    verification = args["verification"]
    main(nx, structured, num_iters, mpi, enzyme, verification)
end
