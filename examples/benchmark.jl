using CUDA
using LULESH
using MPI
using Printf


if length(ARGS) < 2
    printUsage()
    error("Wrong input arguments")
end

if !occursin("-u", ARGS[1]) && !occursin("-s", ARGS[1])
    printUsage()
    error("Wrong input arguments")
end

num_iters = 10 # AD_LENGTH of simulation usually set to "-1" but for debugging set to 10!
if length(ARGS) == 4
    num_iters = parse(Int, ARGS[4])
end

#   bool structured = ( strcmp(argv[1],"-s") == 0 );
structured = occursin("-s", ARGS[1])
structured = true
@show structured
# assume cube subdomain geometry for now (nx)
nx = parse(IndexT, ARGS[2])
# nx = 45
# TODO: change default nr to 11
nr = 1
balance = 1
cost = 1
devicetype = Vector
floattype = Float64
# To disable MPI, pass 'nothing' as the last argument
prob = LuleshProblem(num_iters, structured, nx, nr, balance, cost, Vector, Float64, nothing)

#   Int_t numRanks ;
#   Int_t myRank ;

# #if USE_MPI
#   Domain_member fieldData ;

#   MPI_Init(&argc, &argv) ;
#   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
#   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
# #else
#   numRanks = 1;
#   myRank = 0;
# #endif

#   cuda_init(myRank);


#   Domain *locDom ;

# Set up the mesh and decompose. Assumes regular cubes for now


# TODO: modify this constructor to account for new fields
# TODO: setup communication buffers

domain = Domain(prob)

# #if USE_MPI
#    // copy to the host for mpi transfer
#    locDom->h_nodalMass = locDom->nodalMass;

#    fieldData = &Domain::get_nodalMass;

#    // Initial domain boundary communication
#    CommRecv(*locDom, MSG_COMM_SBN, 1,
#             locDom->sizeX + 1, locDom->sizeY + 1, locDom->sizeZ + 1,
#             true, false) ;
#    CommSend(*locDom, MSG_COMM_SBN, 1, &fieldData,
#             locDom->sizeX + 1, locDom->sizeY + 1, locDom->sizeZ + 1,
#             true, false) ;
#    CommSBN(*locDom, 1, &fieldData) ;

#    // copy back to the device
#    locDom->nodalMass = locDom->h_nodalMass;

#    // End initialization
#    MPI_Barrier(MPI_COMM_WORLD);
# #endif

#   cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#   cudaDeviceSetLimit(cudaLimitMallocHeapSize,1024*1024*1024);

@show domain.time_h
@show domain.stoptime
if getMyRank(prob.comm) == 0
    if (structured)
        @printf("Running until t=%f, Problem size=%dx%dx%d\n", domain.stoptime, nx, nx, nx)
    else
        @printf("Running until t=%f, Problem size=%d \n", domain.stoptime, domain.numElem)
        @warn "Unstructured setup not supported"
    end
end
# timestep to solution
#   cudaProfilerStart();
start = getWtime(prob.comm)
while domain.time_h < domain.stoptime
    # this has been moved after computation of volume forces to hide launch latencies
    timeIncrement!(domain)
    lagrangeLeapFrog(domain)
    # checkErrors(domain, its, myRank)
    if getMyRank(prob.comm) == 0
        @printf("cycle = %d, time = %e, dt=%e\n", domain.cycle, domain.time_h, domain.deltatime_h)
    end
    if domain.cycle == num_iters
        break
    end
end

#   // make sure GPU finished its work
#   cudaDeviceSynchronize();

# Use reduced max elapsed time
elapsed_time = getWtime(prob.comm) - start
elapsed_timeG = comm_max(elapsed_time, prob.comm)

#   cudaProfilerStop();

#   if (myRank == 0)
#     VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, its, nx, numRanks, structured);

# #ifdef SAMI
#   DumpDomain(locDom) ;
# #endif
#   cudaDeviceReset();

MPI.Finalize()