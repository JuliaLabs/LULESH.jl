using CUDA
using LULESH
using MPI


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
# nx = parse(IndexT, ARGS[2])
nx = 45
# TODO: change default nr to 11
nr = 11
balance = 1
cost = 1
devicetype = CuVector
floattype = Float64
prob = LuleshProblem(num_iters, structured, nx, nr, balance, cost, CuVector, Float64, MPI.COMM_WORLD)

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
locDom = NewDomain(prob)
# locDom = NewDomain(numRanks, col, row, plane, nx, side, structured, nr, balance, cost);

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

#   /* timestep to solution */
#   int its=0;

#   if (myRank == 0) {
#     if (structured)
#       printf("Running until t=%f, Problem size=%dx%dx%d\n",locDom->stoptime,nx,nx,nx);
#     else
#       printf("Running until t=%f, Problem size=%d \n",locDom->stoptime,locDom->numElem);
#   }

#   cudaProfilerStart();

# #if USE_MPI
#    double start = MPI_Wtime();
# #else
#    timeval start;
#    gettimeofday(&start, NULL) ;
# #endif

#   while(locDom->time_h < locDom->stoptime)
#   {
#     // this has been moved after computation of volume forces to hide launch latencies
#     //TimeIncrement(locDom) ;

#     LagrangeLeapFrog(locDom) ;

#     checkErrors(locDom,its,myRank);

#     #if LULESH_SHOW_PROGRESS
#      if (myRank == 0)
# 	 printf("cycle = %d, time = %e, dt=%e\n", its+1, double(locDom->time_h), double(locDom->deltatime_h) ) ;
#     #endif
#     its++;
#     if (its == num_iters) break;
#   }

#   // make sure GPU finished its work
#   cudaDeviceSynchronize();

# // Use reduced max elapsed time
#    double elapsed_time;
# #if USE_MPI
#    elapsed_time = MPI_Wtime() - start;
# #else
#    timeval end;
#    gettimeofday(&end, NULL) ;
#    elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
# #endif

#    double elapsed_timeG;
# #if USE_MPI
#    MPI_Reduce(&elapsed_time, &elapsed_timeG, 1, MPI_DOUBLE,
#               MPI_MAX, 0, MPI_COMM_WORLD);
# #else
#    elapsed_timeG = elapsed_time;
# #endif

#   cudaProfilerStop();

#   if (myRank == 0)
#     VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, its, nx, numRanks, structured);

# #ifdef SAMI
#   DumpDomain(locDom) ;
# #endif
#   cudaDeviceReset();

# #if USE_MPI
#    MPI_Finalize() ;
# #endif

#   return 0 ;
# }