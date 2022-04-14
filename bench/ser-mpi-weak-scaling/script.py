#!/usr/bin/python3.8
import os

def printfun(rank, blocklist,itercount):
  for s in blocklist:
    for mode in ["","--enzyme"]:
      os. chdir("/home/ubuntu/enzyme-sc22/LULESH.jl/")
      os.system("./mpiexecjl  -bind-to socket --project -np {}  /home/ubuntu/julia-1.7.2/bin/julia --project examples/benchmark.jl -s {}  --mpi {} > ser-mpi{}_{}_{}.txt".format(rank,s,mode,mode,rank,s))
      os.system(" mv *.txt /home/ubuntu/results_sc2022/lulesh.jl/ser-mpi-weak-scaling/")


itercount=100
printfun(1, [96],itercount)
printfun(8, [96],itercount)
printfun(27, [96],itercount)
printfun(64, [96],itercount)
