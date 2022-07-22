#!/usr/bin/python3.8
import os
import pathlib
scriptdir = pathlib.Path(__file__).parent.resolve()

def printfun(rank, blocklist,itercount):
  os.chdir(str(scriptdir)+"/../../")
  for s in blocklist:
    for mode in ["","--enzyme"]:
      os.system("./mpiexecjl  -bind-to socket --project -np {}  julia --project examples/benchmark.jl -s  --mpi {} {} {} > ser-mpi{}_{}_{}.txt".format(rank,mode,s,itercount, mode,rank,s))
      os.system("mv *.txt bench/ser-mpi-weak-scaling/")
  os.chdir(scriptdir)


itercount=10
printfun(1, [96],itercount)
printfun(8, [96],itercount)
printfun(27, [96],itercount)
printfun(64, [96],itercount)
