#!/usr/bin/python3.8
import os

def printfun(rank, blocklist,itercount,mode):
     for s in blocklist:
        os.system("echo -n {},{}, >> results.txt".format(rank,s))
        os.system("grep  Elapsed ser-mpi{}_{}_{}.txt >> results.txt".format(mode,rank, s))
     os.system("sed -i \"s/Elapsed time         = //g\" results.txt")
     os.system("sed -i \"s/ (s)//g\" results.txt")

os.system("rm results.txt")
itercount=10
for mode in ["","--enzyme"]:
  printfun(1, [192],itercount,mode)
  printfun(8, [96],itercount,mode)
  printfun(27, [64],itercount,mode)
  printfun(64, [48],itercount,mode)
