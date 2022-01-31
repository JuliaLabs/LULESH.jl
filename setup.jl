#!/bin/env julia
import Pkg
Pkg.activate(".")

import Enzyme
dir = joinpath(dirname(dirname(pathof(Enzyme))), "deps")

@info "Building Enzyme from master"

run(`$(Base.julia_cmd()) --project=$(dir) -e 'import Pkg; Pkg.instantiate()'`)
run(`$(Base.julia_cmd()) --project=$(dir) $(dir)/build_local.jl`)

import MPI
MPI.install_mpiexecjl(;destdir=".", force=true)
