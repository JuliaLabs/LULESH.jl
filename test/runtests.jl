
# HACK: work around Pkg.jl#2500
test_project = Base.active_project()
preferences_file = joinpath(dirname(@__DIR__), "LocalPreferences.toml")
test_preferences_file = joinpath(dirname(test_project), "LocalPreferences.toml")
if isfile(preferences_file) && !isfile(test_preferences_file)
    cp(preferences_file, test_preferences_file)
end

using LULESH
using Test
using MPI

if parse(Bool, get(ENV, "ENZYME_CI", "false"))
    using Enzyme_jll
    @info "Testing against" Enzyme_jll.libEnzyme
    args = `--enzyme`
else
    args = ``
end

MPI.install_mpiexecjl(; destdir = ".", force=true)

function run_example(name, args; nranks=nothing)
    example = joinpath(@__DIR__, "..", "examples", name)

    @testset "$(basename(example)) $args $nranks" begin
        cmd = `$(Base.julia_cmd()) --startup-file=no $(example) $(args)`
        if nranks !== nothing
            cmd = `./mpiexecjl --project=$(test_project) -n $(nranks) $(cmd)`
        end
        @debug "Testing $example $nranks" cmd
        @test success(pipeline(cmd, stderr=stderr))
    end
end

@testset "examples" begin
    @testset "benchmark.jl" begin
        run_example("benchmark.jl", `-s 45 $(args)`)
        run_example("benchmark.jl", `-s 45 --mpi $(args)`; nranks=1)
        run_example("benchmark.jl", `-s 45 --mpi $(args)`; nranks=8)
    end
end
