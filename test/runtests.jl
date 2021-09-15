using LULESH
using Test
using MPI

function run_example(name, args; nranks=nothing)
    example = joinpath(@__DIR__, "..", "examples", name)

    @testset "$(basename(example)) $args $nranks" begin
        cmd = `$(Base.julia_cmd()) --startup-file=no $(example) $(args)`
        if nranks !== nothing
            mpiexecjl = joinpath(dirname(pathof(MPI)), "..", "bin", "mpiexecjl")
            cmd = `$(mpiexecjl) -n $(nranks) --oversubscribe $(cmd)`
        end
        @debug "Testing $example" cmd
        @test success(pipeline(cmd, stderr=stderr))
    end
end

@testset "examples" begin
    @testset "benchmark.jl" begin
        run_example("benchmark.jl", `-s 45`)
        run_example("benchmark.jl", `-s 45`; nranks=1)
        run_example("benchmark.jl", `-s 45`; nranks=8)
    end
end
