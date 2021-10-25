using FluxRM

function __get_treedict!(dict, path...)
    next_dict = dict
    for node in path
        next_dict = get!(()-> Dict{String, Any}(), next_dict, node)
    end
    return next_dict
end

function __set_io_path(spec::FluxRM.JobSpec.Jobspec, iotype, name, path)
    system = spec.attributes.system
    if system.shell === nothing
        system.shell = Dict{String, Any}()
    end
    io = __get_treedict!(system.shell, "options", iotype, name)

    io["type"] = "file"
    io["path"] = path
end

function juliaspec(args, dir; num_nodes=1, num_tasks_per_node=8, cores_per_task=1)
    num_tasks = num_nodes*num_tasks_per_node
    cmd = `$(Base.julia_cmd()) $(args)`
    jobspec = FluxRM.JobSpec.from_command(cmd; num_nodes, num_tasks, cores_per_task)
    system = jobspec.attributes.system
    system.cwd = dir
    system.environment = Dict(
        "JULIA_PROJECT" => dir,
        "JULIA_NUM_THREADS" => cores_per_task,
        "JULIA_EXCLUSIVE" => 1
    )
    __set_io_path(jobspec, "output", "stderr", "flux-{{id}}.err")
    __set_io_path(jobspec, "output", "stdout", "flux-{{id}}.out")
    @assert FluxRM.JobSpec.validate(jobspec, 1)

    # FluxRM.JSON3.pretty(FluxRM.JSON3.write(jobspec))

    jobspec
end

# Notes
# - Flux treats hyper-threads as a single core

function nodes()
    rpc = fetch(FluxRM.RPC(Flux(), "resource.status", nodeid=0))
    R = first(rpc.R.execution.R_lite)

    hosts = FluxRM.IDSet(R.rank) 
    Int(length(hosts))
end

const N = nodes()
const workdir =  realpath(joinpath(@__DIR__, "..", "examples"))

@info "Launching Jobs in " workdir

let flux = Flux()
    for i in 0:floor(Int,log2(N))
        n = 2^i
        jobspec = juliaspec(`benchmark.jl -s 45 --mpi`, workdir, num_nodes=n)
        sub = FluxRM.submit(flux, jobspec)
        job = FluxRM.Job(sub)
        @info "Launched" jobid = FluxRM.encode(job) n
    end
end
