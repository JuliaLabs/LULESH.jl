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

function juliaspec(args, dir; num_nodes=1, num_tasks_per_node=6, cores_per_task=6)
    num_tasks = num_nodes*num_tasks_per_node
    cmd = `$(Base.julia_cmd()) $(args)`
    jobspec = FluxRM.JobSpec.from_command(cmd; num_nodes, num_tasks, cores_per_task)
    system = jobspec.attributes.system
    system.cwd = dir
    system.environment = Dict(
        "JULIA_PROJECT" => dir, 
        "OPENBLAS_NUM_THREADS" => "8" # HyperThreads
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

let flux = Flux()

    for i in 0:floor(Int,log2(N))
        n = 2^i
        for psize in (20,)
            jobspec = juliaspec(`-L setup.jl experiment.jl $psize`, realpath("experiment"), num_nodes=n)
            sub = FluxRM.submit(flux, jobspec)
            job = FluxRM.Job(sub)
            @info "Launched" jobid = FluxRM.encode(job) n psize
        end
    end
end
