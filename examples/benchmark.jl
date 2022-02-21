using LLVM
LLVM.clopts("-memdep-block-scan-limit=70000")
LLVM.clopts("-dse-memoryssa-walklimit=10000")
LLVM.clopts("-attributor-max-iterations=128")
LLVM.clopts("-capture-tracking-max-uses-to-explore=256")

using Enzyme

Enzyme.API.printperf!(true)
Enzyme.API.printall!(true)
Enzyme.API.instname!(true)

Enzyme.API.inlineall!(true)
Enzyme.API.maxtypeoffset!(1024)
isdefined(Enzyme.API, :strictAliasing!) && Enzyme.API.strictAliasing!(true)
isdefined(Enzyme.API, :typeWarning!) &&  Enzyme.API.typeWarning!(false)
Enzyme.API.looseTypeAnalysis!(true)

mutable struct Data
   commDataSend::Vector{Float64}
end

function free(buf)
  return nothing
end

mutable struct Something
   x::Int64
end
function Isend()
    req = Something(2)
    finalizer(free, req)
    return req
end

function fooSend(domain, fields, dx)
	 offset = 2
         for field in fields
            for i in 0:(dx-1)
               domain.commDataSend[offset+i + 1] = field[30+3*i + 1]
            end
            offset += 2
         end
         req = Isend()
    return nothing 
end
function foo(domain, domx, dx, dy, dz)
      fields = (domx, domx, domx, domx, domx, domx)
      fooSend(domain, fields,
		dx)

    return nothing
end

function main(enzyme)
     
	domain = Data(Vector{Float64}(undef, 2+11+31))
        shadowDomain = Data(Vector{Float64}(undef, 2+11+31))

   dx = 30 + 1
   dy = 30 + 1
   dz = 30 + 1
	domx = Vector{Float64}(undef, 3*30+31)
	sdomx = Vector{Float64}(undef, 3*30+31)

	if enzyme
            Enzyme.autodiff(foo, Duplicated(domain, shadowDomain), Duplicated(domx, sdomx), dx, dy, dz)
        else
            foo(domain, domx, dx, dy, dz)
        end
end

if !isinteractive()
    main(false)
    @show "ran primal"
    flush(stdout)
    main(true)
end
