# Running the benchmakrs

## Installation

In `benchmarks/`:

```
spack env activate .
spack concretize
despacktivate
spack env activate .
julia --project. -e "import Pkg; Pkg.instantiate()"
```


## Launching
```
spack env activate .
flux start
julia --project=. run.jl
```
