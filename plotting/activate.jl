# Activate and instantiate the plotting environment, silently.
#
# Kept apart from the package: the figures need a plotting stack the model itself does not, and
# nothing in src/ depends on it. The figures read what scripts/run.jl wrote, so this environment
# does not carry the package either.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.instantiate(; io = devnull)
