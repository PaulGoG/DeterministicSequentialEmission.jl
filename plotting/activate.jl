# Activate and instantiate the plotting environment, silently.
#
# Kept apart from the package: the figures need a plotting stack the model itself does not, and
# nothing in src/ depends on it. Below Julia 1.11 the [sources] entry is ignored, so the package
# is developed by path instead.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
if VERSION < v"1.11"
    Pkg.develop(; path = dirname(@__DIR__), io = devnull)
end
Pkg.instantiate(; io = devnull)
