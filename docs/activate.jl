using Pkg
# FissionFragmentsDomain resolves from its repository URL. The git command line honours the
# user's url.insteadOf rewrites and SSH agent, which libgit2 does not; an explicit setting wins.
haskey(ENV, "JULIA_PKG_USE_CLI_GIT") || (ENV["JULIA_PKG_USE_CLI_GIT"] = "true")
Pkg.activate(@__DIR__; io = devnull)
# Resolve before instantiating: the package's own [compat] changes as it grows, and a manifest
# resolved against an older one would otherwise fail to load rather than update itself.
Pkg.resolve(; io = devnull)
Pkg.instantiate(; io = devnull)
