# Pre-commit: format against .JuliaFormatter.toml with the version-bounded formatter, then run the tests.
#
#     julia check.jl                format in place, then test
#     julia check.jl --check        fail on formatting differences, do not rewrite, then test
#     julia check.jl --format-only  fail on formatting differences, do not rewrite, do not test
#
# The same gate the CI workflow applies, so a failure is reproducible locally.

const DIRECTORIES = ("src", "test", "scripts", "plotting", "docs")

let overwrite = !("--check" in ARGS || "--format-only" in ARGS)
    include(joinpath(@__DIR__, "formatter", "activate.jl"))
    using JuliaFormatter: format

    formatted = true
    for directory in DIRECTORIES
        path = joinpath(@__DIR__, directory)
        isdir(path) || continue
        formatted &= format(path; overwrite = overwrite)
    end

    if !formatted
        if overwrite
            @info "formatting applied"
        else
            @error "formatting differs from .JuliaFormatter.toml; run `julia check.jl`"
            exit(1)
        end
    end

    # The test job runs the suite; the formatting job stops here.
    "--format-only" in ARGS && exit(0)

    using Pkg
    # As in activate.jl: the URL dependency is fetched through the git command line.
    haskey(ENV, "JULIA_PKG_USE_CLI_GIT") || (ENV["JULIA_PKG_USE_CLI_GIT"] = "true")
    Pkg.activate(@__DIR__; io = devnull)
    Pkg.test()
end
