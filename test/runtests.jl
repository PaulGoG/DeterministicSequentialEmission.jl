using Test
using TOML: TOML
using Measurements: Measurement, measurement, value, uncertainty
using FissionFragmentsDomain
using DeterministicSequentialEmission

const DATA = joinpath(@__DIR__, "..", "data")
const REFERENCE = joinpath(@__DIR__, "references", "Cf252_sf")

include("physics.jl")

@testset "DeterministicSequentialEmission" begin
    include("test_partition.jl")
    include("test_manifest.jl")
    include("test_sequences.jl")
    include("test_yields.jl")
    include("test_averaging.jl")
    include("test_spectrum.jl")
    include("test_configuration.jl")
    include("test_pipeline.jl")
    include("test_products.jl")

    @testset "quality" begin
        using Aqua
        Aqua.test_all(DeterministicSequentialEmission)

        using ExplicitImports
        @test check_no_implicit_imports(DeterministicSequentialEmission) === nothing
        @test check_no_stale_explicit_imports(DeterministicSequentialEmission) === nothing

        using JET
        # A report counts only where the problem is detected in a frame of this package. The
        # rest is about the internals of what it calls: Measurements.jl keeps per-term
        # derivatives in a dictionary whose iteration JET cannot follow, DrWatson's savename
        # takes an empty `expand` default whose element type JET cannot resolve, and the JET
        # that Julia 1.11 admits (0.9) raises two reports inside Base's twice-precision ranges
        # and its broadcasting machinery that JET 0.12 on Julia 1.13 does not.
        JET.test_package(
            DeterministicSequentialEmission;
            target_modules = (DeterministicSequentialEmission,),
        )
    end
end

isempty(SKIPPED) || @info "testsets skipped for absent input data" SKIPPED
