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
    include("test_averaging.jl")
    include("test_spectrum.jl")
    include("test_configuration.jl")
    include("test_pipeline.jl")

    @testset "quality" begin
        using Aqua
        Aqua.test_all(DeterministicSequentialEmission)

        using ExplicitImports
        @test check_no_implicit_imports(DeterministicSequentialEmission) === nothing
        @test check_no_stale_explicit_imports(DeterministicSequentialEmission) === nothing

        using JET
        using JET: AnyFrameModule
        using Measurements: Measurements
        using DrWatson: DrWatson
        # Measurements.jl keeps per-term derivatives in a dictionary whose iteration JET
        # cannot follow, and DrWatson's savename takes an empty `expand` default whose element
        # type JET cannot resolve; the reports raised there are about those dependencies'
        # internals, not about this package, so they are excluded by frame rather than
        # silenced wholesale.
        JET.test_package(
            DeterministicSequentialEmission;
            ignored_modules = (AnyFrameModule(Measurements), AnyFrameModule(DrWatson)),
        )
    end
end

isempty(SKIPPED) || @info "testsets skipped for absent input data" SKIPPED
