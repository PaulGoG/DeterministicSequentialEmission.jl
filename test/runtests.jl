using Test
using TOML: TOML
using Measurements: Measurement, measurement, value, uncertainty
using DeterministicSequentialEmission

const DATA = joinpath(@__DIR__, "..", "data")
const REFERENCE = joinpath(@__DIR__, "references", "Cf252_sf")

include("physics.jl")

@testset "DeterministicSequentialEmission" begin
    include("test_nuclides.jl")
    include("test_tables.jl")
    include("test_masses.jl")
    include("test_charge.jl")
    include("test_wahl.jl")
    include("test_domain.jl")
    include("test_leveldensity.jl")
    include("test_energetics.jl")
    include("test_partition.jl")
    include("test_sequences.jl")
    include("test_yields.jl")
    include("test_averaging.jl")
    include("test_spectrum.jl")
    include("test_configuration.jl")

    @testset "quality" begin
        using Aqua
        Aqua.test_all(DeterministicSequentialEmission)

        using ExplicitImports
        @test check_no_implicit_imports(DeterministicSequentialEmission) === nothing
        @test check_no_stale_explicit_imports(DeterministicSequentialEmission) === nothing

        using JET
        using JET: AnyFrameModule
        using Measurements: Measurements
        # Measurements.jl keeps per-term derivatives in a dictionary whose iteration JET
        # cannot follow; the reports raised there are about that dependency's internals, not
        # about this package, so they are excluded by frame rather than silenced wholesale.
        JET.test_package(
            DeterministicSequentialEmission;
            ignored_modules = (AnyFrameModule(Measurements),),
        )
    end
end
