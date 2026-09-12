using Test
using Measurements: value, uncertainty
using DeterministicSequentialEmission

const INPUT = joinpath(@__DIR__, "..", "input_data")
const REFERENCE = joinpath(@__DIR__, "references")

include("physics.jl")

@testset "DeterministicSequentialEmission" begin
    include("test_nuclides.jl")
    include("test_tables.jl")
    include("test_masses.jl")
    include("test_charge.jl")
    include("test_domain.jl")
    include("test_energetics.jl")

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
