using Test
using DeterministicSequentialEmission

@testset "DeterministicSequentialEmission" begin
    @testset "quality" begin
        using Aqua
        Aqua.test_all(DeterministicSequentialEmission)

        using ExplicitImports
        @test check_no_implicit_imports(DeterministicSequentialEmission) === nothing
        @test check_no_stale_explicit_imports(DeterministicSequentialEmission) === nothing

        using JET
        JET.test_package(DeterministicSequentialEmission)
    end
end
