@testset "nuclides" begin
    @test Nuclide(98, 252).Z == 98
    @test neutron_number(Nuclide(98, 252)) == 154
    @test NEUTRON == Nuclide(0, 1)
    @test element_symbol(98) == "Cf"
    @test element_symbol(0) == "n"
    @test element_symbol(200) == "Z=200"

    @test_throws ArgumentError Nuclide(-1, 10)
    # A nuclide cannot hold fewer nucleons than protons.
    @test_throws ArgumentError Nuclide(50, 40)
end

@testset "fissioning systems" begin
    cf = spontaneous_fission(Nuclide(98, 252))
    @test cf.spontaneous
    @test cf.compound == cf.target
    @test iszero(cf.incident_energy)
    @test symmetric_mass(cf) == 126.0

    # The compound nucleus carries one more mass unit than the target, and both are kept:
    # the domain is built on the compound, a fission case is named by the target.
    pu = neutron_induced_fission(Nuclide(94, 239), 2.53e-8)
    @test pu.compound == Nuclide(94, 240)
    @test pu.target == Nuclide(94, 239)
    @test !pu.spontaneous
    @test symmetric_mass(pu) == 120.0

    @test_throws ArgumentError FissioningSystem(
        Nuclide(98, 252),
        Nuclide(98, 252),
        1.0,
        true,
    )
    @test_throws ArgumentError FissioningSystem(
        Nuclide(98, 253),
        Nuclide(98, 252),
        0.0,
        true,
    )
    @test_throws ArgumentError neutron_induced_fission(Nuclide(94, 239), -1.0)

    @test occursin("252Cf(SF)", sprint(show, cf))
    @test occursin("239Pu(n,f)", sprint(show, pu))
end
