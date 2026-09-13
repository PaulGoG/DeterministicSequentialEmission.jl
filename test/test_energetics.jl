@testset "compound nucleus excitation" begin
    table = mass_table()

    # Spontaneous fission starts from an unexcited compound nucleus.
    @test value(compound_nucleus_excitation(table, CF252)) == 0.0
    @test uncertainty(compound_nucleus_excitation(table, CF252)) == 0.0

    # Thermal neutron capture leaves essentially Sₙ of the compound nucleus.
    pu = neutron_induced_fission(Nuclide(94, 239), 2.53e-8)
    E = compound_nucleus_excitation(table, pu)
    @test value(E) ≈ value(neutron_separation_energy(table, Nuclide(94, 240))) + 2.53e-8
    @test 5.0 < value(E) < 8.0
end

@testset "total excitation energy" begin
    table = mass_table()
    heavy = Nuclide(52, 134)
    Q = q_value(table, CF252, heavy)
    E_compound = compound_nucleus_excitation(table, CF252)

    # TXE = Q + E*_CN − TKE, and TKE enters without an uncertainty because it is the swept
    # independent variable of the model rather than a measurement.
    TXE = total_excitation_energy(Q, E_compound, 180.0)
    @test value(TXE) ≈ value(Q) - 180.0
    @test uncertainty(TXE) ≈ uncertainty(Q)

    # Above Q the configuration is unphysical and is reported as absent, not as a negative.
    @test total_excitation_energy(Q, E_compound, value(Q) + 1.0) === nothing
end

@testset "energetics of one configuration" begin
    table = mass_table()
    entry = fragmentation_at(134, 52)

    result = energetics(table, CF252, entry, 180.0)
    @test result isa Energetics
    @test result.TKE == 180.0
    @test value(result.TXE) ≈ value(result.Q) - 180.0
    @test occursin("A_H=134", sprint(show, result))

    @test energetics(table, CF252, entry, 400.0) === nothing
end

POLARIZATION_AVAILABLE && @testset "sweeping the domain" begin
    table = mass_table()
    domain = cf252_domain()
    accepted, excluded = sweep_energetics(table, CF252, domain, CF252_TKE)

    # The declared grid: 245 fragmentations × 51 kinetic energies.
    @test length(accepted) + length(excluded) == 245 * 51 == 12495

    # Nothing is dropped silently — every exclusion carries its reason, so the coverage of a
    # run can be audited from its own output.
    summary = exclusion_summary(excluded)
    @test sum(values(summary)) == length(excluded)
    @test issubset(keys(summary), Set([:missing_mass, :nonpositive_q, :nonpositive_txe]))

    # Every 252Cf fragmentation in this range has a tabulated mass and a positive Q, so the
    # only reason a configuration drops out is TXE ≤ 0 at high kinetic energy.
    @test get(summary, :missing_mass, 0) == 0
    @test get(summary, :nonpositive_q, 0) == 0

    # TXE falls monotonically with TKE at fixed Q, so for each fragmentation the excluded
    # kinetic energies must form an upper tail of the grid — never a hole in the middle.
    dropped = Dict{Tuple{Int, Int}, Vector{Float64}}()
    for exclusion in excluded
        key = (exclusion.fragmentation.heavy.A, exclusion.fragmentation.heavy.Z)
        push!(get!(dropped, key, Float64[]), exclusion.TKE)
    end
    for (key, energies) in dropped
        sort!(energies)
        @test energies == collect(CF252_TKE)[(end - length(energies) + 1):end]
    end
    @test !isempty(dropped)

    for result in accepted
        @test value(result.TXE) > 0
        @test uncertainty(result.TXE) > 0
    end
end

@testset "kinetic energy of a fragment" begin
    # KE(A) = TKE (A₀ − A)/A₀: the light fragment takes the larger share, and the two halves
    # of a pair sum back to TKE.
    heavy = Nuclide(52, 134)
    light = complementary_fragment(CF252, heavy)
    KE_heavy = kinetic_energy(CF252, heavy, 180.0)
    KE_light = kinetic_energy(CF252, light, 180.0)
    @test KE_light > KE_heavy
    @test KE_heavy + KE_light ≈ 180.0
    @test kinetic_energy(CF252, Nuclide(49, 126), 180.0) ≈ 90.0
end
