using QuadGK: quadgk

@testset "centre-of-mass spectra are normalized" begin
    for T in (0.4, 0.8, 1.2, 2.2)
        area, _ = quadgk(
            ε -> centre_of_mass_spectrum(ConstantCrossSection(), 0.0, T, ε),
            0,
            Inf;
            rtol = 1e-12,
        )
        @test area ≈ 1 atol = 1e-10

        # First moment of the constant cross section is 2T exactly.
        moment, _ = quadgk(
            ε -> ε * centre_of_mass_spectrum(ConstantCrossSection(), 0.0, T, ε),
            0,
            Inf;
            rtol = 1e-12,
        )
        @test moment ≈ 2T atol = 1e-9
    end

    for (α, T) in ((0.2, 0.8), (0.6, 1.2), (1.22, 1.8))
        area, _ = quadgk(
            ε -> centre_of_mass_spectrum(VariableCrossSection(), α, T, ε),
            0,
            Inf;
            rtol = 1e-12,
        )
        @test area ≈ 1 atol = 1e-10

        # And its first moment is the ⟨ε⟩ the sequential equation was solved with.
        moment, _ = quadgk(
            ε -> ε * centre_of_mass_spectrum(VariableCrossSection(), α, T, ε),
            0,
            Inf;
            rtol = 1e-12,
        )
        @test moment ≈ mean_neutron_energy(VariableCrossSection(), α, T) atol = 1e-9
    end

    # Non-positive arguments give zero rather than a NaN or a negative density.
    @test centre_of_mass_spectrum(ConstantCrossSection(), 0.0, 1.0, 0.0) == 0
    @test centre_of_mass_spectrum(VariableCrossSection(), 0.5, 1.0, -1.0) == 0
    @test centre_of_mass_spectrum(VariableCrossSection(), 0.5, 0.0, 1.0) == 0
end

@testset "the laboratory transformation preserves normalization" begin
    # The kernel width identity makes N(E) inherit unit area from φ exactly, for any E_f. This
    # is the property the transformation rests on, so it is checked rather than assumed.
    cases = [
        (cross_section, α, T, E_f) for
        cross_section in (ConstantCrossSection(), VariableCrossSection()) for
        α in (0.0, 0.35, 1.2) for T in (0.5, 1.1, 1.9) for
        E_f in (0.4, 1.05, 2.0) if !(cross_section isa ConstantCrossSection && α != 0.0)
    ]
    @test length(cases) == 36
    areas = [
        first(
            quadgk(
                E -> laboratory_spectrum(cross_section, α, T, E_f, E),
                0,
                Inf;
                rtol = 1e-11,
            ),
        ) for (cross_section, α, T, E_f) in cases
    ]
    @test maximum(abs(area - 1) for area in areas) <= 1e-8
end

@testset "laboratory moments follow the kinematics" begin
    # ⟨E⟩ = E_f + ⟨ε⟩ and ⟨E²⟩ = ⟨ε²⟩ + (10/3) E_f ⟨ε⟩ + E_f², both exact consequences of
    # isotropic emission from a moving source.
    for (cross_section, α) in
        ((ConstantCrossSection(), 0.0), (VariableCrossSection(), 0.45))
        for T in (0.6, 1.4), E_f in (0.5, 1.6)
            first_moment, _ = quadgk(
                E -> E * laboratory_spectrum(cross_section, α, T, E_f, E),
                0,
                Inf;
                rtol = 1e-11,
            )
            second_moment, _ = quadgk(
                E -> E^2 * laboratory_spectrum(cross_section, α, T, E_f, E),
                0,
                Inf;
                rtol = 1e-11,
            )
            ε₁, _ = quadgk(
                ε -> ε * centre_of_mass_spectrum(cross_section, α, T, ε),
                0,
                Inf;
                rtol = 1e-12,
            )
            ε₂, _ = quadgk(
                ε -> ε^2 * centre_of_mass_spectrum(cross_section, α, T, ε),
                0,
                Inf;
                rtol = 1e-12,
            )
            @test first_moment ≈ E_f + ε₁ rtol = 1e-7
            @test second_moment ≈ ε₂ + (10 / 3) * E_f * ε₁ + E_f^2 rtol = 1e-7
        end
    end
end

@testset "the closed form is the quadrature" begin
    # The substituted form is uniform in s over a band of fixed length, and is the reference
    # the closed form must reproduce. Checked across the whole energy range, not just the peak.
    reference(cross_section, α, T, E_f, E) =
        first(
            quadgk(
                s -> centre_of_mass_spectrum(cross_section, α, T, s^2),
                abs(sqrt(E) - sqrt(E_f)),
                sqrt(E) + sqrt(E_f);
                rtol = 1e-13,
            ),
        ) / (2 * sqrt(E_f))

    cases = [
        (cross_section, α, T, E_f, E) for (cross_section, α) in
        ((ConstantCrossSection(), 0.0), (VariableCrossSection(), 0.7)) for
        T in (0.5, 1.3) for E_f in (0.6, 1.4) for E in (0.01, 0.1, 1.0, 3.0, 8.0, 15.0)
    ]
    @test length(cases) == 48
    relative_error(got, expected) = abs(got - expected) / max(abs(got), abs(expected))
    @test maximum(
        relative_error(laboratory_spectrum(case...), reference(case...)) for case in cases
    ) <= 1e-9

    # The variable form collapses onto the constant one as α → 0.
    for T in (0.7, 1.5), E_f in (0.8,), E in (0.5, 2.0, 6.0)
        @test laboratory_spectrum(VariableCrossSection(), 0.0, T, E_f, E) ≈
              laboratory_spectrum(ConstantCrossSection(), 0.0, T, E_f, E) rtol = 1e-12
    end
end

@testset "the tail is evaluated without cancellation" begin
    # Taken over the lower incomplete gamma the difference cancels catastrophically far out in
    # the tail. These are twenty decades below the peak and must still be right.
    T, E_f, α = 1.0, 1.0, 0.0
    reference(E) =
        first(
            quadgk(
                s -> centre_of_mass_spectrum(ConstantCrossSection(), α, T, s^2),
                abs(sqrt(E) - sqrt(E_f)),
                sqrt(E) + sqrt(E_f);
                rtol = 1e-14,
            ),
        ) / (2 * sqrt(E_f))

    for E in (12.0, 15.0, 20.0, 25.0)
        got = laboratory_spectrum(ConstantCrossSection(), α, T, E_f, E)
        @test got > 0
        @test got ≈ reference(E) rtol = 1e-6
    end

    # Monotonically falling out there, with no cancellation noise turning it negative.
    values =
        [laboratory_spectrum(ConstantCrossSection(), 0.0, 1.0, 1.0, E) for E in 10:2:30]
    @test all(>(0), values)
    @test issorted(values; rev = true)
end

@testset "chain spectra average over the solved sequences" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    cross_section = VariableCrossSection()

    fragmentation = fragmentation_at(134, 52)
    balance = energetics(table, CF252, fragmentation, 180.0)
    method = ExcitationEnergyRatio(
        SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )
    split = partition(method, fragmentation, value(balance.TXE))
    E_light, E_heavy = fragment_energy_per_nucleon(fragmentation, 180.0)

    chain = emission_chain(table, model, cross_section, fragmentation.heavy, split.heavy)
    solved = count(is_solved, chain.sequences)
    @test solved > 0

    # An average of unit-area spectra is a unit-area spectrum.
    area, _ =
        quadgk(E -> chain_spectrum(chain, cross_section, E_heavy, E), 0, Inf; rtol = 1e-10)
    @test area ≈ 1 atol = 1e-7

    cm_area, _ = quadgk(
        ε -> chain_centre_of_mass_spectrum(chain, cross_section, ε),
        0,
        Inf;
        rtol = 1e-10,
    )
    @test cm_area ≈ 1 atol = 1e-7

    # Its first moment is the sequence average the solver already reports.
    cm_moment, _ = quadgk(
        ε -> ε * chain_centre_of_mass_spectrum(chain, cross_section, ε),
        0,
        Inf;
        rtol = 1e-11,
    )
    @test cm_moment ≈ sequence_average(chain, :mean_neutron_energy) rtol = 1e-8

    # And the laboratory first moment is that plus the kinetic energy per nucleon.
    lab_moment, _ = quadgk(
        E -> E * chain_spectrum(chain, cross_section, E_heavy, E),
        0,
        Inf;
        rtol = 1e-11,
    )
    @test lab_moment ≈ E_heavy + cm_moment rtol = 1e-7

    # A chain that solved nothing has no spectrum, rather than a zero-divided one.
    barren = emission_chain(table, model, cross_section, fragmentation.heavy, 0.0)
    @test count(is_solved, barren.sequences) == 0
    @test chain_spectrum(barren, cross_section, E_heavy, 2.0) == 0
    @test chain_centre_of_mass_spectrum(barren, cross_section, 2.0) == 0
end

@testset "pair weighting" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    cross_section = VariableCrossSection()

    fragmentation = fragmentation_at(140, 54)
    balance = energetics(table, CF252, fragmentation, 175.0)
    method = ExcitationEnergyRatio(
        SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )
    split = partition(method, fragmentation, value(balance.TXE))
    light = emission_chain(table, model, cross_section, fragmentation.light, split.light)
    heavy = emission_chain(table, model, cross_section, fragmentation.heavy, split.heavy)

    for weighting in (MultiplicityWeighting(), ExcitationWeighting())
        weights = pair_weights(weighting, light, heavy, split)
        @test weights !== nothing
        @test sum(weights) ≈ 1 atol = 1e-12
        @test all(>=(0), weights)
    end

    # The multiplicity weighting follows ν, the excitation weighting follows E*. They are
    # close by construction but not equal, which is why both are offered.
    ν = pair_weights(MultiplicityWeighting(), light, heavy, split)
    excitation = pair_weights(ExcitationWeighting(), light, heavy, split)
    @test ν[1] ≈
          neutron_multiplicity(light) /
          (neutron_multiplicity(light) + neutron_multiplicity(heavy)) rtol = RTOL
    @test excitation[2] ≈ split.heavy / sum(split) rtol = RTOL

    # A fragmentation that emits nothing has no pair spectrum to weight.
    barren = emission_chain(table, model, cross_section, fragmentation.light, 0.0)
    @test pair_weights(MultiplicityWeighting(), barren, barren, split) === nothing
end

@testset "a pair sum halves the yield at mass symmetry" begin
    # The fragment yield counts both fragments of a symmetric split on the heavy side, so a sum
    # over heavy fragments would give that split twice the weight of an asymmetric one.
    @test DeterministicSequentialEmission._pair_weight(fragmentation_at(126, 49), 3.0) ==
          1.5
    @test DeterministicSequentialEmission._pair_weight(fragmentation_at(126, 50), 3.0) ==
          1.5
    @test DeterministicSequentialEmission._pair_weight(fragmentation_at(127, 50), 3.0) ==
          3.0
end

@testset "total-spectrum weights" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    cross_section = VariableCrossSection()
    fragmentation = fragmentation_at(134, 52)
    balance = energetics(table, CF252, fragmentation, 180.0)
    method = ExcitationEnergyRatio(
        SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )
    split = partition(method, fragmentation, value(balance.TXE))
    light = emission_chain(table, model, cross_section, fragmentation.light, split.light)
    heavy = emission_chain(table, model, cross_section, fragmentation.heavy, split.heavy)
    ν_pair = neutron_multiplicity(light) + neutron_multiplicity(heavy)
    @test ν_pair > 0
    # Per fission the yield stands; per neutron it is scaled by the pair multiplicity.
    @test total_weight(PerFissionAverage(), 2.0, light, heavy) == 2.0
    @test total_weight(PerNeutronAverage(), 2.0, light, heavy) ≈ 2.0 * ν_pair rtol = RTOL
end
