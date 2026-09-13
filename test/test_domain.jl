POLARIZATION_AVAILABLE && @testset "fragmentation domain" begin
    domain = cf252_domain()

    # 49 heavy masses × 5 charge numbers, the configuration of the archived run.
    @test length(domain) == 245
    @test domain.heavy_masses == 126:174
    @test domain.charges_per_mass == 5
    @test mass_range(domain) == 78:174

    # Sorted as whole rows, by heavy mass then heavy charge. The prototype sorted the mass
    # and charge columns independently and then looked values up by the pair, which can
    # report a mass against a charge that never occurred with it.
    keys = [(entry.heavy.A, entry.heavy.Z) for entry in domain]
    @test issorted(keys)
    @test allunique(keys)

    for entry in domain
        @test entry.heavy.A + entry.light.A == 252
        @test entry.heavy.Z + entry.light.Z == 98
        @test entry.probability > 0
    end

    # The probability is a property of the pair: the light fragment's Gaussian argument is
    # the negative of the heavy one's, so both fragments carry the same p(Z,A).
    polarization = cf252_polarization()
    entry = first(domain)
    Zₚ = most_probable_charge(
        CF252,
        entry.heavy.A,
        polarization_at(polarization, entry.heavy.A),
    )
    rms = width_at(polarization, entry.heavy.A)
    @test entry.probability ≈ charge_probability(entry.heavy.Z, Zₚ, rms)
    @test entry.probability ≈ charge_probability(entry.light.Z, 98 - Zₚ, rms)
end

POLARIZATION_AVAILABLE && @testset "symmetric splits" begin
    domain = cf252_domain()
    symmetric = filter(is_symmetric, domain.entries)

    # A_H = 126 is exactly symmetric for 252Cf, so five entries qualify — one per charge.
    @test length(symmetric) == 5
    for entry in symmetric
        @test entry.heavy.A == entry.light.A == 126
    end

    # The fragment view is the nuclide set, each once — what per-nuclide quantities are
    # memoised over. No probability: a nuclide's weight depends on the question being asked.
    listed = fragments(domain)
    @test allunique(listed)
    @test issorted([(nuclide.A, nuclide.Z) for nuclide in listed])
    @test length(listed) == 485
    @test Nuclide(50, 126) in listed
    @test eltype(listed) == Nuclide

    # It agrees with the charge-probability table, the other per-nuclide view.
    @test Set(listed) ==
          Set(keys(fragment_charge_probabilities(domain, cf252_polarization())))
end

POLARIZATION_AVAILABLE && @testset "distinct splits" begin
    domain = cf252_domain()
    splits = distinct_splits(domain)

    # Away from mass symmetry a split is labelled once, so only A₀/2 loses entries. The five
    # charge labels at A = 126 describe three physical splits: {47,51}, {48,50} and {49,49}.
    @test length(splits) == length(domain) - 2
    symmetric = [entry for entry in splits if entry.heavy.A == 126]
    @test length(symmetric) == 3
    @test sort([entry.heavy.Z for entry in symmetric]) == [49, 50, 51]

    # Every physical pair appears exactly once, in one labelling or the other.
    pairs = Set(Set([entry.heavy.Z, entry.light.Z]) for entry in symmetric)
    @test length(pairs) == 3

    # A split whose mirror charge falls outside the retained window has no second labelling and
    # must survive whichever side names it. With the mean polarization Zₚ(126) = 48.5, so the
    # window 46:50 is not symmetric about Z₀/2 = 49 — exactly that case.
    offset = fragmentation_domain(CF252, mean_charge_polarization(), 126:174)
    offset_splits = distinct_splits(offset)
    @test length(offset_splits) == length(offset) - 1
    @test sort([e.heavy.Z for e in offset_splits if e.heavy.A == 126]) == [46, 47, 49, 50]
end

@testset "domain validation" begin
    polarization = mean_charge_polarization()

    # The range is validated against the compound nucleus, which is what splits. The
    # prototype validates against the target, which agrees only when the bound sits at A₀/2.
    pu = neutron_induced_fission(Nuclide(94, 239), 2.53e-8)
    @test symmetric_mass(pu) == 120.0
    @test_throws ArgumentError fragmentation_domain(pu, polarization, 119:160)
    @test length(fragmentation_domain(pu, polarization, 120:160)) == 41 * 5

    @test_throws ArgumentError fragmentation_domain(CF252, polarization, 125:174)
    @test_throws ArgumentError fragmentation_domain(CF252, polarization, 126:252)
    @test_throws ArgumentError fragmentation_domain(
        CF252,
        polarization,
        126:174;
        charges_per_mass = 4,
    )

    # A degenerate single-mass range. The prototype rejects A_H_min >= A_H_max outright.
    @test_throws ArgumentError fragmentation_domain(CF252, polarization, 130:130)

    # Enough charge numbers per mass and the window walks off the end of the charge range,
    # asking for a fragment with no protons or with all of them.
    @test_throws ArgumentError fragmentation_domain(
        CF252,
        polarization,
        126:174;
        charges_per_mass = 199,
    )
end
