@testset "fragmentation domain" begin
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

@testset "symmetric splits" begin
    domain = cf252_domain()
    symmetric = filter(is_symmetric, domain.entries)

    # A_H = 126 is exactly symmetric for 252Cf, so five entries qualify — one per charge.
    @test length(symmetric) == 5
    for entry in symmetric
        @test entry.heavy.A == entry.light.A == 126
    end

    # Every initial fragment appears once in the fragment view, including the symmetric mass,
    # whose two fragments share a nuclide only when their charges also coincide.
    listed = fragments(domain)
    @test allunique(first.(listed))
    @test issorted([(nuclide.A, nuclide.Z) for (nuclide, _) in listed])
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
end
