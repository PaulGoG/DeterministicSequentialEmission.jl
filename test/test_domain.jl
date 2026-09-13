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

    listed = fragments(domain)
    @test allunique(first.(listed))
    @test issorted([(nuclide.A, nuclide.Z) for (nuclide, _) in listed])

    # At the symmetric mass a nuclide is reached by two different fragmentations — once as the
    # heavy fragment, once as the light complement of its charge mirror — with two different
    # probabilities. The fragment view carries their sum, so the census is complete.
    weights = Dict(nuclide => weight for (nuclide, weight) in listed)
    at_126 = Dict(entry.heavy.Z => entry.probability for entry in symmetric)
    for Z in 47:51
        @test weights[Nuclide(Z, 126)] ≈ at_126[Z] + at_126[98 - Z]
    end

    # Concretely: Z = 50 is produced with 0.4058 as a heavy fragment and 0.0143 as the light
    # complement of Z = 48. Keeping only the first, as the previous implementation did, gave it
    # its mirror's weight — a factor of 28 — and lost 40 % of the mass row.
    @test weights[Nuclide(50, 126)] ≈ 0.4201 atol = 1e-4
    @test sum(weights[Nuclide(Z, 126)] for Z in 47:51) ≈ 2 * sum(values(at_126))

    # Away from symmetry each nuclide belongs to exactly one fragmentation, so nothing is
    # summed and the weight is that fragmentation's own.
    @test weights[Nuclide(52, 134)] ≈
          only(e.probability for e in domain.entries if e.heavy == Nuclide(52, 134))
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
