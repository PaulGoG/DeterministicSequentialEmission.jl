@testset "Wahl Zₚ systematics" begin
    # 235-U(n_th,f), the reaction the systematics is referenced to: the differences in eq. (17)
    # all vanish except the excitation term, so the peak parameters reduce to their intercepts.
    U236 = WahlSystematics(92, 236, 6.551)
    @test U236.σZ140 ≈ 0.566
    @test U236.ΔZ140 ≈ -0.487
    @test U236.σZ50 ≈ 0.356
    @test U236.SL50 ≈ 0.191
    @test U236.ΔZmax ≈ 0.699

    # Region boundaries, eq. (10) and Fig. 18, p. 31: the wings start at 77, the near-symmetry
    # region spans roughly 106 to 130, and the Z = 50 crossing sits near 125.
    @test U236.B2 ≈ 77.0
    @test 129 < U236.B4 < 132
    @test 110 < U236.Ba < 113
    @test 123 < U236.Bb < 126
    @test U236.B5 ≈ 236 - U236.B2
    @test U236.B6 ≈ 166.0

    # Charge conservation: the polarization is equal and opposite about the symmetric split, so
    # it vanishes there. Nothing in the construction imposes this — it follows from the straight
    # line between Ba and Bb being centred on A_F/2.
    @test wahl_polarization(U236, 118) ≈ 0.0 atol = 1e-12

    # ΔZ is continuous across every boundary the heavy side crosses: the legs are constructed to
    # meet, and a jump in the centroid of the charge distribution would be unphysical.
    for boundary in (U236.Bb, U236.B4, U236.B5, U236.B6)
        below = wahl_polarization(U236, boundary - 1e-8)
        above = wahl_polarization(U236, boundary + 1e-8)
        @test below ≈ above atol = 1e-6
    end

    # σ_Z is not, and that is the model rather than an implementation slip. Eqs. (12c) and (12g)
    # impose σ_Z(50) as a flat value over the two legs adjoining the Z = 50 crossing while the
    # crossing itself keeps the peak value, and eq. (14c) returns the far wing to σ_Z(B5) rather
    # than continuing the wing slope. Fig. 18a, p. 31, draws both as steps.
    @test wahl_width(U236, U236.B4 - 1e-8) ≈ U236.σZ50
    @test wahl_width(U236, U236.B4 + 1e-8) > U236.σZ50
    @test wahl_width(U236, U236.Bb - 1e-8) > U236.σZ50
    @test wahl_width(U236, U236.B6 - 1e-8) < wahl_width(U236, U236.B6 + 1e-8)
    @test wahl_width(U236, U236.B6 + 1e-8) ≈ wahl_width(U236, U236.B5)

    # The Bb–B4 leg is constructed so that it reaches ΔZ_max exactly at the crossing.
    @test wahl_polarization(U236, U236.Bb) ≈ U236.ΔZmax atol = 1e-6

    # In the peak region the polarization is the straight line of eq. (11a) about A' = 140.
    @test wahl_polarization(U236, 140) ≈ -0.487
    @test wahl_polarization(U236, 150) ≈ -0.487 - 0.0080 * 10
    @test wahl_width(U236, 140) ≈ 0.566
    @test wahl_width(U236, 150) ≈ 0.566 - 0.0038 * 10

    # σ_Z narrows to σ_Z(50) through the near-symmetry region, the signature of the shell
    # closure: the charge is pinned at Z = 50 and the distribution about it is tighter.
    @test wahl_width(U236, (U236.Bb + U236.B4) / 2) ≈ U236.σZ50
    @test U236.σZ50 < U236.σZ140

    # 252-Cf, a heavier and more charged system: every difference in eq. (17) contributes.
    CF252 = WahlSystematics(98, 252, 0.0)
    @test CF252.σZ140 ≈ 0.566 + 0.0064 * 16 + 0.0109 * (0 - 6.551)
    @test CF252.ΔZ140 ≈ -0.487 + 0.0180 * 16 - 0.00203 * 16^2
    @test CF252.σZ50 ≈ 0.356 + 0.060 * 6
    @test wahl_polarization(CF252, 126) ≈ 0.0 atol = 1e-12

    # Outside the fitted ranges the constructor refuses rather than extrapolating quietly.
    @test_throws ArgumentError WahlSystematics(88, 236, 6.5)      # Z_F below 90
    @test_throws ArgumentError WahlSystematics(92, 260, 6.5)      # A_F above 252
    @test_throws ArgumentError WahlSystematics(92, 236, 12.0)     # past the low-energy branch
    @test_throws ArgumentError WahlSystematics(92, 236, -1.0)
    # And the light fragment is not its domain: ΔZ there is the complement, not a lookup.
    @test_throws ArgumentError wahl_polarization(U236, 100)
    @test_throws ArgumentError wahl_width(U236, 100)
end

@testset "Wahl systematics as a charge polarization" begin
    table = mass_table()
    system = neutron_induced_fission(Nuclide(92, 235), 2.53e-8)
    polarization = wahl_charge_polarization(table, system, 118:160)
    @test polarization !== nothing
    @test polarization_at(polarization, 140) ≈
          wahl_polarization(WahlSystematics(table, system), 140)
    @test width_at(polarization, 140) ≈ wahl_width(WahlSystematics(table, system), 140)

    # It is the same type a tabulated file produces, so every consumer is unchanged, and a mass
    # outside the tabulated range still falls back to the conventional means.
    @test polarization isa ChargePolarization
    @test polarization_at(polarization, 300) == -0.5
    @test width_at(polarization, 300) == 0.6
    @test occursin("LA-13928", sprint(show, MIME"text/plain"(), polarization.source))

    # The precursor excitation energy of a thermal capture is Sₙ of the compound nucleus, so the
    # systematics built from the system agrees with the one built from that number by hand.
    excitation = value(compound_nucleus_excitation(table, system))
    @test WahlSystematics(table, system).PE ≈ excitation
    @test WahlSystematics(92, 236, excitation).ΔZ140 ≈ WahlSystematics(table, system).ΔZ140

    # Spontaneous fission enters at zero excitation.
    @test WahlSystematics(table, spontaneous_fission(Nuclide(98, 252))).PE == 0.0

    @test_throws ArgumentError wahl_charge_polarization(
        WahlSystematics(92, 236, 6.5),
        Int[],
    )
end

POLARIZATION_AVAILABLE && @testset "Wahl systematics against the evaluated tables" begin
    # An external check, against numbers this package did not produce. The shipped tables are a
    # per-reaction least-squares fit of the same model — Fig. 18 quotes reduced χ² of 2.9 for
    # that against 7.9 for these systematics — so they agree in structure without coinciding, and
    # what is pinned here is structure and size, not identity. A transcription error in Table 2
    # would break all of it; the internal tests above would not notice.
    table = mass_table()
    systems = (
        (
            "235-U(n_th,f)",
            U235_POLARIZATION_FILE,
            neutron_induced_fission(Nuclide(92, 235), 2.53e-8),
            0.11,
        ),
        (
            "252-Cf(sf)",
            CF252_POLARIZATION_FILE,
            spontaneous_fission(Nuclide(98, 252)),
            0.30,
        ),
    )

    for (name, path, system, tolerance) in systems
        isfile(path) || continue
        evaluated = read_charge_polarization(path)
        systematics = WahlSystematics(table, system)

        # Mean absolute difference over the whole heavy branch. Held loosely: the two parameter
        # sets place the Zₚ = 50 crossing a mass unit or two apart, and that is the steepest
        # feature of the function, so a small horizontal offset costs a large pointwise
        # difference where the physics agrees. 252-Cf is the worse of the two because SL50 is
        # extrapolated 16 mass units past the data it was fitted on.
        differences = [
            abs(wahl_polarization(systematics, A) - ΔZ) for
            (A, ΔZ) in evaluated.polarization
        ]
        @test sum(differences) / length(differences) < tolerance

        # The narrow σ_Z band is the crossing made visible. Both the model and the table must
        # have one, and they must agree on where it is to within a few mass units.
        masses = sort(collect(keys(evaluated.width)))
        widths = [evaluated.width[A] for A in masses]
        threshold = minimum(widths) + 0.25 * (maximum(widths) - minimum(widths))
        tabulated_band = [A for (A, w) in zip(masses, widths) if w < threshold]
        @test !isempty(tabulated_band)
        @test abs(minimum(tabulated_band) - systematics.Bb) < 3.0
        @test abs(
            (minimum(tabulated_band) + maximum(tabulated_band)) / 2 -
            (systematics.Bb + systematics.B4) / 2,
        ) < 10.0

        # Bb is the mass at which ΔZ attains ΔZ_max, so the table's own maximum is where it
        # should be. This is what settles eq. (12): the blended A'_max lands within a mass unit
        # of it, and either estimate taken alone does worse. Zₚ there is recorded alongside,
        # since it is the quantity the interval is named for and the blend does not reach 50
        # exactly except where F₁ vanishes.
        polarizations = [evaluated.polarization[A] for A in masses]
        @test abs(masses[argmax(polarizations)] - systematics.Bb) < 1.5
        crossing =
            systematics.Bb * systematics.Z_F / systematics.A_F + systematics.ΔZmax
        @test 49.0 < crossing <= 50.0 + 1e-9
    end

    # Where A' = 140 falls inside the peak region the model is at its best determined, and the
    # agreement is then a real number rather than a bound. 240-Pu is the case: σ_Z(140) is
    # 0.59141 against a tabulated 0.59131.
    plutonium = joinpath(DATA, "charge_distribution", "Pu239.dat")
    if isfile(plutonium)
        evaluated = read_charge_polarization(plutonium)
        systematics =
            WahlSystematics(table, neutron_induced_fission(Nuclide(94, 239), 2.53e-8))
        @test 140 > systematics.B4
        @test wahl_width(systematics, 140) ≈ width_at(evaluated, 140) atol = 1e-3
    end
end
