@testset "parameterization values are checked against their meaning" begin
    # E*_H/TXE is a fraction of the total excitation energy.
    @test_throws ArgumentError ExcitationEnergyRatio(
        SegmentedCurve([(126, 0.5), (150, 1.4)]),
    )
    @test_throws ArgumentError ExcitationEnergyRatio(
        SegmentedCurve([(126, -0.1), (150, 0.6)]),
    )
    @test ExcitationEnergyRatio(SegmentedCurve([(126, 0.5), (150, 0.6)])) isa
          ExcitationEnergyRatio

    # R_T is a ratio of temperatures, so positive, but it may exceed one.
    model = BackShiftedFermiGas(mass_table())
    @test_throws ArgumentError TemperatureRatio(
        SegmentedCurve([(126, 1.0), (150, -0.3)]),
        model,
    )
    @test TemperatureRatio(SegmentedCurve([(126, 1.0), (150, 1.7)]), model) isa
          TemperatureRatio
end

@testset "segment coverage against a sweep" begin
    # A parameterization fitted over part of the domain leaves the rest unpartitioned, and that
    # exclusion is indistinguishable downstream from one the physics made.
    short = ExcitationEnergyRatio(SegmentedCurve([(126, 0.5), (150, 0.6)]))
    report = segment_coverage(short, 126:174)
    @test report.uncovered == collect(151:174)
    @test report.covered ≈ 25 / 49 rtol = RTOL

    # A fit that starts above the first heavy mass leaves the low end uncovered, as one that
    # stops short leaves the high end.
    late = ExcitationEnergyRatio(SegmentedCurve([(130, 0.5), (174, 0.6)]))
    report = segment_coverage(late, 126:174)
    @test report.uncovered == collect(126:129)
    @test report.covered ≈ 45 / 49 rtol = RTOL

    whole = ExcitationEnergyRatio(
        SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )
    full = segment_coverage(whole, 126:174)
    @test isempty(full.uncovered)
    @test full.covered == 1.0

    # Modelling at scission has no curve; the extra-deformation table sets the reach instead.
    deformation =
        Dict(Nuclide(50, 130) => 1.2, Nuclide(51, 130) => 1.0, Nuclide(52, 132) => 0.8)
    scission = ScissionModelling(deformation, BackShiftedFermiGas(mass_table()))
    reached = segment_coverage(scission, 130:133)
    @test reached.uncovered == [131, 133]
    @test reached.covered == 0.5
end

@testset "excitation energy ratio against the archived run" begin
    table = mass_table()

    # The parameterization the archived run used, transcribed from its own provenance header.
    method = ExcitationEnergyRatio(
        SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )

    reference = read_delimited_table(
        joinpath(REFERENCE, "E_star_vs_A_H_Z_H_TKE_excitation_ratio.dat"),
        TableSpec([:A_H, :Z_H, :TKE, :E_star_L, :E_star_H]; skip = 1),
    )
    masses = column(reference, :A_H, Int)
    charges = column(reference, :Z_H, Int)
    energies = column(reference, :TKE, Float64)
    expected_light = column(reference, :E_star_L, Float64)
    expected_heavy = column(reference, :E_star_H, Float64)
    @test length(reference) == 1200

    fragmentations = [fragmentation_at(A, Z) for (A, Z) in zip(masses, charges)]
    balances = [
        energetics(table, CF252, fragmentation, TKE) for
        (fragmentation, TKE) in zip(fragmentations, energies)
    ]
    @test count(isnothing, balances) == 0
    splits = [
        partition(method, fragmentation, value(balance.TXE)) for
        (fragmentation, balance) in zip(fragmentations, balances)
    ]
    @test count(isnothing, splits) == 0
    # The reference is rounded to eight decimals on write.
    @test maximum(abs(s.light - E) for (s, E) in zip(splits, expected_light)) < 1e-8
    @test maximum(abs(s.heavy - E) for (s, E) in zip(splits, expected_heavy)) < 1e-8
end

guarded(CHARGE_DISTRIBUTION_AVAILABLE, "the partition leaves no fragment negative") do
    @testset "the partition leaves no fragment negative" begin
        table = mass_table()
        domain = cf252_domain()
        method = ExcitationEnergyRatio(
            SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
        )

        splits = ExcitationEnergies[]
        for fragmentation in domain.entries[1:20:end]
            balance = energetics(table, CF252, fragmentation, 180.0)
            balance === nothing && continue
            split = partition(method, fragmentation, value(balance.TXE))
            split === nothing && continue
            push!(splits, split)
        end
        @test all(split.light >= 0 for split in splits)
        @test all(split.heavy >= 0 for split in splits)
    end
end

@testset "a split into two copies of one nuclide divides evenly" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # A_H = A_L = 126 with Z_H = Z_L = 49: the two fragments are the same nuclide, so they have
    # the same level density parameter and the same temperature identically. E*_L = E*_H is a
    # symmetry, and no parameterization of A_H can be relied on to produce it.
    identical = fragmentation_at(126, 49)
    @test identical.heavy == identical.light == Nuclide(49, 126)
    balance = energetics(table, CF252, identical, 160.0)
    TXE = value(balance.TXE)

    # Values chosen to be wrong at symmetry on purpose: R_T = 1.2 is the literature value for
    # 252-Cf and gives E*_H/TXE = 0.41 here, and 0.42 is a plausible ratio parameterization.
    # Both archived thesis runs happen to pin their first segment at R_T = 1 and 0.5, so
    # neither exercises this.
    methods = [
        TemperatureRatio(SegmentedCurve([(126, 1.2), (174, 1.2)]), model),
        ExcitationEnergyRatio(SegmentedCurve([(126, 0.42), (174, 0.42)])),
    ]
    for method in methods
        split = partition(method, identical, TXE)
        @test split.light ≈ split.heavy rtol = RTOL
    end

    # Scission modelling needs no guard: identical fragments give a_L/a_H = 1 and equal
    # extra-deformation energies, so it divides evenly on its own arithmetic.
    scission = ScissionModelling(Dict(Nuclide(49, 126) => 3.0), model)
    split = partition(scission, identical, TXE)
    @test split.light ≈ split.heavy rtol = RTOL
    @test sum(split) ≈ TXE rtol = RTOL

    # A mass-symmetric split with unequal charges is a different matter: the fragments are
    # distinct nuclides, a_L ≠ a_H, and the partition method must be applied as given.
    unequal = fragmentation_at(126, 50)
    @test unequal.light == Nuclide(48, 126)
    split = partition(
        TemperatureRatio(SegmentedCurve([(126, 1.2), (174, 1.2)]), model),
        unequal,
        TXE,
    )
    @test !isapprox(split.light, split.heavy; rtol = RTOL)
end

@testset "temperature ratio against the archived run" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # The second archived run, whose R_T pivots come from its own provenance header. The
    # instructions file shipped with the prototype lists a different tail — (167, 0.69),
    # (174, 0.6) — which is not what produced the published figures.
    method = TemperatureRatio(
        SegmentedCurve([
            (126, 1.0),
            (130, 1.734),
            (134, 1.28),
            (138, 1.28),
            (150, 0.9),
            (167, 0.68),
            (174, 0.59),
        ]),
        model,
    )

    reference = read_delimited_table(
        joinpath(REFERENCE, "E_star_vs_A_H_Z_H_TKE_temperature_ratio.dat"),
        TableSpec([:A_H, :Z_H, :TKE, :E_star_L, :E_star_H]; skip = 1),
    )
    masses = column(reference, :A_H, Int)
    charges = column(reference, :Z_H, Int)
    energies = column(reference, :TKE, Float64)
    expected_light = column(reference, :E_star_L, Float64)
    expected_heavy = column(reference, :E_star_H, Float64)
    @test length(reference) == 1200

    fragmentations = [fragmentation_at(A, Z) for (A, Z) in zip(masses, charges)]
    balances = [
        energetics(table, CF252, fragmentation, TKE) for
        (fragmentation, TKE) in zip(fragmentations, energies)
    ]
    @test count(isnothing, balances) == 0
    splits = [
        partition(method, fragmentation, value(balance.TXE)) for
        (fragmentation, balance) in zip(fragmentations, balances)
    ]
    @test count(isnothing, splits) == 0
    # The reference is rounded to eight decimals on write, so 5e-9 is exact agreement.
    @test maximum(abs(s.light - E) for (s, E) in zip(splits, expected_light)) < 1e-8
    @test maximum(abs(s.heavy - E) for (s, E) in zip(splits, expected_heavy)) < 1e-8
end

@testset "temperature ratio partition" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    method = TemperatureRatio(
        SegmentedCurve([
            (126, 1.0),
            (130, 1.734),
            (134, 1.28),
            (138, 1.28),
            (150, 0.9),
            (167, 0.68),
            (174, 0.59),
        ]),
        model,
    )

    fragmentation = fragmentation_at(134, 52)
    balance = energetics(table, CF252, fragmentation, 180.0)
    split = @inferred Union{ExcitationEnergies, Nothing} partition(
        method,
        fragmentation,
        value(balance.TXE),
    )
    @test split !== nothing

    # R_T > 1 means the light fragment is hotter, so with comparable level density
    # parameters it takes the larger share of TXE.
    R_T = method.curve(134)
    @test R_T > 1
    @test split.light > split.heavy

    # At R_T = 1 the split is fixed by the level density parameters alone, E*_L/E*_H = a_L/a_H.
    unit = TemperatureRatio(SegmentedCurve([(126, 1.0), (174, 1.0)]), model)
    split_unit = partition(unit, fragmentation, value(balance.TXE))
    a_light = level_density_parameter(model, fragmentation.light)
    a_heavy = level_density_parameter(model, fragmentation.heavy)
    @test split_unit.light / split_unit.heavy ≈ a_light / a_heavy rtol = RTOL
end

@testset "scission modelling partition" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    fragmentation = fragmentation_at(134, 52)

    # With no extra deformation the split reduces to thermal equilibrium at scission,
    # E*_L/E*_H = a_L/a_H, the same limit the temperature ratio reaches at R_T = 1.
    bare = ScissionModelling(
        Dict(fragmentation.light => 0.0, fragmentation.heavy => 0.0),
        model,
    )
    split = partition(bare, fragmentation, 30.0)
    a_light = level_density_parameter(model, fragmentation.light)
    a_heavy = level_density_parameter(model, fragmentation.heavy)
    @test split.light / split.heavy ≈ a_light / a_heavy rtol = RTOL
    @test sum(split) ≈ 30.0 rtol = RTOL

    # Extra deformation energy is returned to its own fragment on top of the thermal share.
    deformed = ScissionModelling(
        Dict(fragmentation.light => 2.0, fragmentation.heavy => 5.0),
        model,
    )
    split = partition(deformed, fragmentation, 30.0)
    @test sum(split) ≈ 30.0 rtol = RTOL
    @test split.heavy > partition(bare, fragmentation, 30.0).heavy

    # Not evaluable where a deformation energy is missing, or where the deformation exceeds
    # the energy available to share.
    @test partition(
        ScissionModelling(Dict{Nuclide, Float64}(), model),
        fragmentation,
        30.0,
    ) === nothing
    @test partition(deformed, fragmentation, 3.0) === nothing
end
