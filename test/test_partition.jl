@testset "jointed segments" begin
    segments = JointedSegments([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)])
    @test length(segments) == 4

    # Endpoints are reproduced exactly, and interpolation is linear between them.
    @test segments(126) ≈ 0.5
    @test segments(130) ≈ 0.14
    @test segments(150) ≈ 0.607
    @test segments(174) ≈ 0.865
    @test segments(140) ≈ 0.14 + (0.607 - 0.14) * (140 - 130) / (150 - 130)

    # Undefined above the last point; the published parameterizations end at A_H_max.
    @test segments(175) === nothing

    # Below the first point the first segment is extended rather than clamped.
    @test segments(125) ≈ 0.5 + (0.14 - 0.5) * (125 - 126) / (130 - 126)

    # A segment crossing zero leaves the function undefined there rather than negative.
    crossing = JointedSegments([(120, 1.0), (140, -1.0)])
    @test crossing(130) === nothing
    @test crossing(125) > 0

    @test_throws ArgumentError JointedSegments([(126, 0.5)])
    @test_throws ArgumentError JointedSegments([(150, 0.5), (126, 0.1)])
    @test_throws ArgumentError JointedSegments([(126, 0.5), (126, 0.1)])
end

@testset "reading a tabulated parameterization" begin
    # The layout FissionTemperatureRatio.jl writes: a header line, then A_H, value and an
    # uncertainty column that the partition discards because it carries none downstream.
    mktempdir() do dir
        path = joinpath(dir, "R_T_parameterized.csv")
        write(
            path,
            """
            A_H,value,uncertainty
            126,1.0,0.0
            127,1.205662,0.007435
            128,1.471892,0.02157
            129,1.932797,0.061908
            """,
        )
        segments = read_segments(path)
        @test length(segments) == 4

        # Tabulated at every integer mass, so evaluation reproduces each value exactly rather
        # than interpolating across it.
        @test segments(126) == 1.0
        @test segments(127) == 1.205662
        @test segments(128) == 1.471892
        @test segments(129) == 1.932797

        # The range is the fit's, not the configured domain's: a set that stops short leaves
        # the partition undefined above its last point instead of extrapolating.
        @test segments(130) === nothing

        # The uncertainty column is optional, as it is for every other reader here.
        bare = joinpath(dir, "bare.csv")
        write(bare, "A_H,value\n126,1.0\n127,1.2\n")
        @test read_segments(bare)(127) == 1.2
    end
end

@testset "excitation energy ratio against the archived run" begin
    table = mass_table()

    # The parameterization the archived run used, transcribed from its own provenance header.
    method = ExcitationEnergyRatio(
        JointedSegments([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )

    reference = read_table(
        joinpath(REFERENCE, "CF52_excitation_PARAM.dat"),
        TableSpec([:A_H, :Z_H, :TKE, :E_excit_L, :E_excit_H]; skip = 1),
    )
    masses = column(reference, :A_H, Int)
    charges = column(reference, :Z_H, Int)
    energies = column(reference, :TKE, Float64)
    expected_light = column(reference, :E_excit_L, Float64)
    expected_heavy = column(reference, :E_excit_H, Float64)
    @test length(reference) == 1200

    worst = 0.0
    checked = 0
    for row in eachindex(masses)
        fragmentation = fragmentation_at(masses[row], charges[row])
        balance = energetics(table, CF252, fragmentation, energies[row])
        @test balance !== nothing
        balance === nothing && continue
        split = partition(method, fragmentation, value(balance.TXE))
        @test split !== nothing
        split === nothing && continue
        worst = max(worst, abs(split.light - expected_light[row]))
        worst = max(worst, abs(split.heavy - expected_heavy[row]))
        checked += 1
    end
    @test checked == 1200
    # The reference is rounded to eight decimals on write.
    @test worst < 1e-8
end

POLARIZATION_AVAILABLE && @testset "partition conserves TXE" begin
    table = mass_table()
    domain = cf252_domain()
    method = ExcitationEnergyRatio(
        JointedSegments([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )

    for fragmentation in domain.entries[1:20:end]
        balance = energetics(table, CF252, fragmentation, 180.0)
        balance === nothing && continue
        split = partition(method, fragmentation, value(balance.TXE))
        split === nothing && continue
        @test sum(split) ≈ value(balance.TXE)
        @test split.light >= 0
        @test split.heavy >= 0
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
        TemperatureRatio(JointedSegments([(126, 1.2), (174, 1.2)]), model),
        ExcitationEnergyRatio(JointedSegments([(126, 0.42), (174, 0.42)])),
    ]
    for method in methods
        split = partition(method, identical, TXE)
        @test split.light ≈ split.heavy
        @test sum(split) ≈ TXE
    end

    # Scission modelling needs no guard: identical fragments give a_L/a_H = 1 and equal
    # extra-deformation energies, so it divides evenly on its own arithmetic.
    scission = ScissionModelling(Dict(Nuclide(49, 126) => 3.0), model)
    split = partition(scission, identical, TXE)
    @test split.light ≈ split.heavy
    @test sum(split) ≈ TXE

    # A mass-symmetric split with unequal charges is a different matter: the fragments are
    # distinct nuclides, a_L ≠ a_H, and the prescription must be applied as given.
    unequal = fragmentation_at(126, 50)
    @test unequal.light == Nuclide(48, 126)
    split = partition(
        TemperatureRatio(JointedSegments([(126, 1.2), (174, 1.2)]), model),
        unequal,
        TXE,
    )
    @test !(split.light ≈ split.heavy)
    @test sum(split) ≈ TXE
end

@testset "temperature ratio against the archived run" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # The second archived run, whose R_T pivots come from its own provenance header. The
    # instructions file shipped with the prototype lists a different tail — (167, 0.69),
    # (174, 0.6) — which is not what produced the published figures.
    method = TemperatureRatio(
        JointedSegments([
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

    reference = read_table(
        joinpath(REFERENCE, "CF52_excitation_RT.dat"),
        TableSpec([:A_H, :Z_H, :TKE, :E_excit_L, :E_excit_H]; skip = 1),
    )
    masses = column(reference, :A_H, Int)
    charges = column(reference, :Z_H, Int)
    energies = column(reference, :TKE, Float64)
    expected_light = column(reference, :E_excit_L, Float64)
    expected_heavy = column(reference, :E_excit_H, Float64)
    @test length(reference) == 1200

    worst = 0.0
    checked = 0
    for row in eachindex(masses)
        fragmentation = fragmentation_at(masses[row], charges[row])
        balance = energetics(table, CF252, fragmentation, energies[row])
        @test balance !== nothing
        balance === nothing && continue
        split = partition(method, fragmentation, value(balance.TXE))
        @test split !== nothing
        split === nothing && continue
        worst = max(worst, abs(split.light - expected_light[row]))
        worst = max(worst, abs(split.heavy - expected_heavy[row]))
        checked += 1
    end
    @test checked == 1200
    # The reference is rounded to eight decimals on write, so 5e-9 is exact agreement.
    @test worst < 1e-8
end

@testset "temperature ratio partition" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    method = TemperatureRatio(
        JointedSegments([
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
    split = partition(method, fragmentation, value(balance.TXE))
    @test split !== nothing
    @test sum(split) ≈ value(balance.TXE)

    # R_T > 1 means the light fragment is hotter, so with comparable level density
    # parameters it takes the larger share of TXE.
    R_T = method.segments(134)
    @test R_T > 1
    @test split.light > split.heavy

    # At R_T = 1 the split is fixed by the level density parameters alone, E*_L/E*_H = a_L/a_H.
    unit = TemperatureRatio(JointedSegments([(126, 1.0), (174, 1.0)]), model)
    split_unit = partition(unit, fragmentation, value(balance.TXE))
    a_light = level_density_parameter(model, fragmentation.light)
    a_heavy = level_density_parameter(model, fragmentation.heavy)
    @test split_unit.light / split_unit.heavy ≈ a_light / a_heavy
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
    @test split.light / split.heavy ≈ a_light / a_heavy
    @test sum(split) ≈ 30.0

    # Extra deformation energy is returned to its own fragment on top of the thermal share.
    deformed = ScissionModelling(
        Dict(fragmentation.light => 2.0, fragmentation.heavy => 5.0),
        model,
    )
    split = partition(deformed, fragmentation, 30.0)
    @test sum(split) ≈ 30.0
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
