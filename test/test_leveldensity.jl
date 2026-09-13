@testset "shell correction file layout" begin
    corrections = read_shell_corrections(joinpath(INPUT, "Aux files", "SZSN.GC"))

    # Pinned to the file's own header, `n S(N) S(Z)`. The prototype declares the two columns
    # the other way round, and because S(Z) is looked up at the proton number while S(N) is
    # looked up at the neutron number, the swap is not absorbed by their sum.
    header = split(strip(first(readlines(joinpath(INPUT, "Aux files", "SZSN.GC")))))
    @test header[2] == "S(N)"
    @test header[3] == "S(Z)"
    @test corrections.neutron[11] ≈ 6.80
    @test corrections.proton[11] ≈ -2.91
end

@testset "back-shifted Fermi gas against the archived run" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    reference = read_table(
        joinpath(REFERENCE, "CF52_level_density.dat"),
        TableSpec([:A_H, :Z_H, :a_L, :a_H]; skip = 1),
    )
    heavy_masses = column(reference, :A_H, Int)
    heavy_charges = column(reference, :Z_H, Int)
    a_light = column(reference, :a_L, Float64)
    a_heavy = column(reference, :a_H, Float64)
    @test length(reference) == 244

    worst = 0.0
    for row in eachindex(heavy_masses)
        heavy = Nuclide(heavy_charges[row], heavy_masses[row])
        light = complementary_fragment(CF252, heavy)
        computed_heavy = level_density_parameter(model, heavy)
        computed_light = level_density_parameter(model, light)
        @test computed_heavy !== nothing
        @test computed_light !== nothing
        (computed_heavy === nothing || computed_light === nothing) && continue
        worst = max(worst, abs(computed_heavy - a_heavy[row]))
        worst = max(worst, abs(computed_light - a_light[row]))
    end
    # The reference is rounded to eight decimals on write.
    @test worst < 1e-8
end

@testset "BSFG pinned to von Egidy & Bucurescu" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # Eq. (19) of the 2009 paper, coefficients and all. Reconstructing a from S′ by hand must
    # reproduce what the implementation returns, which pins the additive form: the variant
    # p₁(1 + p₂ S′)A^p₃ would differ by a factor p₁ on the S′ term.
    for nuclide in (Nuclide(52, 134), Nuclide(42, 104), Nuclide(50, 132))
        S = shell_correction(model, nuclide)
        @test S !== nothing
        @test level_density_parameter(model, nuclide) ≈
              (0.199 + 0.0096 * S) * nuclide.A^0.869
        @test !isapprox(
            level_density_parameter(model, nuclide),
            0.199 * (1 + 0.0096 * S) * nuclide.A^0.869,
        )
    end

    # Eq. (12) of the 2009 paper: the pairing term enters as 0.5 P_a′ with
    # P_a′ = ½[M(A+2,Z+1) − 2M(A,Z) + M(A−2,Z−1)], i.e. a quarter of the second difference,
    # and with no (−1)^Z alternation — that factor distinguishes P_a′ from the P_a of the
    # 2005 paper and of Audi's tables, and it is P_a′ that enters eq. (19).
    nuclide = Nuclide(52, 134)
    Δ = value(mass_excess(table, nuclide))
    Δ_up = value(mass_excess(table, Nuclide(53, 136)))
    Δ_down = value(mass_excess(table, Nuclide(51, 132)))
    pairing = (Δ_up - 2Δ + Δ_down) / 4
    δW₀ =
        liquid_drop_energy(model, nuclide) - (
            nuclide.Z * value(mass_excess(table, Nuclide(1, 1))) +
            neutron_number(nuclide) * value(mass_excess(table, NEUTRON)) - Δ
        )
    @test shell_correction(model, nuclide) ≈ δW₀ + pairing

    # S = M_exp − M_LD (2005, eq. 7), so a nucleus more bound than the liquid drop has S < 0.
    # Shell closures are exactly those nuclei, and the suppressed a is what puts the sawtooth
    # into E*(A) and hence into ν(A).
    @test shell_correction(model, Nuclide(50, 132)) < -8.0      # doubly magic
    @test shell_correction(model, Nuclide(42, 104)) > 0.0       # mid-shell
    @test abs(shell_correction(model, Nuclide(50, 120))) < 1.0  # magic Z, mid-shell N

    # a/A ≈ 1/8 away from closures, and well below it at one.
    @test level_density_parameter(model, Nuclide(42, 104)) / 104 > 0.11
    @test level_density_parameter(model, Nuclide(50, 132)) / 132 < 0.06
end

@testset "level density parameter behaviour" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)

    # Roughly A/8 MeV⁻¹ across the fragment range, the familiar magnitude for these nuclei.
    for nuclide in (Nuclide(52, 134), Nuclide(46, 118), Nuclide(38, 96))
        a = level_density_parameter(model, nuclide)
        @test a !== nothing
        @test 0.05 * nuclide.A < a < 0.2 * nuclide.A
    end

    # The N = 82, Z = 50 closure suppresses the parameter relative to its neighbours, which is
    # the origin of the sawtooth in E*(A) and hence in ν(A).
    magic = level_density_parameter(model, Nuclide(50, 132))
    neighbour = level_density_parameter(model, Nuclide(50, 138))
    @test magic < neighbour

    # Unavailable rather than an error where the neighbouring masses are not tabulated.
    @test level_density_parameter(model, Nuclide(1, 1)) === nothing

    corrections = read_shell_corrections(joinpath(INPUT, "Aux files", "SZSN.GC"))
    gilbert = GilbertCameron(corrections)
    a_gc = level_density_parameter(gilbert, Nuclide(52, 134))
    @test a_gc !== nothing
    @test 0.05 * 134 < a_gc < 0.2 * 134

    # Outside the tabulated nucleon numbers the systematics simply has nothing to say.
    @test level_density_parameter(gilbert, Nuclide(2, 4)) === nothing
end
