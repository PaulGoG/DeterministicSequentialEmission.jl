# The fragment yield built from the staged ²⁵²Cf joint distribution, against the archived thesis
# run. The table has no citable origin, so this regression lives with the package that holds it.

guarded(
    CHARGE_DISTRIBUTION_AVAILABLE && YIELD_AVAILABLE,
    "fragment yields against the archived run",
) do
    @testset "fragment yields against the archived run" begin
        distribution = cf252_charge_distribution()
        domain = cf252_domain()
        measured = read_mass_energy_yield(CF252_YIELD_FILE)
        yields = build_fragment_yield(domain, distribution, measured)

        # The archived matrix, cell for cell.
        @test length(yields) == 48985
        @test sum(yields) ≈ 200 rtol = RTOL

        # Each nuclide takes the charge probability of its own mass. Summing a nuclide's
        # contributions instead, as the fragment census does, is wrong here and shows up at the
        # symmetric mass as a factor of 28.
        probabilities = fragment_charge_probabilities(domain, distribution)
        @test probabilities[Nuclide(50, 126)] ≈ 0.405781 atol = 1e-6
        @test probabilities[Nuclide(48, 126)] ≈ 0.014324 atol = 1e-6

        # Mirror symmetry away from the symmetric mass: the light fragment's Gaussian is the
        # heavy one's reflected, so the two agree without being computed the same way.
        @test probabilities[Nuclide(52, 134)] ≈ probabilities[Nuclide(46, 118)] rtol = RTOL

        for (dimension, file, column_names) in (
            (:mass, "Y_vs_A.dat", [:A, :Y]),
            (:charge, "Y_vs_Z.dat", [:Z, :Y]),
            (:kinetic_energy, "Y_vs_TKE.dat", [:TKE, :Y]),
        )
            reference = read_delimited_table(
                joinpath(REFERENCE, file),
                TableSpec(column_names; skip = 1),
            )
            arguments = column(reference, column_names[1], Float64)
            expected = column(reference, column_names[2], Float64)
            got = Dict(yield_over(yields, dimension))

            @test length(got) == length(reference)
            @test issubset(arguments, keys(got))
            # Reproduced in double precision, not to a tolerance chosen to pass.
            @test maximum(
                abs(get(got, x, NaN) - y) for (x, y) in zip(arguments, expected)
            ) < 1e-12

            # Every projection carries the whole distribution.
            @test sum(values(got)) ≈ 200 rtol = RTOL
        end

        # ⟨TKE⟩(A) against the reference.
        reference = read_delimited_table(
            joinpath(REFERENCE, "TKE_vs_A_H.dat"),
            TableSpec([:A_H, :TKE]; skip = 1),
        )
        masses = column(reference, :A_H, Int)
        expected = column(reference, :TKE, Float64)
        got = Dict(mean_kinetic_energy_by_mass(yields))
        @test issubset(masses, keys(got))
        @test maximum(abs(get(got, A, NaN) - y) for (A, y) in zip(masses, expected)) < 1e-10

        # The measurement is tabulated every 1 MeV and the archived run swept every 2, so half
        # the yield sits at kinetic energies the model never solved at. Pinned because it is a
        # property of the reference configuration that any average over it has to account for.
        report = coverage(yields, CF252_TKE)
        @test length(report.tabulated) == 101
        @test length(report.swept) == 51
        @test report.reached≈0.5 atol = 1e-4
        @test report.missed≈0.5 atol = 1e-4
    end
end
