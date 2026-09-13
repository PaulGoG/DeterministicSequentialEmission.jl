@testset "experimental yield input" begin
    mktempdir() do dir
        path = joinpath(dir, "Y.dat")
        write(
            path,
            """
            A TKE Value σ
            140 170.0 1.5 0.1
            140 171.0 2.5 0.2
            141 170.0 3.0 NaN
            141 171.0 NaN NaN
            """,
        )
        measured = read_experimental_yield(path)

        # A cell with no finite yield is dropped: an unmeasured (A, TKE) is not a yield of zero.
        @test length(measured) == 3
        @test measured.values[(140, 171.0)] == 2.5
        @test !haskey(measured.values, (141, 171.0))
        @test measured.masses == [140, 141]
        @test measured.energies == [170.0, 171.0]
        @test occursin("2 masses", sprint(show, measured))

        # The uncertainty column is optional, as everywhere else.
        bare = joinpath(dir, "bare.dat")
        write(bare, "A TKE Value\n140 170.0 1.5\n")
        @test length(read_experimental_yield(bare)) == 1

        repeated = joinpath(dir, "repeated.dat")
        write(repeated, "A TKE Value\n140 170.0 1.5\n140 170.0 2.0\n")
        @test_throws ArgumentError read_experimental_yield(repeated)

        empty = joinpath(dir, "empty.dat")
        write(empty, "A TKE Value\n140 170.0 NaN\n")
        @test_throws ArgumentError read_experimental_yield(empty)
    end
end

@testset "coverage reports what a sweep reaches" begin
    yields = FragmentYields(
        Dict((140, 54, 170.0) => 3.0, (140, 54, 171.0) => 1.0, (140, 54, 172.0) => 4.0),
        8.0,
    )
    report = coverage(yields, 170.0:2.0:172.0)
    @test report.reached ≈ 7 / 8
    @test report.missed ≈ 1 / 8
    @test report.tabulated == [170.0, 171.0, 172.0]
    @test report.swept == [170.0, 172.0]
end

POLARIZATION_AVAILABLE &&
    YIELD_AVAILABLE &&
    @testset "fragment yields against the archived run" begin
        polarization = cf252_polarization()
        domain = cf252_domain()
        measured = read_experimental_yield(CF252_YIELD_FILE)
        yields = build_fragment_yields(domain, polarization, measured)

        # The archived matrix, cell for cell.
        @test length(yields) == 48985
        @test sum(yields) ≈ 200

        # Each nuclide takes the charge probability of its own mass. Summing a nuclide's
        # contributions instead, as the fragment census does, is wrong here and shows up at the
        # symmetric mass as a factor of 28.
        probabilities = fragment_charge_probabilities(domain, polarization)
        @test probabilities[Nuclide(50, 126)] ≈ 0.405781 atol = 1e-6
        @test probabilities[Nuclide(48, 126)] ≈ 0.014324 atol = 1e-6

        # Mirror symmetry away from the symmetric mass: the light fragment's Gaussian is the
        # heavy one's reflected, so the two agree without being computed the same way.
        @test probabilities[Nuclide(52, 134)] ≈ probabilities[Nuclide(46, 118)]

        for (dimension, file, column_names) in (
            (:mass, "CF52_Y_A.dat", [:A, :Y]),
            (:charge, "CF52_Y_Z.dat", [:Z, :Y]),
            (:kinetic_energy, "CF52_Y_TKE.dat", [:TKE, :Y]),
        )
            reference =
                read_table(joinpath(REFERENCE, file), TableSpec(column_names; skip = 1))
            arguments = column(reference, column_names[1], Float64)
            expected = column(reference, column_names[2], Float64)
            got = Dict(yield_over(yields, dimension))

            @test length(got) == length(reference)
            worst = 0.0
            for row in eachindex(arguments)
                @test haskey(got, arguments[row])
                haskey(got, arguments[row]) || continue
                worst = max(worst, abs(got[arguments[row]] - expected[row]))
            end
            # Reproduced in double precision, not to a tolerance chosen to pass.
            @test worst < 1e-12

            # Every projection carries the whole distribution.
            @test sum(values(got)) ≈ 200
        end

        # ⟨TKE⟩(A) against the archive.
        reference = read_table(
            joinpath(REFERENCE, "CF52_TKE_AH.dat"),
            TableSpec([:A_H, :TKE]; skip = 1),
        )
        masses = column(reference, :A_H, Int)
        expected = column(reference, :TKE, Float64)
        got = Dict(mean_kinetic_energy_by_mass(yields))
        worst = 0.0
        for row in eachindex(masses)
            @test haskey(got, masses[row])
            haskey(got, masses[row]) || continue
            worst = max(worst, abs(got[masses[row]] - expected[row]))
        end
        @test worst < 1e-10

        # The measurement is tabulated every 1 MeV and the archived run swept every 2, so half
        # the yield sits at kinetic energies the model never solved at. Pinned because it is a
        # property of the reference configuration that any average over it has to account for.
        report = coverage(yields, CF252_TKE)
        @test length(report.tabulated) == 101
        @test length(report.swept) == 51
        @test report.reached≈0.5 atol = 1e-4
        @test report.missed≈0.5 atol = 1e-4
    end

@testset "yield construction validation" begin
    @test_throws ArgumentError yield_over(
        FragmentYields(Dict((1, 1, 1.0) => 1.0), 1.0),
        :temperature,
    )
end
