# A run writes the charge distribution it was weighted by, beside the ones it could have taken
# instead, since the choice has a measurable effect on ν(A). The fixture takes the Zₚ model
# ("wahl") and the shipped mass table only, so the run needs no unshipped input.

@testset "a run records the charge distribution it used" begin
    mktempdir() do root
        data = joinpath(root, "data")
        mkpath(joinpath(data, "Cf252_sf"))
        cp(MASS_EXCESS_FILE, joinpath(data, "MASSES"))
        write(
            joinpath(data, "Cf252_sf", "R_T.csv"),
            "A_H,R_T,R_T_uncertainty\n126,1.0,0.0\n174,0.6,0.0\n",
        )
        write(
            joinpath(data, "Cf252_sf", "r_nu.csv"),
            "A_H,r_nu,r_nu_uncertainty\n126,0.5,0.0\n174,0.5,0.0\n",
        )
        write(
            joinpath(data, "Cf252_sf", "manifest_probe.toml"),
            """
            [system]
            label = "Cf252_sf"
            notation = "²⁵²Cf(sf)"
            target_A = 252
            target_Z = 98
            channel = "sf"
            reaction = "0,f"
            incident_energy_MeV = 0.0
            compound_A = 252
            compound_Z = 98

            [run]
            ordinate = "R_T"
            abscissa = ["A_H"]

            [[segmented_curve]]
            label = "systematic_trend"
            kind = "systematic_trend"
            temperature_ratio_file = "R_T.csv"
            multiplicity_ratio_pivots_file = "r_nu.csv"
            """,
        )
        write(
            joinpath(root, "Cf252_sf.toml"),
            """
            [system]
            target_A = 252
            target_Z = 98
            channel = "sf"
            incident_energy = 0.0

            [fragmentation]
            charges_per_mass = 3
            heavy_mass_min = 126
            heavy_mass_max = 140
            charge_distribution_file = "wahl"

            [energetics]
            mass_excess_file = "MASSES"
            kinetic_energy_min = 170.0
            kinetic_energy_max = 190.0
            kinetic_energy_step = 20.0

            [level_density]
            model = "BSFG"

            [evaporation]
            cross_section = "variable"
            maximum_temperature = 2.0

            [partition]
            method = "temperature_ratio"
            subdirectory = "Cf252_sf"
            curve_label = "systematic_trend"

            [spectrum]
            enabled = true
            energy_min = 0.2
            energy_max = 20.0
            points = 24

            [output]
            significant_digits = 8
            """,
        )

        configuration =
            load_configuration(joinpath(root, "Cf252_sf.toml"); data_directory = data)
        result = run_pipeline(configuration)
        destination = joinpath(root, "results")
        @test write_results(result, destination) == destination
        # A directory that already holds a run is refused, not overwritten.
        @test_throws ArgumentError write_results(result, destination)

        used = joinpath(destination, "charge_distribution_vs_A.dat")
        per_reaction =
            joinpath(destination, "charge_distribution_vs_A_wahl_per_reaction.dat")
        systematics = joinpath(destination, "charge_distribution_vs_A_wahl_systematics.dat")
        mean = joinpath(destination, "charge_distribution_vs_A_mean.dat")
        @test isfile(used)
        # ²⁵²Cf(sf) is one of the four reactions Table A of the 1988 evaluation covers, so the
        # first layer is written too.
        @test isfile(per_reaction)
        @test isfile(systematics)
        @test isfile(mean)

        # Each is the same three-column form an evaluation is supplied in, over exactly the
        # heavy masses the run swept, so the four are directly comparable.
        for path in (used, per_reaction, systematics, mean)
            table = read_charge_distribution(path)
            @test sort(collect(keys(table.ΔZ))) == collect(126:140)
            @test first(readlines(path)) == "A dZ sigma_Z"
        end

        # What was written is what the run was weighted by, reduced to effective Gaussians, to
        # the eight significant digits of the output.
        written = read_charge_distribution(used)
        weighted = effective_charge_distribution(configuration.charge_distribution, 126:140)
        @test maximum(
            abs(charge_polarization(written, A) - charge_polarization(weighted, A)) for
            A in 126:140
        ) <= 1e-7
        @test maximum(
            abs(charge_dispersion(written, A) - charge_dispersion(weighted, A)) for
            A in 126:140
        ) <= 1e-7

        # The fallback carries no mass dependence at all, which is the whole point of drawing it
        # beside the other two.
        fallback = read_charge_distribution(mean)
        @test all(charge_polarization(fallback, A) == -0.5 for A in 126:140)
        @test all(charge_dispersion(fallback, A) == 0.6 for A in 126:140)

        # And the systematics is the 2002 Wahl model, not a copy of either.
        model = read_charge_distribution(systematics)
        reference = effective_charge_distribution(
            WahlSystematics(configuration.mass_excess_table, configuration.system),
            126:140,
        )
        @test charge_polarization(model, 130) ≈ charge_polarization(reference, 130) atol =
            1e-7

        # And the per-reaction layer is neither of the other two: Table A's parameters are
        # narrower than the systematics estimates them to be for this reaction.
        fitted = read_charge_distribution(per_reaction)
        @test charge_dispersion(fitted, 130) < charge_dispersion(model, 130)
        @test charge_polarization(fitted, 130) != charge_polarization(model, 130)

        # The laboratory spectrum is written on the configured grid. It is normalised over the
        # half-line, so its trapezoidal area over [0.2, 20] MeV falls short of unity by the
        # low-energy tail, offset in part by the discretization error of the coarse grid.
        spectrum_file = joinpath(destination, "spectrum_vs_E.dat")
        @test isfile(spectrum_file)
        spectrum = read_delimited_table(spectrum_file, TableSpec([:E, :N]; skip = 1))
        E = column(spectrum, :E, Float64)
        N = column(spectrum, :N, Float64)
        @test length(E) == 24
        area = sum((E[i] - E[i - 1]) * (N[i] + N[i - 1]) / 2 for i in 2:length(E))
        @test 0.85 <= area <= 1.02
        @test isfinite(result.mean_neutron_lab_energy)
        @test 1.0 < result.mean_neutron_lab_energy < 3.0

        # The provenance record is the script's to write, not the library's.
        @test !isfile(joinpath(destination, "metadata.toml"))

        # What the library supplies of it names the resolved curve, so the partition input
        # is traceable, and is TOML-representable as it stands.
        metadata = run_metadata(result)
        @test Set(keys(metadata)) == Set(["system", "model", "result"])
        @test metadata["system"]["label"] == "Cf252_sf"
        @test metadata["result"]["configurations_solved"] == result.configurations_solved
        @test metadata["model"]["partition_curve"] == configuration.method.curve.source
        @test endswith(metadata["model"]["partition_curve"], "R_T.csv")
        parsed = TOML.parse(sprint(io -> TOML.print(io, metadata; sorted = true)))
        @test parsed["system"]["label"] == "Cf252_sf"
        @test isequal(parsed["result"]["mean_E_MeV"], metadata["result"]["mean_E_MeV"])
        @test parsed["model"]["heavy_masses"] == [126, 140]
    end
end

@testset "the input yield's mean kinetic energy" begin
    yields = FragmentYield(Dict((140, 54, 170.0) => 1.0, (96, 38, 180.0) => 3.0), 4.0)
    @test input_mean_kinetic_energy(yields) ≈ 177.5 rtol = RTOL
end
