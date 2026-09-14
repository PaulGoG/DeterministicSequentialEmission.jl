# A run writes the charge distribution it was weighted by, beside the two it could have taken
# instead. Which one a run used is a choice with a measurable effect on ν(A), and before this it
# was recorded only as a path in a configuration.

CHARGE_DISTRIBUTION_AVAILABLE &&
    @testset "a run records the charge distribution it used" begin
        mktempdir() do root
            data = joinpath(root, "data")
            mkpath(joinpath(data, "Cf252_sf"))
            cp(MASS_EXCESS_FILE, joinpath(data, "MASSES"))
            cp(
                CF252_CHARGE_DISTRIBUTION_FILE,
                joinpath(data, "Cf252_sf", "charge_distribution_vs_A.dat"),
            )
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
                charge_distribution_file = "Cf252_sf/charge_distribution_vs_A.dat"

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

                [output]
                significant_digits = 8
                """,
            )

            configuration =
                load_configuration(joinpath(root, "Cf252_sf.toml"); data_directory = data)
            result = run_pipeline(configuration; write = false)
            destination = write_results(result; directory = joinpath(root, "results"))

            used = joinpath(destination, "charge_distribution_vs_A.dat")
            systematics =
                joinpath(destination, "charge_distribution_vs_A_wahl_systematics.dat")
            mean = joinpath(destination, "charge_distribution_vs_A_mean.dat")
            @test isfile(used)
            @test isfile(systematics)
            @test isfile(mean)

            # Each is the same three-column form an evaluation is supplied in, over exactly the
            # heavy masses the run swept, so the three are directly comparable.
            for path in (used, systematics, mean)
                table = read_charge_distribution(path)
                @test sort(collect(keys(table.ΔZ))) == collect(126:140)
                @test first(readlines(path)) == "A dZ sigma_Z"
            end

            # What was written is what the run was weighted by, not a fresh read of the file.
            written = read_charge_distribution(used)
            for A in 126:140
                @test charge_polarization(written, A) ≈
                      charge_polarization(configuration.charge_distribution, A)
                @test charge_dispersion(written, A) ≈
                      charge_dispersion(configuration.charge_distribution, A)
            end

            # The fallback carries no mass dependence at all, which is the whole point of drawing it
            # beside the other two.
            fallback = read_charge_distribution(mean)
            @test all(charge_polarization(fallback, A) == -0.5 for A in 126:140)
            @test all(charge_dispersion(fallback, A) == 0.6 for A in 126:140)

            # And the systematics is the Wahl model, not a copy of either.
            model = read_charge_distribution(systematics)
            reference =
                WahlSystematics(configuration.mass_excess_table, configuration.system)
            @test charge_polarization(model, 130) ≈ wahl_polarization(reference, 130) atol =
                1e-7
            @test charge_dispersion(model, 130) ≈ wahl_dispersion(reference, 130) atol =
                1e-7

            # The run metadata names the resolved curve, so the partition input is traceable too.
            metadata = read(joinpath(destination, "metadata.toml"), String)
            @test occursin("partition_curve", metadata)
            @test occursin("R_T.csv", metadata)
        end
    end
