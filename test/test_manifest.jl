# The manifest is what makes the quantity a matter of record rather than of staging. These
# exercise every refusal, because a manifest that is read leniently is worth no more than the
# hand-staged path it replaced.

"Write a manifest with `replacements` applied to the fixture, and the two tables it names."
function _manifest(directory::AbstractString, replacements::Pair{String, String}...)
    text = """
    [system]
    label = "U235_nth"
    notation = "²³⁵U(nth,f)"
    target_A = 235
    target_Z = 92
    channel = "nth"
    reaction = "n,f"
    incident_energy_MeV = 2.53e-8
    compound_A = 236
    compound_Z = 92

    [run]
    ordinate = "R_T"
    abscissa = ["A_H"]
    columns = ["A_H", "R_T", "R_T_uncertainty"]

    [[segmented_curve]]
    label = "Goeoek_2018"
    kind = "dataset"
    temperature_ratio_file = "R_T_vs_A_H_segmented_Goeoek_2018_run.csv"
    multiplicity_ratio_pivots_file = "r_nu_vs_A_H_pivots_Goeoek_2018_run.csv"

    [[segmented_curve]]
    label = "systematic_trend"
    kind = "systematic_trend"
    temperature_ratio_file = "R_T_vs_A_H_segmented_systematic_trend_run.csv"
    multiplicity_ratio_pivots_file = "r_nu_vs_A_H_pivots_systematic_trend_run.csv"
    """
    for (from, to) in replacements
        occursin(from, text) || error("the fixture has no $(repr(from)) to replace")
        text = replace(text, from => to)
    end
    # R_T rises through unity; r_nu sits at one half. Both are plausible floats in the same
    # column, which is the whole reason the manifest exists.
    for label in ("Goeoek_2018", "systematic_trend")
        write(
            joinpath(directory, "R_T_vs_A_H_segmented_$(label)_run.csv"),
            "A_H,R_T,R_T_uncertainty\n118,1.0,0.0\n140,1.2,0.01\n160,1.4,0.02\n",
        )
        write(
            joinpath(directory, "r_nu_vs_A_H_pivots_$(label)_run.csv"),
            "A_H,r_nu,r_nu_uncertainty\n118,0.5,0.0\n140,0.5,0.01\n160,0.5,0.02\n",
        )
    end
    path = joinpath(directory, "manifest_run.toml")
    write(path, text)
    return path
end

@testset "the partition takes its curve through the manifest" begin
    mktempdir() do root
        system = joinpath(root, "U235_nth")
        mkpath(system)
        _manifest(system)
        isfile(MASS_EXCESS_FILE) || return nothing
        cp(MASS_EXCESS_FILE, joinpath(root, "MASSES"))

        text = """
        [system]
        target_A = 235
        target_Z = 92
        channel = "nth"
        incident_energy = 2.53e-8

        [fragmentation]
        charges_per_mass = 5
        heavy_mass_min = 118
        heavy_mass_max = 160
        charge_distribution_file = "mean"

        [energetics]
        mass_excess_file = "MASSES"
        kinetic_energy_min = 150.0
        kinetic_energy_max = 190.0
        kinetic_energy_step = 20.0

        [level_density]
        model = "BSFG"

        [evaporation]
        cross_section = "variable"
        maximum_temperature = 2.0

        [partition]
        method = "temperature_ratio"
        subdirectory = "U235_nth"
        curve_label = "systematic_trend"

        [output]
        significant_digits = 8
        """
        path = joinpath(root, "U235_nth.toml")
        write(path, text)

        configuration = load_configuration(path; data_directory = root)
        @test configuration.method isa TemperatureRatio
        # Resolved through the manifest, and to the temperature-ratio table.
        @test occursin("R_T_vs_A_H_segmented", configuration.method.curve.source)
        @test !occursin("pivots", configuration.method.curve.source)
        @test configuration.method.curve(140) ≈ 1.2 rtol = RTOL

        # A label the manifest does not list is refused, naming the key at fault.
        thrown = try
            write(
                path,
                replace(
                    text,
                    "curve_label = \"systematic_trend\"" => "curve_label = \"absent\"",
                ),
            )
            load_configuration(path; data_directory = root)
            nothing
        catch err
            err
        end
        @test thrown isa ArgumentError
        @test occursin("partition.curve_label", thrown.msg)

        # A system whose own directory holds no staged run is refused rather than falling back.
        write(
            path,
            replace(text, "subdirectory = \"U235_nth\"" => "subdirectory = \"U235_nres\""),
        )
        @test_throws ArgumentError load_configuration(path; data_directory = root)

        # But a run of one system may take its curve from another system's manifest. The
        # resonance-channel case does exactly that, and it must not be turned into an error.
        mkpath(joinpath(root, "U235_nres"))
        write(
            path,
            replace(
                text,
                "target_A = 235\ntarget_Z = 92\nchannel = \"nth\"\nincident_energy = 2.53e-8" => "target_A = 235\ntarget_Z = 92\nchannel = \"nres\"\nincident_energy = 5.8013e-4",
            ),
        )
        crossed = load_configuration(path; data_directory = root)
        @test system_label(crossed.system) == "U235_nres"
        @test occursin("U235_nth", crossed.method.curve.source)
    end
end

@testset "a curve extracted on another domain is refused" begin
    isfile(MASS_EXCESS_FILE) || return nothing
    mktempdir() do root
        cp(MASS_EXCESS_FILE, joinpath(root, "MASSES"))
        system = joinpath(root, "U235_nth")
        mkpath(system)
        # The fixture's domain as this configuration builds it: BSFG, five charges, the mean
        # charge distribution, the mass table staged as MASSES.
        record(; kwargs...) = begin
            fields = merge(
                (
                    level_density_model = "BSFG",
                    ratio_averaging = "charge_resolved",
                    charges_per_mass = 5,
                    charge_model = charge_model_label(mean_charge_distribution()),
                    mass_table = "MASSES",
                ),
                values(kwargs),
            )
            "[domain]\n" *
            join(("$k = $(repr(v))" for (k, v) in pairs(fields)), "\n") *
            "\npackage_version = \"0.1.4\"\n\n[run]"
        end
        text = """
        [system]
        target_A = 235
        target_Z = 92
        channel = "nth"
        incident_energy = 2.53e-8

        [fragmentation]
        charges_per_mass = 5
        heavy_mass_min = 118
        heavy_mass_max = 160
        charge_distribution_file = "mean"

        [energetics]
        mass_excess_file = "MASSES"
        kinetic_energy_min = 150.0
        kinetic_energy_max = 190.0
        kinetic_energy_step = 20.0

        [level_density]
        model = "BSFG"

        [evaporation]
        cross_section = "variable"
        maximum_temperature = 2.0

        [partition]
        method = "temperature_ratio"
        subdirectory = "U235_nth"
        curve_label = "systematic_trend"

        [output]
        significant_digits = 8
        """
        path = joinpath(root, "U235_nth.toml")
        write(path, text)
        load_with(fields...) = begin
            rm(joinpath(system, "manifest_run.toml"); force = true)
            _manifest(system, "[run]" => record(; fields...))
            load_configuration(path; data_directory = root)
        end

        # The same domain loads without a word.
        @test_logs load_with() isa Any
        @test load_with().method isa TemperatureRatio

        # Each recorded property that differs from this configuration is refused by name.
        for (key, value) in (
            :level_density_model => "GC",
            :charges_per_mass => 7,
            :charge_model => "Wahl1988(U235T, Z_F = 92, A_F = 236)",
            :mass_table => "mass_excess_ame2016.dat",
        )
            thrown = try
                load_with(key => value)
                nothing
            catch err
                err
            end
            @test thrown isa ArgumentError
            @test occursin("[domain] $key", thrown.msg)
        end

        # An effective-ratio inversion loads, with a warning that it does not invert this
        # partition exactly.
        @test_logs (:warn, r"effective level density ratio") match_mode = :any load_with(
            :ratio_averaging => "ratio_of_means",
        )
    end
end
