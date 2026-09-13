# The validator is the promise that a pipeline cannot start from a configuration it cannot
# honour. These exercise it on files written here, so they need no input data, and each asserts
# that the message names the key at fault — a rejection that does not say which key is nearly as
# unhelpful as no rejection at all.

const MINIMAL_CONFIGURATION = """
[system]
target_A = 252
target_Z = 98
channel = "sf"
incident_energy = 0.0

[fragmentation]
charges_per_mass = 5
heavy_mass_min = 126
heavy_mass_max = 174
charge_distribution_file = "mean"

[energetics]
mass_excess_file = "MASSES"
kinetic_energy_min = 130.0
kinetic_energy_max = 230.0
kinetic_energy_step = 10.0

[level_density]
model = "BSFG"

[evaporation]
cross_section = "variable"
maximum_temperature = 2.0

[partition]
method = "excitation_ratio"
segmented_curve_file = "SEGMENTS"

[output]
significant_digits = 8
"""

"Write a configuration with `replacements` applied, and load it against `directory`."
function _configuration(directory::AbstractString, replacements::Pair{String, String}...)
    text = MINIMAL_CONFIGURATION
    for (from, to) in replacements
        occursin(from, text) || error("the fixture has no $(repr(from)) to replace")
        text = replace(text, from => to)
    end
    path = joinpath(directory, "probe.toml")
    write(path, text)
    return load_configuration(path; data_directory = directory)
end

function _rejects(
    directory::AbstractString,
    key::AbstractString,
    replacements::Pair{String, String}...,
)
    thrown = try
        _configuration(directory, replacements...)
        nothing
    catch err
        err
    end
    thrown isa ArgumentError || return false
    return occursin(key, thrown.msg)
end

@testset "run configuration" begin
    mktempdir() do directory
        # The two files the fixture names. Segments are the comma-delimited form read_segmented_curve
        # takes; the mass table is the shipped evaluation where it is present.
        write(
            joinpath(directory, "SEGMENTS"),
            "A_H,R_T,R_T_uncertainty\n126,0.5,0.0\n150,0.607,0.0\n174,0.865,0.0\n",
        )
        if isfile(MASS_EXCESS_FILE)
            cp(MASS_EXCESS_FILE, joinpath(directory, "MASSES"))
        else
            write(
                joinpath(directory, "MASSES"),
                "Z A element mass_excess mass_excess_uncertainty\n",
            )
        end

        if isfile(MASS_EXCESS_FILE)
            configuration = _configuration(directory)
            @test configuration isa Configuration
            @test configuration.system.compound == Nuclide(98, 252)
            @test is_spontaneous(configuration.system)
            @test configuration.heavy_masses == 126:174
            @test configuration.charges_per_mass == 5
            @test configuration.kinetic_energies == 130.0:10.0:230.0
            @test configuration.model isa BackShiftedFermiGas
            @test configuration.cross_section isa VariableCrossSection
            @test configuration.method isa ExcitationEnergyRatio
            @test configuration.fragment_yield === nothing        # no [yield] section in the fixture
            @test configuration.spectrum_grid === nothing  # nor a [spectrum] one
            @test configuration.significant_digits == 8
            # The output directory is not a key: it follows from the system section.
            @test system_label(configuration.system) == "Cf252_sf"
            @test occursin("probe.toml", configuration.source)

            # A spectrum grid is logarithmic and spans exactly what it was asked for.
            with_spectrum = _configuration(
                directory,
                "[output]" => "[spectrum]\nenabled = true\nenergy_min = 0.1\nenergy_max = 10.0\npoints = 64\n\n[output]",
            )
            grid = with_spectrum.spectrum_grid
            @test grid !== nothing
            @test length(grid) == 64
            @test first(grid) ≈ 0.1
            @test last(grid) ≈ 10.0
            @test issorted(grid)
            @test grid[2] / grid[1] ≈ grid[end] / grid[end - 1]

            # "wahl" builds the systematics in place of a table, over exactly the domain's masses.
            systematic = _configuration(
                directory,
                "charge_distribution_file = \"mean\"" => "charge_distribution_file = \"wahl\"",
            )
            @test occursin("LA-13928", systematic.charge_distribution.source)
            @test haskey(systematic.charge_distribution.ΔZ, 126)
            @test haskey(systematic.charge_distribution.ΔZ, 174)
        end

        # Every rejection names its key.
        @test _rejects(
            directory,
            "system.channel",
            "channel = \"sf\"" => "channel = \"p,f\"",
        )
        @test _rejects(
            directory,
            "system.incident_energy",
            "incident_energy = 0.0" => "incident_energy = 2.53e-8",
        )
        @test _rejects(
            directory,
            "fragmentation.charges_per_mass",
            "charges_per_mass = 5" => "charges_per_mass = 4",
        )
        @test _rejects(
            directory,
            "fragmentation.heavy_mass_min",
            "heavy_mass_min = 126" => "heavy_mass_min = 100",
        )
        @test _rejects(
            directory,
            "fragmentation.heavy_mass_max",
            "heavy_mass_max = 174" => "heavy_mass_max = 260",
        )
        @test _rejects(
            directory,
            "energetics.kinetic_energy_min",
            "kinetic_energy_min = 130.0" => "kinetic_energy_min = 240.0",
        )
        @test _rejects(
            directory,
            "energetics.kinetic_energy_step",
            "kinetic_energy_step = 10.0" => "kinetic_energy_step = -1.0",
        )
        @test _rejects(
            directory,
            "level_density.model",
            "model = \"BSFG\"" => "model = \"Fermi\"",
        )
        @test _rejects(
            directory,
            "evaporation.cross_section",
            "cross_section = \"variable\"" => "cross_section = \"flat\"",
        )
        @test _rejects(
            directory,
            "partition.method",
            "method = \"excitation_ratio\"" => "method = \"guess\"",
        )
        @test _rejects(
            directory,
            "output.significant_digits",
            "significant_digits = 8" => "significant_digits = 40",
        )

        # A named input that is not there says which file, not which key alone.
        @test _rejects(
            directory,
            "does not exist",
            "segmented_curve_file = \"SEGMENTS\"" => "segmented_curve_file = \"ABSENT\"",
        )

        # A missing section, and a missing required key within one.
        @test_throws ArgumentError begin
            path = joinpath(directory, "headless.toml")
            write(path, "[system]\ntarget_A = 252\n")
            load_configuration(path; data_directory = directory)
        end
        @test_throws ArgumentError load_configuration(joinpath(directory, "absent.toml"))

        # A type error is a type error, not a silent coercion.
        @test _rejects(
            directory,
            "system.target_A",
            "target_A = 252" => "target_A = \"252\"",
        )
    end
end

@testset "shipped configurations" begin
    # Each is either loadable or refused for the one reason a bare clone permits: an input file
    # the licence does not let this package redistribute.
    directory = joinpath(dirname(@__DIR__), "config")
    files = sort(filter(name -> endswith(name, ".toml"), readdir(directory; join = true)))
    @test !isempty(files)
    for file in files
        outcome = try
            load_configuration(file)
        catch err
            err isa ArgumentError || rethrow()
            @test occursin("does not exist", err.msg)
            nothing
        end
        outcome === nothing || @test outcome isa Configuration
        # A configuration is named for the system it declares, so a file cannot drift from
        # what it runs. The name is read straight from the file rather than from the loaded
        # configuration, which a bare clone may not have.
        document = TOML.parsefile(file)
        declared = system_label(
            document["system"]["channel"] == "sf" ?
            spontaneous_fission(
                Nuclide(document["system"]["target_Z"], document["system"]["target_A"]),
            ) :
            neutron_induced_fission(
                Nuclide(document["system"]["target_Z"], document["system"]["target_A"]),
                document["system"]["incident_energy"],
                document["system"]["channel"],
            ),
        )
        @test basename(file) == declared * ".toml"
    end
end

@testset "run identifiers" begin
    # Every key a savename uses has an abbreviation, and no two keys share a token: one place
    # to look, and a reader can get from a directory name back to a configuration.
    tokens = collect(values(RUN_IDENTIFIER_ABBREVIATIONS))
    @test allunique(tokens)
    @test all(!isempty, tokens)

    mktempdir() do directory
        write(
            joinpath(directory, "SEGMENTS"),
            "A_H,R_T,R_T_uncertainty\n126,0.5,0.0\n150,0.607,0.0\n174,0.865,0.0\n",
        )
        isfile(MASS_EXCESS_FILE) || return nothing
        cp(MASS_EXCESS_FILE, joinpath(directory, "MASSES"))

        configuration = _configuration(directory)
        identifier = run_identifier(configuration)
        # Every token is `<abbreviation>=<value>`, and every abbreviation is in the table.
        abbreviations = Set(values(RUN_IDENTIFIER_ABBREVIATIONS))
        for piece in split(identifier, "_")
            key = first(split(piece, "="))
            occursin("=", piece) && @test key in abbreviations
        end
        @test occursin("ldm=BSFG", identifier)
        @test occursin("AHmin=126", identifier)
        @test occursin("AHmax=174", identifier)
        @test occursin("part=excitation_ratio", identifier)

        # Two runs differing only in a covered key are two identifiers, not one.
        wider = _configuration(directory, "charges_per_mass = 5" => "charges_per_mass = 7")
        @test run_identifier(wider) != identifier

        # The incident energy is covered: it changes the compound excitation, and two runs of
        # one system on one channel may differ in nothing else.
        # Capturing a neutron raises A₀ to 253, so the symmetric split moves off 126.
        slow = _configuration(
            directory,
            "channel = \"sf\"" => "channel = \"nfast\"",
            "incident_energy = 0.0" => "incident_energy = 1.0",
            "heavy_mass_min = 126" => "heavy_mass_min = 127",
        )
        fast = _configuration(
            directory,
            "channel = \"sf\"" => "channel = \"nfast\"",
            "incident_energy = 0.0" => "incident_energy = 2.0",
            "heavy_mass_min = 126" => "heavy_mass_min = 127",
        )
        @test run_identifier(slow) != run_identifier(fast)
    end
end
