# The validator is the promise that a pipeline cannot start from a configuration it cannot
# honour. These exercise it on files written here, so they need no input data, and each asserts
# that the message names the key at fault — a rejection that does not say which key is nearly as
# unhelpful as no rejection at all.

const MINIMAL_CONFIGURATION = """
[system]
target_A = 252
target_Z = 98
reaction = "0,f"
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
prescription = "BSFG"

[evaporation]
cross_section = "variable"
maximum_temperature = 2.0

[partition]
method = "excitation_ratio"
segments_file = "SEGMENTS"

[output]
subdirectory = "probe"
digits = 8
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
        # The two files the fixture names. Segments are the comma-delimited form read_segments
        # takes; the mass table is the shipped evaluation where it is present.
        write(
            joinpath(directory, "SEGMENTS"),
            "A_H,value,uncertainty\n126,0.5,0.0\n150,0.607,0.0\n174,0.865,0.0\n",
        )
        if isfile(MASS_EXCESS_FILE)
            cp(MASS_EXCESS_FILE, joinpath(directory, "MASSES"))
        else
            write(joinpath(directory, "MASSES"), "Z A element Δ σΔ\n")
        end

        if isfile(MASS_EXCESS_FILE)
            configuration = _configuration(directory)
            @test configuration isa Configuration
            @test configuration.system.compound == Nuclide(98, 252)
            @test configuration.system.spontaneous
            @test configuration.heavy_masses == 126:174
            @test configuration.charges_per_mass == 5
            @test configuration.kinetic_energies == 130.0:10.0:230.0
            @test configuration.model isa BackShiftedFermiGas
            @test configuration.cross_section isa VariableCrossSection
            @test configuration.method isa ExcitationEnergyRatio
            @test configuration.yields === nothing        # no [yield] section in the fixture
            @test configuration.spectrum_grid === nothing  # nor a [spectrum] one
            @test configuration.subdirectory == "probe"
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
            @test occursin("LA-13928", systematic.polarization.source)
            @test haskey(systematic.polarization.polarization, 126)
            @test haskey(systematic.polarization.polarization, 174)
        end

        # Every rejection names its key.
        @test _rejects(
            directory,
            "system.reaction",
            "reaction = \"0,f\"" => "reaction = \"p,f\"",
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
            "level_density.prescription",
            "prescription = \"BSFG\"" => "prescription = \"Fermi\"",
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
        @test _rejects(directory, "output.digits", "digits = 8" => "digits = 40")
        @test _rejects(
            directory,
            "output.subdirectory",
            "subdirectory = \"probe\"" => "subdirectory = \"a/b\"",
        )

        # A named input that is not there says which file, not which key alone.
        @test _rejects(
            directory,
            "does not exist",
            "segments_file = \"SEGMENTS\"" => "segments_file = \"ABSENT\"",
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
    end
end
