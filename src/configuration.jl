"""
    Configuration

A complete run: the fissioning system, the fragmentation domain, the level density and
evaporation prescriptions, the TXE partition, the yield weights, and where the results go.

Built by [`load_configuration`](@ref) from a TOML file. Every field is validated on load, so a
pipeline cannot start from a configuration it is unable to honour, and the values that reach the
model are the ones the file states rather than defaults filled in silently along the way.
"""
struct Configuration
    system::FissioningSystem
    masses::MassTable
    polarization::ChargePolarization
    heavy_masses::UnitRange{Int}
    charges_per_mass::Int
    kinetic_energies::StepRangeLen{Float64}
    model::LevelDensityModel
    cross_section::EvaporationCrossSection
    maximum_temperature::Float64
    method::PartitionMethod
    yields::Union{FragmentYields, Nothing}
    spectrum_grid::Union{Vector{Float64}, Nothing}
    maxwellian_temperature::Float64
    digits::Int
    subdirectory::String
    source::String
    data_directory::String
end

function Base.show(io::IO, configuration::Configuration)
    print(
        io,
        "Configuration(",
        configuration.subdirectory,
        ", A_H ∈ ",
        configuration.heavy_masses,
        ", TKE ∈ ",
        first(configuration.kinetic_energies),
        ":",
        step(configuration.kinetic_energies),
        ":",
        last(configuration.kinetic_energies),
        " MeV)",
    )
    return nothing
end

# Reading one key with its type, its bounds and its allowed values, so that a malformed
# configuration names the offending key rather than failing somewhere downstream.

function _section(document::AbstractDict, name::AbstractString)
    haskey(document, name) ||
        throw(ArgumentError("the configuration has no [$name] section"))
    section = document[name]
    section isa AbstractDict ||
        throw(ArgumentError("[$name] must be a table, got $(typeof(section))"))
    return section
end

function _value(
    section::AbstractDict,
    key::AbstractString,
    ::Type{T},
    label::AbstractString,
    default = nothing,
) where {T}
    if !haskey(section, key)
        default === nothing && throw(ArgumentError("$label is required"))
        return default::T
    end
    raw = section[key]
    raw isa T ||
        throw(ArgumentError("$label must be $(T), got $(repr(raw)) of type $(typeof(raw))"))
    return raw
end

function _choice(
    section::AbstractDict,
    key::AbstractString,
    allowed,
    label::AbstractString,
    default = nothing,
)
    raw = _value(section, key, String, label, default)
    raw in allowed ||
        throw(ArgumentError("$label must be one of $(join(map(repr, allowed), ", ")), \
                       got $(repr(raw))"))
    return raw
end

function _bounded(
    value::Real,
    label::AbstractString;
    min = -Inf,
    max = Inf,
    exclusive_min::Bool = false,
)
    isfinite(value) || throw(ArgumentError("$label must be finite, got $value"))
    lower_ok = exclusive_min ? value > min : value >= min
    lower_ok || throw(
        ArgumentError("$label must be $(exclusive_min ? "above" : "at least") $min, \
                       got $value"),
    )
    value <= max || throw(ArgumentError("$label must be at most $max, got $value"))
    return value
end

function _input(
    data_directory::AbstractString,
    relative::AbstractString,
    label::AbstractString,
)
    path = isabspath(relative) ? relative : joinpath(data_directory, relative)
    isfile(path) || throw(ArgumentError("$label names a file that does not exist: $path"))
    return path
end

"""
    load_configuration(path; data_directory) -> Configuration

Read and validate a run configuration.

`data_directory` is the root every input path in the file is taken relative to, `data/` beside
the project by default, so a configuration is portable between checkouts and a data tree can
live anywhere. Absolute paths in the file are honoured as written.
"""
function load_configuration(
    path::AbstractString;
    data_directory::AbstractString = joinpath(PACKAGE_ROOT, "data"),
)
    isfile(path) || throw(ArgumentError("no configuration at $path"))
    document = TOML.parsefile(path)

    system = _load_system(document)
    energetics_section = _section(document, "energetics")
    masses = read_mass_table(
        _input(
            data_directory,
            _value(
                energetics_section,
                "mass_excess_file",
                String,
                "energetics.mass_excess_file",
            ),
            "energetics.mass_excess_file",
        ),
    )
    heavy_masses, charges_per_mass, polarization =
        _load_fragmentation(document, system, masses, data_directory)
    kinetic_energies = _load_kinetic_energies(energetics_section)
    model = _load_level_density(document, masses, data_directory)
    cross_section, maximum_temperature = _load_evaporation(document)
    method = _load_partition(document, model, data_directory)
    yields = _load_yields(
        document,
        system,
        polarization,
        heavy_masses,
        charges_per_mass,
        kinetic_energies,
        data_directory,
    )
    spectrum_grid, maxwellian = _load_spectrum(document)
    digits, subdirectory = _load_output(document)

    return Configuration(
        system,
        masses,
        polarization,
        heavy_masses,
        charges_per_mass,
        kinetic_energies,
        model,
        cross_section,
        maximum_temperature,
        method,
        yields,
        spectrum_grid,
        maxwellian,
        digits,
        subdirectory,
        abspath(path),
        abspath(data_directory),
    )
end

function _load_system(document::AbstractDict)
    section = _section(document, "system")
    A = Int(
        _bounded(
            _value(section, "target_A", Int, "system.target_A"),
            "system.target_A";
            min = 2,
            max = 400,
        ),
    )
    Z = Int(
        _bounded(
            _value(section, "target_Z", Int, "system.target_Z"),
            "system.target_Z";
            min = 1,
            max = 118,
        ),
    )
    Z < A || throw(ArgumentError("system.target_Z = $Z must be below system.target_A = $A"))
    reaction = _choice(section, "reaction", ("0,f", "n,f"), "system.reaction")
    energy = Float64(
        _bounded(
            _value(section, "incident_energy", Real, "system.incident_energy", 0.0),
            "system.incident_energy";
            min = 0,
            max = 20,
        ),
    )
    target = Nuclide(Z, A)
    if reaction == "0,f"
        iszero(energy) || throw(
            ArgumentError("system.reaction = \"0,f\" is spontaneous fission and requires \
                           system.incident_energy = 0, got $energy MeV"),
        )
        return spontaneous_fission(target)
    end
    return neutron_induced_fission(target, energy)
end

function _load_fragmentation(
    document::AbstractDict,
    system::FissioningSystem,
    masses::MassTable,
    data_directory::AbstractString,
)
    section = _section(document, "fragmentation")
    charges = Int(
        _bounded(
            _value(section, "charges_per_mass", Int, "fragmentation.charges_per_mass"),
            "fragmentation.charges_per_mass";
            min = 1,
            max = 15,
        ),
    )
    isodd(charges) || throw(
        ArgumentError(
            "fragmentation.charges_per_mass must be odd so the charges are symmetric \
                       about Zₚ, got $charges",
        ),
    )
    lower = _value(section, "heavy_mass_min", Int, "fragmentation.heavy_mass_min")
    upper = _value(section, "heavy_mass_max", Int, "fragmentation.heavy_mass_max")
    lower < upper ||
        throw(ArgumentError("fragmentation.heavy_mass_min = $lower must be below \
                       fragmentation.heavy_mass_max = $upper"))
    2 * lower >= system.compound.A || throw(
        ArgumentError(
            "fragmentation.heavy_mass_min = $lower is below the symmetric split \
                       A₀/2 = $(system.compound.A / 2) and would label light fragments heavy",
        ),
    )
    upper < system.compound.A || throw(
        ArgumentError(
            "fragmentation.heavy_mass_max = $upper must be below the compound mass \
                       A₀ = $(system.compound.A)",
        ),
    )

    source = _value(
        section,
        "charge_distribution_file",
        String,
        "fragmentation.charge_distribution_file",
    )
    polarization = if source == "wahl"
        # The Zₚ systematics rather than a table: the route for a system with no evaluated
        # ΔZ(A), rms(A) of its own. Tabulated over exactly the masses the domain will ask for.
        systematics = WahlSystematics(masses, system)
        systematics === nothing && throw(
            ArgumentError(
                "fragmentation.charge_distribution_file = \"wahl\" needs Sₙ of the \
                           compound nucleus, which the mass table cannot form",
            ),
        )
        wahl_charge_polarization(systematics, lower:upper)
    elseif source == "mean"
        mean_charge_polarization()
    else
        read_charge_polarization(
            _input(data_directory, source, "fragmentation.charge_distribution_file"),
        )
    end
    return lower:upper, charges, polarization
end

function _load_kinetic_energies(section::AbstractDict)
    lower = Float64(
        _bounded(
            _value(section, "kinetic_energy_min", Real, "energetics.kinetic_energy_min"),
            "energetics.kinetic_energy_min";
            min = 0,
            exclusive_min = true,
        ),
    )
    upper = Float64(
        _bounded(
            _value(section, "kinetic_energy_max", Real, "energetics.kinetic_energy_max"),
            "energetics.kinetic_energy_max";
            min = 0,
            exclusive_min = true,
        ),
    )
    step = Float64(
        _bounded(
            _value(section, "kinetic_energy_step", Real, "energetics.kinetic_energy_step"),
            "energetics.kinetic_energy_step";
            min = 0,
            max = 50,
            exclusive_min = true,
        ),
    )
    lower < upper ||
        throw(ArgumentError("energetics.kinetic_energy_min = $lower must be below \
                       energetics.kinetic_energy_max = $upper"))
    return lower:step:upper
end

function _load_level_density(
    document::AbstractDict,
    masses::MassTable,
    data_directory::AbstractString,
)
    section = _section(document, "level_density")
    prescription = _choice(
        section,
        "prescription",
        ("BSFG", "GC"),
        "level_density.prescription",
        "BSFG",
    )
    prescription == "BSFG" && return BackShiftedFermiGas(masses)
    corrections = read_shell_corrections(
        _input(
            data_directory,
            _value(
                section,
                "shell_correction_file",
                String,
                "level_density.shell_correction_file",
            ),
            "level_density.shell_correction_file",
        ),
    )
    deformed =
        _value(section, "deformed_branch", Bool, "level_density.deformed_branch", true)
    return GilbertCameron(corrections; deformed_branch = deformed)
end

function _load_evaporation(document::AbstractDict)
    section = _section(document, "evaporation")
    kind = _choice(
        section,
        "cross_section",
        ("variable", "constant"),
        "evaporation.cross_section",
        "variable",
    )
    maximum_temperature = Float64(
        _bounded(
            _value(
                section,
                "maximum_temperature",
                Real,
                "evaporation.maximum_temperature",
                2.0,
            ),
            "evaporation.maximum_temperature";
            min = 0,
            max = 10,
            exclusive_min = true,
        ),
    )
    cross_section = kind == "variable" ? VariableCrossSection() : ConstantCrossSection()
    return cross_section, maximum_temperature
end

function _load_partition(
    document::AbstractDict,
    model::LevelDensityModel,
    data_directory::AbstractString,
)
    section = _section(document, "partition")
    method = _choice(
        section,
        "method",
        ("temperature_ratio", "excitation_ratio", "scission"),
        "partition.method",
    )
    if method == "scission"
        deformation = read_extra_deformation(
            _input(
                data_directory,
                _value(section, "deformation_file", String, "partition.deformation_file"),
                "partition.deformation_file",
            ),
        )
        return ScissionModelling(deformation, model)
    end
    segments = read_segments(
        _input(
            data_directory,
            _value(section, "segments_file", String, "partition.segments_file"),
            "partition.segments_file",
        ),
    )
    method == "temperature_ratio" && return TemperatureRatio(segments, model)
    return ExcitationEnergyRatio(segments)
end

function _load_yields(
    document::AbstractDict,
    system::FissioningSystem,
    polarization::ChargePolarization,
    heavy_masses::UnitRange{Int},
    charges_per_mass::Integer,
    kinetic_energies,
    data_directory::AbstractString,
)
    haskey(document, "yield") || return nothing
    section = _section(document, "yield")
    total = Float64(
        _bounded(
            _value(section, "total", Real, "yield.total", DEFAULT_YIELD_TOTAL),
            "yield.total";
            min = 0,
            exclusive_min = true,
        ),
    )
    joint = haskey(section, "file")
    marginal = haskey(section, "mass_yield_file")
    joint || marginal || return nothing
    joint &&
        marginal &&
        throw(
            ArgumentError(
                "yield.file and yield.mass_yield_file are alternatives: the second \
                       reconstructs a joint distribution from its marginals and is for systems \
                       where the first does not exist",
            ),
        )

    experimental = if joint
        read_experimental_yield(
            _input(
                data_directory,
                _value(section, "file", String, "yield.file"),
                "yield.file",
            ),
        )
    else
        width = Float64(
            _bounded(
                _value(
                    section,
                    "kinetic_energy_width",
                    Real,
                    "yield.kinetic_energy_width",
                    10.0,
                ),
                "yield.kinetic_energy_width";
                min = 0,
                max = 100,
                exclusive_min = true,
            ),
        )
        factorized_yield(
            _input(
                data_directory,
                _value(section, "mass_yield_file", String, "yield.mass_yield_file"),
                "yield.mass_yield_file",
            ),
            _input(
                data_directory,
                _value(section, "kinetic_energy_file", String, "yield.kinetic_energy_file"),
                "yield.kinetic_energy_file",
            ),
            kinetic_energies,
            system.compound.A;
            width = width,
            total = total,
        )
    end

    domain = fragmentation_domain(
        system,
        polarization,
        heavy_masses;
        charges_per_mass = charges_per_mass,
    )
    return build_fragment_yields(domain, polarization, experimental; total = total)
end

function _load_spectrum(document::AbstractDict)
    haskey(document, "spectrum") || return nothing, 1.42
    section = _section(document, "spectrum")
    maxwellian = Float64(
        _bounded(
            _value(
                section,
                "maxwellian_temperature",
                Real,
                "spectrum.maxwellian_temperature",
                1.42,
            ),
            "spectrum.maxwellian_temperature";
            min = 0,
            max = 5,
            exclusive_min = true,
        ),
    )
    _value(section, "enabled", Bool, "spectrum.enabled", true) || return nothing, maxwellian
    lower = Float64(
        _bounded(
            _value(section, "energy_min", Real, "spectrum.energy_min", 0.05),
            "spectrum.energy_min";
            min = 0,
            max = 100,
            exclusive_min = true,
        ),
    )
    upper = Float64(
        _bounded(
            _value(section, "energy_max", Real, "spectrum.energy_max", 20.0),
            "spectrum.energy_max";
            min = 0,
            max = 100,
            exclusive_min = true,
        ),
    )
    lower < upper || throw(ArgumentError("spectrum.energy_min = $lower must be below \
                       spectrum.energy_max = $upper"))
    points = Int(
        _bounded(
            _value(section, "points", Int, "spectrum.points", 240),
            "spectrum.points";
            min = 8,
            max = 4000,
        ),
    )
    # Logarithmic, because the spectrum is compared over four decades and a linear grid resolves
    # none of the sub-MeV shape.
    return 10 .^ range(log10(lower), log10(upper); length = points), maxwellian
end

function _load_output(document::AbstractDict)
    section = _section(document, "output")
    digits = Int(
        _bounded(
            _value(section, "digits", Int, "output.digits", 8),
            "output.digits";
            min = 1,
            max = 15,
        ),
    )
    subdirectory = _value(section, "subdirectory", String, "output.subdirectory")
    isempty(subdirectory) && throw(ArgumentError("output.subdirectory must not be empty"))
    occursin('/', subdirectory) || occursin('\\', subdirectory) ?
    throw(ArgumentError("output.subdirectory is one directory name, not a path, \
                       got $(repr(subdirectory))")) : nothing
    return digits, subdirectory
end
