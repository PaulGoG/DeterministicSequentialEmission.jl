"""
    Configuration

A complete run: the fissioning system, the fragmentation domain, the level density and
evaporation models, the TXE partition, and the yield weights. Built by
[`load_configuration`](@ref) from a TOML file; every field is validated on load, and a section or
key outside [`CONFIGURATION_KEYS`](@ref) is refused. The output location is not a field: a run
writes to the directory named by [`system_label`](@ref), under the identifier built by
[`run_identifier`](@ref).
"""
struct Configuration
    system::FissioningSystem
    mass_excess_table::MassExcessTable
    charge_distribution::ChargeDistribution
    heavy_masses::UnitRange{Int}
    charges_per_mass::Int
    kinetic_energies::StepRangeLen{Float64}
    model::LevelDensityModel
    cross_section::EvaporationCrossSection
    maximum_temperature::Float64
    method::PartitionMethod
    fragment_yield::Union{FragmentYield, Nothing}
    spectrum_grid::Union{Vector{Float64}, Nothing}
    maxwellian_temperature::Float64
    spectrum_pair_weighting::SpectrumWeighting
    spectrum_average::SpectrumAverage
    significant_digits::Int
    source::String
    data_directory::String
end

function Base.show(io::IO, configuration::Configuration)
    print(
        io,
        "Configuration(",
        system_label(configuration.system),
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

"""
    CONFIGURATION_KEYS

The sections a configuration may hold and the keys each may carry. Any other section or key is
refused on load, so a misspelt optional key cannot silently take its default.
"""
const CONFIGURATION_KEYS = Dict(
    "system" => ("target_A", "target_Z", "channel", "incident_energy"),
    "energetics" => (
        "mass_excess_file",
        "kinetic_energy_min",
        "kinetic_energy_max",
        "kinetic_energy_step",
    ),
    "fragmentation" => (
        "charges_per_mass",
        "heavy_mass_min",
        "heavy_mass_max",
        "charge_distribution_file",
    ),
    "level_density" => ("model", "shell_correction_file", "deformed_branch"),
    "evaporation" => ("cross_section", "maximum_temperature"),
    "partition" => (
        "method",
        "subdirectory",
        "curve_label",
        "excitation_ratio_file",
        "deformation_energy_file",
    ),
    "yield" => (
        "total",
        "mass_energy_yield_file",
        "mass_yield_file",
        "mean_kinetic_energy_file",
        "kinetic_energy_dispersion",
        "kinetic_energy_dispersion_file",
    ),
    "spectrum" => (
        "enabled",
        "maxwellian_temperature",
        "energy_min",
        "energy_max",
        "points",
        "pair_weighting",
        "average",
    ),
    "output" => ("significant_digits",),
)

function _reject_unknown_keys(document::AbstractDict)
    for (name, section) in document
        haskey(CONFIGURATION_KEYS, name) || throw(
            ArgumentError(
                "the configuration has an unknown section [$name]; expected one of \
                 $(join(sort!(collect(keys(CONFIGURATION_KEYS))), ", "))",
            ),
        )
        section isa AbstractDict || continue
        allowed = CONFIGURATION_KEYS[name]
        for key in keys(section)
            key in allowed || throw(
                ArgumentError("$name.$key is not a configuration key; [$name] takes \
                     $(join(allowed, ", "))"),
            )
        end
    end
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

Read and validate a run configuration. Relative input paths in the file are resolved against
`data_directory`, `data/` beside the project by default; absolute paths are used as written.
"""
function load_configuration(
    path::AbstractString;
    data_directory::AbstractString = joinpath(PACKAGE_ROOT, "data"),
)
    isfile(path) || throw(ArgumentError("no configuration at $path"))
    document = TOML.parsefile(path)
    _reject_unknown_keys(document)

    system = _load_system(document)
    energetics_section = _section(document, "energetics")
    mass_excess_table = read_mass_excess_table(
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
    heavy_masses, charges_per_mass, charge_distribution =
        _load_fragmentation(document, system, mass_excess_table, data_directory)
    kinetic_energies = _load_kinetic_energies(energetics_section)
    model = _load_level_density(document, mass_excess_table, data_directory)
    cross_section, maximum_temperature = _load_evaporation(document)
    method = _load_partition(document, model, data_directory)
    fragment_yield = _load_yield(
        document,
        system,
        charge_distribution,
        heavy_masses,
        charges_per_mass,
        kinetic_energies,
        data_directory,
    )
    spectrum_grid, maxwellian, spectrum_pair_weighting, spectrum_average =
        _load_spectrum(document)
    significant_digits = _load_output(document)

    return Configuration(
        system,
        mass_excess_table,
        charge_distribution,
        heavy_masses,
        charges_per_mass,
        kinetic_energies,
        model,
        cross_section,
        maximum_temperature,
        method,
        fragment_yield,
        spectrum_grid,
        maxwellian,
        spectrum_pair_weighting,
        spectrum_average,
        significant_digits,
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
    channel = _choice(section, "channel", CHANNELS, "system.channel")
    energy = Float64(
        _bounded(
            _value(section, "incident_energy", Real, "system.incident_energy", 0.0),
            "system.incident_energy";
            min = 0,
            max = 20,
        ),
    )
    target = Nuclide(Z, A)
    if channel == "sf"
        iszero(energy) || throw(
            ArgumentError("system.channel = \"sf\" is spontaneous fission and requires \
 system.incident_energy = 0, got $energy MeV"),
        )
        return spontaneous_fission(target)
    end
    return neutron_induced_fission(target, energy, channel)
end

function _load_fragmentation(
    document::AbstractDict,
    system::FissioningSystem,
    mass_excess_table::MassExcessTable,
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
    distribution = if source == "wahl"
        # The Zₚ model, for a system without an evaluated ΔZ(A), σ_Z(A), tabulated over the
        # domain's heavy masses; the resolved layer is recorded in the distribution's source.
        build_charge_distribution(mass_excess_table, system, lower:upper)
    elseif source == "mean"
        mean_charge_distribution()
    else
        read_charge_distribution(
            _input(data_directory, source, "fragmentation.charge_distribution_file"),
        )
    end
    return lower:upper, charges, distribution
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
    mass_excess_table::MassExcessTable,
    data_directory::AbstractString,
)
    section = _section(document, "level_density")
    model = _choice(section, "model", LEVEL_DENSITY_MODELS, "level_density.model", "BSFG")
    model == "BSFG" && return BackShiftedFermiGas(mass_excess_table)
    corrections = read_shell_correction_table(
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
                _value(
                    section,
                    "deformation_energy_file",
                    String,
                    "partition.deformation_energy_file",
                ),
                "partition.deformation_energy_file",
            ),
        )
        return ScissionModelling(deformation, model)
    end
    if method == "excitation_ratio"
        # `E*_H/TXE (A_H)` has no producing run record, so the file is named directly.
        curve = read_segmented_curve(
            _input(
                data_directory,
                _value(
                    section,
                    "excitation_ratio_file",
                    String,
                    "partition.excitation_ratio_file",
                ),
                "partition.excitation_ratio_file",
            );
            spec = EXCITATION_RATIO_CURVE_SPEC,
        )
        return ExcitationEnergyRatio(curve)
    end

    # `R_T(A_H)` is read through the producing run's manifest, not a configured path: that run
    # writes `R_T` and `r_ν` side by side in the same layout, and the manifest identifies each.
    subdirectory = _value(section, "subdirectory", String, "partition.subdirectory")
    isempty(subdirectory) &&
        throw(ArgumentError("partition.subdirectory must not be empty"))
    manifest = read_temperature_ratio_manifest(
        staged_manifest(joinpath(data_directory, subdirectory)),
    )
    label = _value(section, "curve_label", String, "partition.curve_label")
    labels = curve_labels(manifest)
    label in labels || throw(
        ArgumentError(
            "partition.curve_label = $(repr(label)) is not in $(manifest.source), which \
 lists $(join(map(repr, labels), ", "))",
        ),
    )
    curve = read_segmented_curve(temperature_ratio_path(manifest, label))
    return TemperatureRatio(curve, model)
end

function _load_yield(
    document::AbstractDict,
    system::FissioningSystem,
    charge_distribution::ChargeDistribution,
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
    joint = haskey(section, "mass_energy_yield_file")
    marginal = haskey(section, "mass_yield_file")
    joint || marginal || return nothing
    joint &&
        marginal &&
        throw(
            ArgumentError(
                "yield.mass_energy_yield_file and yield.mass_yield_file are alternatives: the \
 second reconstructs a joint distribution from its marginals and is for \
 systems where the first does not exist",
            ),
        )

    measured = if joint
        read_mass_energy_yield(
            _input(
                data_directory,
                _value(
                    section,
                    "mass_energy_yield_file",
                    String,
                    "yield.mass_energy_yield_file",
                ),
                "yield.mass_energy_yield_file",
            ),
        )
    else
        # One width unless a measured table is given; the width is then the table's fallback
        # outside its range.
        default = Float64(
            _bounded(
                _value(
                    section,
                    "kinetic_energy_dispersion",
                    Real,
                    "yield.kinetic_energy_dispersion",
                    DEFAULT_KINETIC_ENERGY_DISPERSION,
                ),
                "yield.kinetic_energy_dispersion";
                min = 0,
                max = 100,
                exclusive_min = true,
            ),
        )
        dispersion = if haskey(section, "kinetic_energy_dispersion_file")
            read_kinetic_energy_dispersion(
                _input(
                    data_directory,
                    _value(
                        section,
                        "kinetic_energy_dispersion_file",
                        String,
                        "yield.kinetic_energy_dispersion_file",
                    ),
                    "yield.kinetic_energy_dispersion_file",
                );
                default = default,
            )
        else
            uniform_kinetic_energy_dispersion(default)
        end
        factorized_yield(
            _input(
                data_directory,
                _value(section, "mass_yield_file", String, "yield.mass_yield_file"),
                "yield.mass_yield_file",
            ),
            _input(
                data_directory,
                _value(
                    section,
                    "mean_kinetic_energy_file",
                    String,
                    "yield.mean_kinetic_energy_file",
                ),
                "yield.mean_kinetic_energy_file",
            ),
            kinetic_energies,
            system.compound.A;
            dispersion = dispersion,
            total = total,
        )
    end

    domain = fragmentation_domain(
        system,
        charge_distribution,
        heavy_masses;
        charges_per_mass = charges_per_mass,
    )
    return build_fragment_yield(domain, charge_distribution, measured; total = total)
end

function _load_spectrum(document::AbstractDict)
    haskey(document, "spectrum") ||
        return nothing, 1.42, MultiplicityWeighting(), PerNeutronAverage()
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
    pair_weighting = _spectrum_pair_weighting(
        _choice(
            section,
            "pair_weighting",
            ("multiplicity", "excitation"),
            "spectrum.pair_weighting",
            "multiplicity",
        ),
    )
    average = _spectrum_average(
        _choice(
            section,
            "average",
            ("per_fission", "per_neutron"),
            "spectrum.average",
            "per_neutron",
        ),
    )
    _value(section, "enabled", Bool, "spectrum.enabled", true) ||
        return nothing, maxwellian, pair_weighting, average
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
    # Logarithmic grid: the spectrum spans four decades, and a linear grid does not resolve the
    # sub-MeV shape.
    grid = 10 .^ range(log10(lower), log10(upper); length = points)
    return grid, maxwellian, pair_weighting, average
end

_spectrum_pair_weighting(kind::String) =
    kind == "excitation" ? ExcitationWeighting() : MultiplicityWeighting()
_spectrum_average(kind::String) =
    kind == "per_neutron" ? PerNeutronAverage() : PerFissionAverage()

function _load_output(document::AbstractDict)
    section = _section(document, "output")
    return Int(
        _bounded(
            _value(section, "significant_digits", Int, "output.significant_digits", 8),
            "output.significant_digits";
            min = 1,
            max = 15,
        ),
    )
end
