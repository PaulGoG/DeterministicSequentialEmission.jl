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
    charge_distribution::ChargeModel
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
    yield_symmetrized::Bool
    yield_width::String
    yield_inputs::Dict{String, Any}
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
        "symmetrize",
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
    mass_excess_source = _value(
        energetics_section,
        "mass_excess_file",
        String,
        "energetics.mass_excess_file",
    )
    mass_excess_table = read_mass_excess_table(
        mass_excess_source == "ame2020" ? String(AME2020_MASS_EXCESS_FILE) :
        _input(data_directory, mass_excess_source, "energetics.mass_excess_file"),
    )
    heavy_masses, charges_per_mass, charge_distribution =
        _load_fragmentation(document, system, mass_excess_table, data_directory)
    kinetic_energies = _load_kinetic_energies(energetics_section)
    model = _load_level_density(document, mass_excess_table, data_directory)
    cross_section, maximum_temperature = _load_evaporation(document)
    method = _load_partition(
        document,
        model,
        data_directory;
        charges_per_mass = charges_per_mass,
        charge_distribution = charge_distribution,
        mass_excess_table = mass_excess_table,
    )
    fragment_yield, yield_symmetrized, yield_width, yield_inputs = _load_yield(
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
        yield_symmetrized,
        yield_width,
        yield_inputs,
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
        # Wahl's Zₚ model: the reaction's own 1988 parameters, else the 2002 systematics, else
        # the means. The model serves the domain directly; the layer is recorded in the metadata.
        charge_model(mass_excess_table, system)
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
    # "gilbert_cameron" is the table shipped with the domain package, Table III of the 1965 paper;
    # any other value is a path under the data directory.
    shell_source = _value(
        section,
        "shell_correction_file",
        String,
        "level_density.shell_correction_file",
        "gilbert_cameron",
    )
    corrections = read_shell_correction_table(
        shell_source == "gilbert_cameron" ?
        String(GILBERT_CAMERON_SHELL_CORRECTION_FILE) :
        _input(data_directory, shell_source, "level_density.shell_correction_file"),
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
    data_directory::AbstractString;
    charges_per_mass::Integer,
    charge_distribution::ChargeModel,
    mass_excess_table::MassExcessTable,
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
    _check_manifest_domain(
        manifest,
        model,
        charges_per_mass,
        charge_distribution,
        mass_excess_table,
    )
    curve = read_segmented_curve(temperature_ratio_path(manifest, label))
    return TemperatureRatio(curve, model)
end

# An R_T(A_H) curve inverts the partition on the domain it was extracted on. Where the run
# record says which domain that was, a different one here is refused: the curve would be applied
# to fragmentations, level densities or masses it was not derived from. A record written before
# the [domain] table existed can only be taken on trust.
function _check_manifest_domain(
    manifest::TemperatureRatioManifest,
    model::LevelDensityModel,
    charges_per_mass::Integer,
    charge_distribution::ChargeModel,
    mass_excess_table::MassExcessTable,
)
    domain = manifest.domain
    if domain === nothing
        @warn "$(manifest.source) carries no [domain] record, so the domain its R_T(A_H) was \
               extracted on cannot be checked against this configuration's" maxlog = 1
        return nothing
    end
    here = (
        level_density_model = level_density_label(model),
        deformed_branch = model isa GilbertCameron && model.deformed_branch,
        charges_per_mass = Int(charges_per_mass),
        charge_model = charge_model_label(charge_distribution),
        mass_table = basename(mass_excess_table.source),
    )
    for key in keys(here)
        recorded = getfield(domain, key)
        recorded == here[key] || throw(
            ArgumentError(
                "partition.subdirectory: $(manifest.source) was extracted with [domain] \
                 $key = $(repr(recorded)), and this configuration has $(repr(here[key])); \
                 the curve does not invert this partition",
            ),
        )
    end
    domain.ratio_averaging == "charge_resolved" || @warn "$(manifest.source) inverted \
        R_T(A_H) through an effective level density ratio, ratio_averaging = \
        $(repr(domain.ratio_averaging)); this partition applies it charge by charge, and the \
        round trip misses by up to 1e-2 in R_T" maxlog = 1
    return nothing
end

function _load_yield(
    document::AbstractDict,
    system::FissioningSystem,
    charge_distribution::ChargeModel,
    heavy_masses::UnitRange{Int},
    charges_per_mass::Integer,
    kinetic_energies,
    data_directory::AbstractString,
)
    haskey(document, "yield") || return nothing, false, "none", Dict{String, Any}()
    section = _section(document, "yield")
    # The two fragments of a split are one event, so pre-neutron data obey
    # Y(A) = Y(A₀ − A), and likewise ⟨TKE⟩ and σ_TKE. Measured tables that break the identity are
    # made to satisfy it by averaging the complements.
    symmetrize = _value(section, "symmetrize", Bool, "yield.symmetrize", true)
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
    joint || marginal || return nothing, false, "none", Dict{String, Any}()
    joint &&
        marginal &&
        throw(
            ArgumentError(
                "yield.mass_energy_yield_file and yield.mass_yield_file are alternatives: the \
 second reconstructs a joint distribution from its marginals and is for \
 systems where the first does not exist",
            ),
        )

    width = "joint"
    inputs = Dict{String, Any}()
    measured = if joint
        table = read_mass_energy_yield(
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
        inputs["mass_energy_yield_file"] = _value(
            section,
            "mass_energy_yield_file",
            String,
            "yield.mass_energy_yield_file",
        )
        symmetrize ? symmetrized_yield(table, system.compound.A) : table
    else
        # The width of P(TKE|A) is the one the experiment that measured ⟨TKE⟩(A) reported, in
        # the dispersion file. A stated single width is used only where the configuration names
        # one. Otherwise no width is assumed, and each mass is placed at its mean TKE.
        stated =
            haskey(section, "kinetic_energy_dispersion") ?
            Float64(
                _bounded(
                    _value(
                        section,
                        "kinetic_energy_dispersion",
                        Real,
                        "yield.kinetic_energy_dispersion",
                    ),
                    "yield.kinetic_energy_dispersion";
                    min = 0,
                    max = 100,
                    exclusive_min = true,
                ),
            ) : nothing
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
                default = stated,
            )
        elseif stated !== nothing
            uniform_kinetic_energy_dispersion(stated)
        else
            nothing
        end
        mass_path = _input(
            data_directory,
            _value(section, "mass_yield_file", String, "yield.mass_yield_file"),
            "yield.mass_yield_file",
        )
        energy_path = _input(
            data_directory,
            _value(
                section,
                "mean_kinetic_energy_file",
                String,
                "yield.mean_kinetic_energy_file",
            ),
            "yield.mean_kinetic_energy_file",
        )
        inputs["mass_yield_file"] =
            _value(section, "mass_yield_file", String, "yield.mass_yield_file")
        inputs["mean_kinetic_energy_file"] = _value(
            section,
            "mean_kinetic_energy_file",
            String,
            "yield.mean_kinetic_energy_file",
        )
        reconstructed = factorized_yield(
            mass_path,
            energy_path,
            kinetic_energies,
            system.compound.A;
            dispersion = dispersion,
            total = total,
            symmetrize = symmetrize,
        )
        width = if dispersion === nothing
            "none"
        elseif dispersion.source == "uniform"
            "stated"
        else
            inputs["kinetic_energy_dispersion_file"] = _value(
                section,
                "kinetic_energy_dispersion_file",
                String,
                "yield.kinetic_energy_dispersion_file",
            )
            inputs["masses_without_width"] = _masses_without_width(
                reconstructed,
                dispersion,
                system.compound.A,
                symmetrize,
            )
            _width_provenance(basename(dispersion.source), basename(energy_path))
        end
        reconstructed
    end

    domain = fragmentation_domain(
        system,
        charge_distribution,
        heavy_masses;
        charges_per_mass = charges_per_mass,
    )
    return build_fragment_yield(domain, charge_distribution, measured; total = total),
    symmetrize,
    width,
    inputs
end

# Whether a width file and a mean-TKE file come from one experiment, judged by the EXFOR entry
# that begins a retrieved file's name, `<accession>_<Author>_<year>.dat`. A width may be borrowed
# from another experiment on the same system; the run record says which it was.
function _width_provenance(width_file::AbstractString, energy_file::AbstractString)
    entry(name) = (m = match(r"^(\d{5})\d{3,}_", name); m === nothing ? nothing : m[1])
    width_entry, energy_entry = entry(width_file), entry(energy_file)
    (width_entry === nothing || energy_entry === nothing) && return "unverified"
    return width_entry == energy_entry ? "same experiment" : "borrowed"
end

# The heavy masses of the reconstructed yield that the width table reaches neither directly nor,
# when the identities are imposed, through the complement; those sit at their mean TKE. Taken
# from the masses the reconstruction placed, so a mass standing for an unmeasured complement
# counts and one without a mean TKE does not.
function _masses_without_width(
    reconstructed::MassEnergyYield,
    dispersion::KineticEnergyDispersion,
    compound_mass::Integer,
    symmetrize::Bool,
)
    dispersion.default === nothing || return Int[]
    widths = dispersion.σ_TKE
    return [
        A for A in reconstructed.masses if 2A >= compound_mass &&
            !haskey(widths, A) &&
            !(symmetrize && haskey(widths, compound_mass - A))
    ]
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
