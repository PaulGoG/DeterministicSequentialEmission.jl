"""
    Curve

A tabulated one-dimensional result, `Vector{Tuple{Float64,Float64}}`. Abscissa and value, which
is what every written table holds and what reads back without either side having to know the
abscissa in advance.
"""
const Curve = Vector{Tuple{Float64, Float64}}

"""
    RunResult

What one run produced: the domain it swept, the averaged observables against fragment mass and
against total kinetic energy, the totals, and the prompt neutron spectrum where one was asked
for. [`write_results`](@ref) writes it; [`run_pipeline`](@ref) does both.

Curves are [`Curve`](@ref)s.
"""
struct RunResult
    configuration::Configuration
    configurations_solved::Int
    configurations_attempted::Int
    yield_reached::Float64
    segment_coverage::Float64
    multiplicity::Curve
    temperature::Curve
    excitation::Curve
    residual::Curve
    mean_neutron_energy::Curve
    multiplicity_by_energy::Curve
    excitation_by_energy::Curve
    spectrum::Curve
    total_multiplicity::Float64
    mean_kinetic_energy::Float64
    mean_neutron_lab_energy::Float64

    # Written out rather than left to the default constructor: with the fields declared and no
    # generic path into the type, every argument is concrete at the call site and nothing has to
    # be converted from an abstract array on the way in.
    function RunResult(
        configuration::Configuration,
        configurations_solved::Integer,
        configurations_attempted::Integer,
        yield_reached::Real,
        segment_coverage::Real,
        multiplicity::Curve,
        temperature::Curve,
        excitation::Curve,
        residual::Curve,
        mean_neutron_energy::Curve,
        multiplicity_by_energy::Curve,
        excitation_by_energy::Curve,
        spectrum::Curve,
        total_multiplicity::Real,
        mean_kinetic_energy::Real,
        mean_neutron_lab_energy::Real,
    )
        return new(
            configuration,
            Int(configurations_solved),
            Int(configurations_attempted),
            Float64(yield_reached),
            Float64(segment_coverage),
            multiplicity,
            temperature,
            excitation,
            residual,
            mean_neutron_energy,
            multiplicity_by_energy,
            excitation_by_energy,
            spectrum,
            Float64(total_multiplicity),
            Float64(mean_kinetic_energy),
            Float64(mean_neutron_lab_energy),
        )
    end
end

function Base.show(io::IO, result::RunResult)
    print(
        io,
        "RunResult(",
        system_label(result.configuration.system),
        ", ⟨ν⟩ = ",
        round(2 * result.total_multiplicity; digits = 4),
        " per fission, ",
        result.configurations_solved,
        " configurations)",
    )
    return nothing
end

"""
    run_pipeline(configuration; write = true) -> RunResult

Solve every configuration the domain and the kinetic energy grid span, contract the results
against the yield weights, and write the tables.

The sweep is deterministic: no sampling, no seeds. A configuration is dropped only where the
physics excludes it — a non-positive `TXE`, a partition the segments do not cover, a fragment
the mass table cannot reach — and the counts of attempted and solved configurations are
reported so that the coverage of a run can be audited from its own output.
"""
function run_pipeline(configuration::Configuration; write::Bool = true)
    domain = fragmentation_domain(
        configuration.system,
        configuration.charge_distribution,
        configuration.heavy_masses;
        charges_per_mass = configuration.charges_per_mass,
    )
    weights = if configuration.fragment_yield === nothing
        weights_from(domain, configuration.charge_distribution, configuration.kinetic_energies)
    else
        weights_from(configuration.fragment_yield)
    end
    reached =
        configuration.fragment_yield === nothing ? 1.0 :
        coverage(configuration.fragment_yield, configuration.kinetic_energies).reached

    ν = FragmentQuantity{Float64}()
    temperature = FragmentQuantity{Float64}()
    excitation = FragmentQuantity{Float64}()
    residual = FragmentQuantity{Float64}()
    neutron_energy = FragmentQuantity{Float64}()
    chains = Vector{
        Tuple{Fragmentation, Float64, EmissionChain, EmissionChain, ExcitationEnergies},
    }()

    attempted = 0
    solved = 0
    for entry in domain.entries, TKE in configuration.kinetic_energies
        attempted += 1
        balance =
            energetics(configuration.mass_excess_table, configuration.system, entry, TKE)
        balance === nothing && continue
        split = partition(configuration.method, entry, value(balance.TXE))
        split === nothing && continue
        light = emission_chain(
            configuration.mass_excess_table,
            configuration.model,
            configuration.cross_section,
            entry.light,
            split.light;
            maximum_temperature = configuration.maximum_temperature,
        )
        heavy = emission_chain(
            configuration.mass_excess_table,
            configuration.model,
            configuration.cross_section,
            entry.heavy,
            split.heavy;
            maximum_temperature = configuration.maximum_temperature,
        )
        solved += 1
        for (fragment, energy, chain) in
            ((entry.light, split.light, light), (entry.heavy, split.heavy, heavy))
            key = (fragment.A, fragment.Z, Float64(TKE))
            ν[key] = neutron_multiplicity(chain)
            excitation[key] = energy
            for (store, quantity) in (
                (temperature, :temperature),
                (residual, :residual_energy),
                (neutron_energy, :mean_neutron_energy),
            )
                averaged = sequence_average(chain, quantity)
                averaged === nothing || (store[key] = averaged)
            end
        end
        push!(chains, (entry, Float64(TKE), light, heavy, split))
    end

    solved > 0 || throw(
        ArgumentError(
            "no configuration of $(system_label(configuration.system)) could be solved; \
                       check the kinetic energy grid against the TXE the system has to share",
        ),
    )

    # Contracted and narrowed to plain pairs in one step. The weights of a run are yields, which
    # carry no uncertainty of their own here, so the result is `Float64`; saying so keeps the
    # result type concrete rather than leaving it to inference across the contraction.
    curve(store) = _pairs(average(store, weights; over = (:charge, :kinetic_energy)))
    by_energy(store) = _pairs(average(store, weights; over = (:mass, :charge)))

    spectrum, mean_lab = _run_spectrum(configuration, chains, weights)
    covered = segment_coverage(configuration.method, configuration.heavy_masses).covered
    # `solved > 0` above guarantees this, but the contraction cannot say so in its return type.
    total = total_average(ν, weights)
    total === nothing && throw(
        ArgumentError(
            "no configuration of $(system_label(configuration.system)) carries both a \
                       multiplicity and a weight",
        ),
    )

    result = RunResult(
        configuration,
        solved,
        attempted,
        reached,
        covered,
        curve(ν),
        curve(temperature),
        curve(excitation),
        curve(residual),
        curve(neutron_energy),
        by_energy(ν),
        by_energy(excitation),
        spectrum,
        Float64(total),
        _mean_kinetic_energy(weights),
        mean_lab,
    )
    write && write_results(result)
    return result
end

# Narrow a contraction to plain pairs. The return annotation is load-bearing: it fixes the type
# at the function boundary, so the result struct is built from concrete vectors rather than from
# whatever inference makes of a comprehension over a parametric result.
function _pairs(contracted)::Curve
    out = Curve()
    for (x, y) in contracted
        push!(out, (Float64(x), Float64(value(y))))
    end
    return out
end

# The laboratory spectrum of the whole system: every solved pair contributes its two chains,
# weighted by the yield of the fragmentation and, within a pair, by how the neutrons divide
# between the fragments. Normalised to unit area by construction, so the ratio to a Maxwellian
# is an absolute statement.
function _run_spectrum(
    configuration::Configuration,
    chains,
    weights::FragmentQuantity{W},
) where {W}
    grid = configuration.spectrum_grid
    grid === nothing && return Curve(), NaN

    records = Tuple{
        Float64,
        EmissionChain,
        EmissionChain,
        Tuple{Float64, Float64},
        Float64,
        Float64,
    }[]
    total = 0.0
    for (entry, TKE, light, heavy, split) in chains
        weight = get(weights.values, (entry.heavy.A, entry.heavy.Z, TKE), nothing)
        weight === nothing && continue
        w = Float64(value(weight))
        w > 0 || continue
        share = pair_weights(MultiplicityWeighting(), light, heavy, split)
        share === nothing && continue
        E_light, E_heavy = fragment_energy_per_nucleon(entry, TKE)
        push!(records, (w, light, heavy, share, E_light, E_heavy))
        total += w
    end
    isempty(records) && return Curve(), NaN

    function spectrum(E)
        accumulated = 0.0
        for (w, light, heavy, share, E_light, E_heavy) in records
            accumulated +=
                w * (
                    share[1] *
                    chain_spectrum(light, configuration.cross_section, E_light, E) +
                    share[2] *
                    chain_spectrum(heavy, configuration.cross_section, E_heavy, E)
                )
        end
        return accumulated / total
    end

    values = Curve(undef, length(grid))
    for (index, E) in pairs(grid)
        values[index] = (E, spectrum(E))
    end
    # The first moment over the tabulated grid rather than over the whole half-line: what is
    # reported is a property of what was written, so a reader can reproduce it from the file.
    norm = _trapezoid(values, (E, N) -> N)
    mean = norm > 0 ? _trapezoid(values, (E, N) -> E * N) / norm : NaN
    return values, mean
end

function _trapezoid(values::Curve, integrand)
    length(values) < 2 && return 0.0
    accumulated = 0.0
    for index in 2:length(values)
        E₀, N₀ = values[index - 1]
        E₁, N₁ = values[index]
        accumulated += (E₁ - E₀) * (integrand(E₀, N₀) + integrand(E₁, N₁)) / 2
    end
    return accumulated
end

function _mean_kinetic_energy(weights::FragmentQuantity{W}) where {W}
    weighted = 0.0
    total = 0.0
    for ((_, _, TKE), w) in weights.values
        weight = Float64(value(w))
        weighted += TKE * weight
        total += weight
    end
    return total > 0 ? weighted / total : NaN
end

"""
    RUN_IDENTIFIER_ABBREVIATIONS

The token each configuration key contributes to a run identifier.

A token *is* the key it came from wherever that is short enough to carry; where it is not, the
abbreviation is drawn from here and from nowhere else, so there is one place to look and a
reader can get from a directory name back to a configuration without reading the source.
[`run_identifier`](@ref) refuses a key with no entry rather than inventing a token.
"""
const RUN_IDENTIFIER_ABBREVIATIONS = Dict(
    "incident_energy" => "En",
    "heavy_mass_min" => "AHmin",
    "heavy_mass_max" => "AHmax",
    "charges_per_mass" => "nZ",
    "kinetic_energy_min" => "TKEmin",
    "kinetic_energy_max" => "TKEmax",
    "kinetic_energy_step" => "TKEstep",
    "model" => "ldm",
    "deformed_branch" => "def",
    "cross_section" => "sigma",
    "maximum_temperature" => "Tmax",
    "method" => "part",
    "total" => "Ytot",
)

"""
    run_identifier(configuration) -> String

The token that names one run of one fissioning system, built from the configuration keys that
change the result.

Every token is a configuration key abbreviated through [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref),
and a key with no entry there throws rather than being given a token on the spot.

**It covers the scalar and enumerated keys, and it says which rather than claiming all of them.**
`deformed_branch` appears only under Gilbert-Cameron and `total` only where a yield was loaded,
because those are the runs each has an effect on. `incident_energy` is covered: it changes the
compound excitation and so every quantity downstream of it, and two runs of one system on one
channel may differ in nothing else.

Four things are deliberately **not** tokens and are written in full into the metadata beside the
result instead. The input files — the mass evaluation, the charge distribution, the segmented
curve, the yields — change the result but do not reduce to a savename token. The spectrum grid
and `maxwellian_temperature` change the tabulation of an already-solved spectrum rather than the
physics, though the grid does move the first moment reported in `totals.dat`, which is a property
of what was written. `significant_digits` changes how a number is rendered. And the system is not
a token either: it names the directory the identifier sits in.

# Examples

```jldoctest
julia> RUN_IDENTIFIER_ABBREVIATIONS["heavy_mass_max"]
"AHmax"

julia> RUN_IDENTIFIER_ABBREVIATIONS["charges_per_mass"]
"nZ"
```
"""
function run_identifier(configuration::Configuration)
    entries = Pair{String, Any}[
        "incident_energy" => configuration.system.incident_energy,
        "heavy_mass_min" => first(configuration.heavy_masses),
        "heavy_mass_max" => last(configuration.heavy_masses),
        "charges_per_mass" => configuration.charges_per_mass,
        "kinetic_energy_min" => first(configuration.kinetic_energies),
        "kinetic_energy_max" => last(configuration.kinetic_energies),
        "kinetic_energy_step" => step(configuration.kinetic_energies),
        "model" => _model_token(configuration.model),
        "cross_section" => _cross_section_token(configuration.cross_section),
        "maximum_temperature" => configuration.maximum_temperature,
        "method" => _method_token(configuration.method),
    ]
    if configuration.model isa GilbertCameron
        push!(entries, "deformed_branch" => configuration.model.deformed_branch)
    end
    fragment_yield = configuration.fragment_yield
    fragment_yield === nothing || push!(entries, "total" => fragment_yield.total)

    tokens = String[]
    for (key, value) in entries
        haskey(RUN_IDENTIFIER_ABBREVIATIONS, key) || throw(
            ArgumentError(
                "the run identifier uses the configuration key $(repr(key)), which has no \
                 entry in RUN_IDENTIFIER_ABBREVIATIONS",
            ),
        )
        push!(tokens, string(RUN_IDENTIFIER_ABBREVIATIONS[key], "=", _token(value)))
    end
    sort!(tokens)
    return join(tokens, "_")
end

_token(x::Integer) = string(x)
_token(x::Bool) = string(x)
_token(x::AbstractString) = String(x)
# Trailing zeros carry no information and would make two spellings of one grid.
_token(x::Real) = isinteger(x) ? string(Int(round(x))) : string(x)

_model_token(::BackShiftedFermiGas) = "BSFG"
_model_token(::GilbertCameron) = "GC"
_cross_section_token(::VariableCrossSection) = "variable"
_cross_section_token(::ConstantCrossSection) = "constant"
_method_token(::TemperatureRatio) = "temperature_ratio"
_method_token(::ExcitationEnergyRatio) = "excitation_ratio"
_method_token(::ScissionModelling) = "scission"

"""
    write_results(result; directory) -> String

Write a run to `directory/<system>/<run>`: one two-column file per curve, a table of totals, and
a `metadata.toml` recording the configuration, the package version, the commit and the platform.

The system is [`system_label`](@ref) and the run is [`run_identifier`](@ref), so two runs of one
system that differ in anything the identifier covers land in different directories rather than
one displacing the other, and a directory name states what produced it.

Nothing is overwritten in place. A directory that already holds a run is moved aside with a
numbered suffix first, so a second run of the same configuration cannot destroy the first.
"""
function write_results(
    result::RunResult;
    directory::AbstractString = joinpath(PACKAGE_ROOT, "results"),
)
    configuration = result.configuration
    destination = joinpath(
        directory,
        system_label(configuration.system),
        run_identifier(configuration),
    )
    _preserve(destination)
    mkpath(destination)

    digits = configuration.significant_digits
    curves = (
        ("nu_vs_A.dat", "A nu", result.multiplicity),
        ("T_vs_A.dat", "A T", result.temperature),
        ("E_star_vs_A.dat", "A E_star", result.excitation),
        ("E_r_vs_A.dat", "A E_r", result.residual),
        ("eps_vs_A.dat", "A eps", result.mean_neutron_energy),
        ("nu_vs_TKE.dat", "TKE nu", result.multiplicity_by_energy),
        ("E_star_vs_TKE.dat", "TKE E_star", result.excitation_by_energy),
        ("spectrum_vs_E.dat", "E spectrum", result.spectrum),
    )
    for (name, header, values) in curves
        isempty(values) && continue
        open(joinpath(destination, name), "w") do io
            println(io, header)
            for (x, y) in values
                println(io, _significant(x, digits), " ", _significant(y, digits))
            end
        end
    end

    fragment_yield = configuration.fragment_yield
    if fragment_yield !== nothing
        for (name, header, dimension) in (
            ("Y_vs_A.dat", "A Y", :mass),
            ("Y_vs_Z.dat", "Z Y", :charge),
            ("Y_vs_TKE.dat", "TKE Y", :kinetic_energy),
        )
            projection = yield_over(fragment_yield, dimension)
            open(joinpath(destination, name), "w") do io
                println(io, header)
                for (x, Y) in projection
                    println(
                        io,
                        _significant(Float64(x), digits),
                        " ",
                        _significant(Y, digits),
                    )
                end
            end
        end
    end

    _write_charge_distributions(destination, configuration, digits)

    open(joinpath(destination, "totals.dat"), "w") do io
        println(io, "quantity value")
        for (name, quantity) in (
            ("nu_per_fragment", result.total_multiplicity),
            ("nu_per_fission", 2 * result.total_multiplicity),
            ("mean_TKE", result.mean_kinetic_energy),
            ("mean_E", result.mean_neutron_lab_energy),
            ("configurations_solved", Float64(result.configurations_solved)),
            ("configurations_attempted", Float64(result.configurations_attempted)),
            ("yield_reached", result.yield_reached),
            ("segment_coverage", result.segment_coverage),
        )
            println(io, name, " ", _significant(quantity, digits))
        end
    end

    open(joinpath(destination, "metadata.toml"), "w") do io
        _write_metadata(io, result)
    end
    return destination
end

"""
    _write_charge_distributions(destination, configuration, digits)

Write the isobaric charge distribution the run used, beside the two it could have used instead.

Which `ΔZ(A)`, `σ_Z(A)` a run was weighted by is a choice with a measurable effect on `ν(A)`, and
until it is written down the only record of it is a path in a configuration. Three tables, each in
the same `A dZ sigma_Z` form an evaluation is supplied in, so the comparison needs nothing but a
reader:

- `charge_distribution_vs_A.dat` — what the run used, whatever `charge_distribution_file` selected.
- `charge_distribution_vs_A_wahl_systematics.dat` — the Wahl `Zₚ` systematics over the same masses,
  where [`is_wahl_applicable`](@ref) says the system is inside the fitted range. Marked as
  generated in its name, as a systematics table always is.
- `charge_distribution_vs_A_mean.dat` — the conventional means, the fallback that carries no mass
  dependence at all.

Where the run already took the systematics, the first two coincide, and that is the honest
picture rather than a redundancy to suppress.
"""
function _write_charge_distributions(
    destination::AbstractString,
    configuration::Configuration,
    digits::Integer,
)
    masses = configuration.heavy_masses
    _write_charge_distribution(
        joinpath(destination, "charge_distribution_vs_A.dat"),
        configuration.charge_distribution,
        masses,
        digits,
    )
    if is_wahl_applicable(configuration.mass_excess_table, configuration.system)
        systematics = wahl_charge_distribution(
            configuration.mass_excess_table,
            configuration.system,
            masses,
        )
        systematics === nothing || _write_charge_distribution(
            joinpath(destination, "charge_distribution_vs_A_wahl_systematics.dat"),
            systematics,
            masses,
            digits,
        )
    end
    _write_charge_distribution(
        joinpath(destination, "charge_distribution_vs_A_mean.dat"),
        mean_charge_distribution(),
        masses,
        digits,
    )
    return nothing
end

function _write_charge_distribution(
    path::AbstractString,
    distribution::ChargeDistribution,
    masses,
    digits::Integer,
)
    open(path, "w") do io
        println(io, "A dZ sigma_Z")
        for A in masses
            println(
                io,
                A,
                " ",
                _significant(charge_polarization(distribution, A), digits),
                " ",
                _significant(charge_dispersion(distribution, A), digits),
            )
        end
    end
    return path
end

# Never overwrite a run. The previous one moves to <name>#1, <name>#2, … in the manner of
# DrWatson's safesave, without taking the dependency for one function.
function _preserve(destination::AbstractString)
    isdir(destination) || return nothing
    index = 1
    while ispath("$(destination)#$(index)")
        index += 1
    end
    mv(destination, "$(destination)#$(index)")
    return nothing
end

_significant(x::Real, digits::Integer) =
    isfinite(x) ? string(round(Float64(x); sigdigits = digits)) : string(Float64(x))

function _write_metadata(io::IO, result::RunResult)
    configuration = result.configuration
    println(
        io,
        "# Written by DeterministicSequentialEmission.jl. Every field records what the",
    )
    println(
        io,
        "# run actually used, so a result can be traced to a configuration, a commit",
    )
    println(io, "# and a machine without consulting anything outside this directory.")
    println(io)
    println(io, "[run]")
    println(io, "identifier = ", repr(run_identifier(configuration)))
    println(io, "configuration = ", repr(configuration.source))
    println(io, "data_directory = ", repr(configuration.data_directory))
    println(io, "written = ", repr(string(Dates.now())))
    println(io, "package_version = ", repr(string(_package_version())))
    println(io, "commit = ", repr(_git_commit()))
    println(io)
    println(io, "[system]")
    println(io, "label = ", repr(system_label(configuration.system)))
    println(io, "notation = ", repr(system_notation(configuration.system)))
    println(io, "target_A = ", configuration.system.target.A)
    println(io, "target_Z = ", configuration.system.target.Z)
    println(io, "channel = ", repr(configuration.system.channel))
    println(io, "reaction = ", repr(reaction(configuration.system)))
    println(io, "incident_energy_MeV = ", configuration.system.incident_energy)
    println(io, "compound_A = ", configuration.system.compound.A)
    println(io, "compound_Z = ", configuration.system.compound.Z)
    println(io)
    println(io, "[model]")
    println(
        io,
        "heavy_masses = [",
        first(configuration.heavy_masses),
        ", ",
        last(configuration.heavy_masses),
        "]",
    )
    println(io, "charges_per_mass = ", configuration.charges_per_mass)
    println(
        io,
        "kinetic_energies_MeV = [",
        first(configuration.kinetic_energies),
        ", ",
        step(configuration.kinetic_energies),
        ", ",
        last(configuration.kinetic_energies),
        "]",
    )
    println(io, "level_density = ", repr(sprint(show, configuration.model)))
    println(io, "cross_section = ", repr(sprint(show, configuration.cross_section)))
    println(io, "maximum_temperature_MeV = ", configuration.maximum_temperature)
    println(io, "partition = ", repr(sprint(show, configuration.method)))
    # Which tabulation the partition actually used. A curve staged by hand used to leave no
    # trace of the run that produced it; the file name a producing run writes carries its own
    # run identifier, so recording the resolved path here closes that gap.
    println(io, "partition_curve = ", repr(_curve_source(configuration.method)))
    println(io, "charge_distribution = ", repr(configuration.charge_distribution.source))
    println(io, "mass_excess = ", repr(configuration.mass_excess_table.source))
    println(io, "significant_digits = ", configuration.significant_digits)
    println(io, "maxwellian_temperature_MeV = ", configuration.maxwellian_temperature)
    grid = configuration.spectrum_grid
    println(
        io,
        "spectrum_grid_MeV = ",
        grid === nothing ? repr("none") :
        string("[", first(grid), ", ", last(grid), ", ", length(grid), "]"),
    )
    fragment_yield = configuration.fragment_yield
    println(
        io,
        "yield = ",
        repr(fragment_yield === nothing ? "none" : sprint(show, fragment_yield)),
    )
    println(io)
    println(io, "[result]")
    println(io, "nu_per_fission = ", 2 * result.total_multiplicity)
    println(io, "mean_TKE_MeV = ", result.mean_kinetic_energy)
    println(io, "mean_E_MeV = ", result.mean_neutron_lab_energy)
    println(io, "configurations_solved = ", result.configurations_solved)
    println(io, "configurations_attempted = ", result.configurations_attempted)
    println(io, "yield_reached = ", result.yield_reached)
    println(io, "segment_coverage = ", result.segment_coverage)
    println(io)
    println(io, "[platform]")
    println(io, "julia = ", repr(string(VERSION)))
    println(io, "hostname = ", repr(gethostname()))
    println(io, "cpu = ", repr(_cpu_model()))
    println(io, "cpu_threads = ", Sys.CPU_THREADS)
    println(io, "threads = ", Threads.nthreads())
    println(io, "total_memory_gb = ", round(Sys.total_memory() / 2^30; digits = 2))
    return nothing
end

_curve_source(method::TemperatureRatio) =
    isempty(method.curve.source) ? "supplied in memory" : method.curve.source
_curve_source(method::ExcitationEnergyRatio) =
    isempty(method.curve.source) ? "supplied in memory" : method.curve.source
_curve_source(::ScissionModelling) = "none; the partition is modelled at scission"

function _package_version()
    project = joinpath(PACKAGE_ROOT, "Project.toml")
    isfile(project) || return v"0.0.0"
    document = TOML.parsefile(project)
    return VersionNumber(get(document, "version", "0.0.0"))
end

function _git_commit()
    repository = PACKAGE_ROOT
    isdir(joinpath(repository, ".git")) || return "not a git checkout"
    try
        commit = readchomp(`git -C $repository rev-parse HEAD`)
        dirty = !isempty(readchomp(`git -C $repository status --porcelain`))
        return dirty ? commit * "-dirty" : commit
    catch
        return "unavailable"
    end
end

function _cpu_model()
    info = Sys.cpu_info()
    isempty(info) && return "unknown"
    return String(first(info).model)
end
