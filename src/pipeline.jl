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
        result.configuration.subdirectory,
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
        configuration.polarization,
        configuration.heavy_masses;
        charges_per_mass = configuration.charges_per_mass,
    )
    weights = if configuration.yields === nothing
        weights_from(domain, configuration.polarization, configuration.kinetic_energies)
    else
        weights_from(configuration.yields)
    end
    reached =
        configuration.yields === nothing ? 1.0 :
        coverage(configuration.yields, configuration.kinetic_energies).reached

    ν = ConfigurationValues{Float64}()
    temperature = ConfigurationValues{Float64}()
    excitation = ConfigurationValues{Float64}()
    residual = ConfigurationValues{Float64}()
    neutron_energy = ConfigurationValues{Float64}()
    chains = Vector{
        Tuple{Fragmentation, Float64, EmissionChain, EmissionChain, ExcitationEnergies},
    }()

    attempted = 0
    solved = 0
    for entry in domain.entries, TKE in configuration.kinetic_energies
        attempted += 1
        balance = energetics(configuration.masses, configuration.system, entry, TKE)
        balance === nothing && continue
        split = partition(configuration.method, entry, value(balance.TXE))
        split === nothing && continue
        light = emission_chain(
            configuration.masses,
            configuration.model,
            configuration.cross_section,
            entry.light,
            split.light;
            maximum_temperature = configuration.maximum_temperature,
        )
        heavy = emission_chain(
            configuration.masses,
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
            "no configuration of $(configuration.subdirectory) could be solved; \
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
        ArgumentError("no configuration of $(configuration.subdirectory) carries both a \
                       multiplicity and a weight"),
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
    weights::ConfigurationValues{W},
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

function _mean_kinetic_energy(weights::ConfigurationValues{W}) where {W}
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
    write_results(result; directory) -> String

Write a run to `directory/<subdirectory>`: one two-column file per curve, a table of totals, and
a `metadata.toml` recording the configuration, the package version, the commit and the platform.

Nothing is overwritten in place. A directory that already holds a run is moved aside with a
numbered suffix first, so a second run of the same configuration cannot destroy the first.
"""
function write_results(
    result::RunResult;
    directory::AbstractString = joinpath(PACKAGE_ROOT, "results"),
)
    destination = joinpath(directory, result.configuration.subdirectory)
    _preserve(destination)
    mkpath(destination)

    digits = result.configuration.digits
    curves = (
        ("multiplicity_A.dat", "A nu", result.multiplicity),
        ("temperature_A.dat", "A T", result.temperature),
        ("excitation_A.dat", "A E_excitation", result.excitation),
        ("residual_energy_A.dat", "A E_residual", result.residual),
        ("neutron_energy_A.dat", "A epsilon", result.mean_neutron_energy),
        ("multiplicity_TKE.dat", "TKE nu", result.multiplicity_by_energy),
        ("excitation_TKE.dat", "TKE E_excitation", result.excitation_by_energy),
        ("spectrum_E.dat", "E N", result.spectrum),
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

    yields = result.configuration.yields
    if yields !== nothing
        for (name, dimension) in (
            ("yield_A.dat", :mass),
            ("yield_Z.dat", :charge),
            ("yield_TKE.dat", :kinetic_energy),
        )
            projection = yield_over(yields, dimension)
            open(joinpath(destination, name), "w") do io
                println(io, dimension === :kinetic_energy ? "TKE Y" : "$(dimension) Y")
                for (x, y) in projection
                    println(
                        io,
                        _significant(Float64(x), digits),
                        " ",
                        _significant(y, digits),
                    )
                end
            end
        end
    end

    open(joinpath(destination, "totals.dat"), "w") do io
        println(io, "quantity value")
        for (name, quantity) in (
            ("multiplicity_per_fragment", result.total_multiplicity),
            ("multiplicity_per_fission", 2 * result.total_multiplicity),
            ("mean_total_kinetic_energy", result.mean_kinetic_energy),
            ("mean_neutron_laboratory_energy", result.mean_neutron_lab_energy),
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
    println(io, "configuration = ", repr(configuration.source))
    println(io, "data_directory = ", repr(configuration.data_directory))
    println(io, "written = ", repr(string(Dates.now())))
    println(io, "package_version = ", repr(string(_package_version())))
    println(io, "commit = ", repr(_git_commit()))
    println(io)
    println(io, "[system]")
    println(io, "compound = ", repr(sprint(show, configuration.system.compound)))
    println(io, "target = ", repr(sprint(show, configuration.system.target)))
    println(io, "spontaneous = ", configuration.system.spontaneous)
    println(io, "incident_energy_mev = ", configuration.system.incident_energy)
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
        "kinetic_energies_mev = [",
        first(configuration.kinetic_energies),
        ", ",
        step(configuration.kinetic_energies),
        ", ",
        last(configuration.kinetic_energies),
        "]",
    )
    println(io, "level_density = ", repr(sprint(show, configuration.model)))
    println(io, "cross_section = ", repr(sprint(show, configuration.cross_section)))
    println(io, "maximum_temperature_mev = ", configuration.maximum_temperature)
    println(io, "partition = ", repr(sprint(show, configuration.method)))
    println(io, "charge_distribution = ", repr(configuration.polarization.source))
    yields = configuration.yields
    println(io, "yields = ", repr(yields === nothing ? "none" : sprint(show, yields)))
    println(io)
    println(io, "[result]")
    println(io, "multiplicity_per_fission = ", 2 * result.total_multiplicity)
    println(io, "mean_total_kinetic_energy_mev = ", result.mean_kinetic_energy)
    println(io, "mean_neutron_laboratory_energy_mev = ", result.mean_neutron_lab_energy)
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
