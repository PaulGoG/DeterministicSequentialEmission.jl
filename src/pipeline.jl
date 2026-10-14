"""
    Curve

A tabulated one-dimensional result, `Vector{Tuple{Float64,Float64}}` of abscissa and value
pairs, the content of every written two-column table.
"""
const Curve = Vector{Tuple{Float64, Float64}}

"""
    RunResult

What one run produced: the domain it swept, the averaged observables against fragment mass and
against total kinetic energy, the totals, the prompt neutron spectrum where one was asked for,
and the post-neutron product distributions where they were. [`run_pipeline`](@ref) produces it;
[`write_results`](@ref) writes it.

Curves are [`Curve`](@ref)s. The products, empty when `[products]` is disabled, are
`product_mass_yield` `(A_p, Y)`, `product_yield` `(Z, A_p, Y)` sorted by `Z` and then `A_p`,
`product_neutron_yield` `(N_p, Y)`, `product_kinetic_energy` `(A_p, E_K_p)`,
`product_kinetic_energy_by_charge` `(Z, E_K_p)` and `product_total_kinetic_energy`
`(A_p, TKE_p)` keyed by the heavy product mass; `mean_product_mass_loss`, `⟨A − A_p⟩` per
fragment, and `mean_product_total_kinetic_energy` are their totals, `NaN` when disabled.
`solved_mean_kinetic_energy` and `solved_kinetic_energy_dispersion` are the mean and the pooled
within-mass standard deviation of `TKE` over the configurations the sweep solved, the effective
`⟨TKE⟩` and `σ_TKE` of the run; `mean_kinetic_energy` is the mean of the input yield over every
cell it holds, solved or not. `gamma_energy` is the energy a fragment's neutrons leave,
[`gamma_energy`](@ref), against fragment mass; `gamma_energy_total_by_mass` and
`gamma_energy_total_by_energy` are the sum over the two fragments of a split against the mass of
either and against `TKE`, and `total_gamma_energy` its mean per fragment.
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
    product_mass_yield::Curve
    product_yield::Vector{Tuple{Int, Int, Float64}}
    product_neutron_yield::Curve
    product_kinetic_energy::Curve
    product_kinetic_energy_by_charge::Curve
    product_total_kinetic_energy::Curve
    mean_product_mass_loss::Float64
    mean_product_total_kinetic_energy::Float64
    solved_mean_kinetic_energy::Float64
    solved_kinetic_energy_dispersion::Float64
    gamma_energy::Curve
    gamma_energy_total_by_mass::Curve
    gamma_energy_total_by_energy::Curve
    total_gamma_energy::Float64

    # Explicit inner constructor: concrete argument types, no generic `convert` path into the
    # type.
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
        product_mass_yield::Curve,
        product_yield::Vector{Tuple{Int, Int, Float64}},
        product_neutron_yield::Curve,
        product_kinetic_energy::Curve,
        product_kinetic_energy_by_charge::Curve,
        product_total_kinetic_energy::Curve,
        mean_product_mass_loss::Real,
        mean_product_total_kinetic_energy::Real,
        solved_mean_kinetic_energy::Real,
        solved_kinetic_energy_dispersion::Real,
        gamma_energy::Curve,
        gamma_energy_total_by_mass::Curve,
        gamma_energy_total_by_energy::Curve,
        total_gamma_energy::Real,
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
            product_mass_yield,
            product_yield,
            product_neutron_yield,
            product_kinetic_energy,
            product_kinetic_energy_by_charge,
            product_total_kinetic_energy,
            Float64(mean_product_mass_loss),
            Float64(mean_product_total_kinetic_energy),
            Float64(solved_mean_kinetic_energy),
            Float64(solved_kinetic_energy_dispersion),
            gamma_energy,
            gamma_energy_total_by_mass,
            gamma_energy_total_by_energy,
            Float64(total_gamma_energy),
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
    run_pipeline(configuration) -> RunResult

Solve every configuration the domain and the kinetic energy grid span, contract the results
against the yield weights, and return the result without writing it; [`write_results`](@ref)
writes the tables, and `scripts/run.jl` adds the provenance. The sweep is deterministic, with
no sampling. A configuration is dropped where `TXE` is non-positive, the partition input does
not cover it, or the mass table does not reach a fragment; the counts of attempted and solved
configurations are reported with the result. The post-neutron product distributions are
contracted from the same chains, normalised to the total of the input yield, `[yield] total`,
or to `DEFAULT_YIELD_TOTAL` (200) where none was loaded.
"""
function run_pipeline(configuration::Configuration)
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
    γ = FragmentQuantity{Float64}()
    γ_pair = FragmentQuantity{Float64}()
    chains = Vector{
        Tuple{Fragmentation, Float64, EmissionChain, EmissionChain, ExcitationEnergies},
    }()

    attempted = 0
    solved = 0
    solved_keys = Set{Tuple{Int, Int, Float64}}()
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
        # The energy the two fragments' neutrons leave, booked against either fragment so that a
        # contraction over the mass of one gives the sum per split.
        pair_gamma = gamma_energy(light, split.light) + gamma_energy(heavy, split.heavy)
        for (fragment, energy, chain) in
            ((entry.light, split.light, light), (entry.heavy, split.heavy, heavy))
            key = (fragment.A, fragment.Z, Float64(TKE))
            push!(solved_keys, key)
            ν[key] = neutron_multiplicity(chain)
            γ[key] = gamma_energy(chain, energy)
            γ_pair[key] = pair_gamma
            excitation[key] = energy
            for (store, quantity, weighting) in (
                (temperature, :temperature, UniformSequenceWeighting()),
                (residual, :residual_energy, UniformSequenceWeighting()),
                (neutron_energy, :mean_neutron_energy, configuration.sequence_weighting),
            )
                averaged = sequence_average(chain, quantity, weighting)
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

    # Contract and narrow to `Float64` pairs in one step; the yield weights carry no uncertainty
    # here, and the narrowing keeps the result type concrete.
    curve(store) = _pairs(average(store, weights; over = (:charge, :kinetic_energy)))
    by_energy(store) = _pairs(average(store, weights; over = (:mass, :charge)))

    spectrum, mean_lab = _run_spectrum(configuration, chains, weights)
    covered = segment_coverage(configuration.method, configuration.heavy_masses).covered
    solved_mean, solved_dispersion = _solved_kinetic_energy(weights, solved_keys)
    # Guaranteed by `solved > 0` above; the check narrows the return type of the contraction.
    total = total_average(ν, weights)
    total === nothing && throw(
        ArgumentError(
            "no configuration of $(system_label(configuration.system)) carries both a \
                       multiplicity and a weight",
        ),
    )
    # Every fragment that carries a multiplicity carries a γ energy, so this contraction exists.
    total_γ = something(total_average(γ, weights))

    products = if configuration.products_enabled
        _product_distributions(
            chains,
            weights,
            configuration.system.compound.A,
            configuration.product_count,
            configuration.fragment_yield === nothing ? DEFAULT_YIELD_TOTAL :
            configuration.fragment_yield.total,
        )
    else
        (
            mass_yield = Curve(),
            yield = Tuple{Int, Int, Float64}[],
            neutron_yield = Curve(),
            kinetic_energy = Curve(),
            kinetic_energy_by_charge = Curve(),
            total_kinetic_energy = Curve(),
            mass_loss = NaN,
            mean_total_kinetic_energy = NaN,
        )
    end

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
        _pairs(
            average(
                neutron_energy,
                _neutron_weights(configuration.sequence_weighting, weights, ν);
                over = (:charge, :kinetic_energy),
            ),
        ),
        by_energy(ν),
        by_energy(excitation),
        spectrum,
        Float64(total),
        _mean_kinetic_energy(weights),
        mean_lab,
        products.mass_yield,
        products.yield,
        products.neutron_yield,
        products.kinetic_energy,
        products.kinetic_energy_by_charge,
        products.total_kinetic_energy,
        products.mass_loss,
        products.mean_total_kinetic_energy,
        solved_mean,
        solved_dispersion,
        curve(γ),
        curve(γ_pair),
        by_energy(γ_pair),
        Float64(total_γ),
    )
    return result
end

# Mean TKE and pooled within-mass standard deviation of TKE over the solved fragment cells,
# weighted as `_mean_kinetic_energy` weights every cell: the moments the run's averages carry.
function _solved_kinetic_energy(
    weights::FragmentQuantity{W},
    solved::Set{Tuple{Int, Int, Float64}},
) where {W}
    sums = Dict{Int, NTuple{3, Float64}}()
    for (key, w) in weights.values
        key in solved || continue
        weight = Float64(value(w))
        TKE = key[3]
        S₀, S₁, S₂ = get(sums, key[1], (0.0, 0.0, 0.0))
        sums[key[1]] = (S₀ + weight, S₁ + weight * TKE, S₂ + weight * TKE^2)
    end
    total = sum(first, values(sums); init = 0.0)
    total > 0 || return NaN, NaN
    mean = sum(S[2] for S in values(sums)) / total
    within = sum(S[3] - S[2]^2 / S[1] for S in values(sums) if S[1] > 0) / total
    return mean, sqrt(max(within, 0.0))
end

# Narrow a contraction to plain pairs. The return annotation fixes the type at the function
# boundary, so `RunResult` is built from concrete vectors.
function _pairs(contracted)::Curve
    out = Curve()
    for (x, y) in contracted
        push!(out, (Float64(x), Float64(value(y))))
    end
    return out
end

"""
    _pair_weight(fragmentation, yield) -> Float64

The weight of one fragmentation in a sum over pairs: the fragment yield of its heavy side times
`split_weight` of the domain package, which halves every mass-symmetric entry so that the
heavy-side sum is exactly half the sum over all fragments.
"""
_pair_weight(fragmentation::Fragmentation, yield::Float64) =
    yield * split_weight(fragmentation)

# Weights of the contraction of ⟨ε⟩ over charge and kinetic energy. Under the emission
# probability weighting the mean is taken over the emitted neutrons, so a fragment counts
# with its yield times its multiplicity; the uniform weighting keeps the yield alone.
_neutron_weights(
    ::UniformSequenceWeighting,
    weights::FragmentQuantity,
    ::FragmentQuantity,
) = weights

function _neutron_weights(
    ::EmissionProbabilityWeighting,
    weights::FragmentQuantity,
    ν::FragmentQuantity{Float64},
)
    scaled = FragmentQuantity{Float64}()
    for (key, weight) in weights.values
        multiplicity = get(ν.values, key, nothing)
        multiplicity === nothing && continue
        scaled[key] = Float64(value(weight)) * multiplicity
    end
    return scaled
end

# Laboratory spectrum of the system: each solved pair contributes its two chains, weighted by
# `_pair_weight` scaled by `total_weight`, and within the pair by `pair_weights`. Unit area by
# construction.
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
        w = _pair_weight(entry, Float64(value(weight)))
        w > 0 || continue
        share = pair_weights(configuration.spectrum_pair_weighting, light, heavy, split)
        share === nothing && continue
        w = total_weight(configuration.spectrum_average, w, light, heavy)
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
                    share[1] * chain_spectrum(
                        light,
                        configuration.cross_section,
                        E_light,
                        E,
                        configuration.sequence_weighting,
                    ) +
                    share[2] * chain_spectrum(
                        heavy,
                        configuration.cross_section,
                        E_heavy,
                        E,
                        configuration.sequence_weighting,
                    )
                )
        end
        return accumulated / total
    end

    values = Curve(undef, length(grid))
    for (index, E) in pairs(grid)
        values[index] = (E, spectrum(E))
    end
    # First moment over the tabulated grid, not the half-line, so it is reproducible from the
    # written file.
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

The token each configuration key contributes to a run identifier. This table is the single
source of abbreviations, so a directory name maps back to configuration keys.
[`run_identifier`](@ref) throws for a key with no entry. One entry, `kinetic_energy_bound`, is
not a key: it stands for `kinetic_energy_limit` and `kinetic_energy_moments` of `[yield]`
together, `"grid"` or the moments restored under the physical bound.
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
    "pair_weighting" => "pw",
    "average" => "avg",
    "sequence_weighting" => "sw",
    "count" => "pc",
    "total" => "Ytot",
    "symmetrize" => "Ysym",
    "kinetic_energy_bound" => "TKEb",
)

"""
    run_parameters(configuration) -> Dict{String, Any}

The configuration values that change the result, keyed by the token each contributes to the
run identifier.

Every key is abbreviated through [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref); a key with no entry
throws. The scalar and enumerated keys are covered, including `incident_energy`;
`deformed_branch` appears only under Gilbert-Cameron, and `total`, `symmetrize` and
`kinetic_energy_bound` only where a yield was loaded. `kinetic_energy_bound` is not a key of its
own: it is `"grid"` where `kinetic_energy_limit = "grid"`, and the `kinetic_energy_moments`
restored where the bound is physical, so the two keys take one token.
Not tokens, and recorded in the run metadata instead: the input files, the spectrum grid and
`maxwellian_temperature` (which affect the tabulation, and through the grid the first moment in
`totals.dat`, but not the solution), and `significant_digits`. The system names the directory
containing the identifier.
"""
function run_parameters(configuration::Configuration)
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
        "pair_weighting" => _pair_weighting_token(configuration.spectrum_pair_weighting),
        "average" => _average_token(configuration.spectrum_average),
        "sequence_weighting" => _sequence_weighting_token(configuration.sequence_weighting),
        "count" => _product_count_token(configuration.product_count),
    ]
    if configuration.model isa GilbertCameron
        push!(entries, "deformed_branch" => configuration.model.deformed_branch)
    end
    fragment_yield = configuration.fragment_yield
    if fragment_yield !== nothing
        push!(entries, "total" => fragment_yield.total)
        push!(entries, "symmetrize" => configuration.yield_symmetrized)
        push!(
            entries,
            "kinetic_energy_bound" => (
                configuration.yield_kinetic_energy_limit == "grid" ? "grid" :
                configuration.yield_kinetic_energy_moments
            ),
        )
    end

    parameters = Dict{String, Any}()
    for (key, value) in entries
        haskey(RUN_IDENTIFIER_ABBREVIATIONS, key) || throw(
            ArgumentError(
                "the run identifier uses the configuration key $(repr(key)), which has no \
                 entry in RUN_IDENTIFIER_ABBREVIATIONS",
            ),
        )
        parameters[RUN_IDENTIFIER_ABBREVIATIONS[key]] = value
    end
    return parameters
end

"""
    run_identifier(configuration) -> String

The token that names one run of one fissioning system: DrWatson's `savename` over
[`run_parameters`](@ref), tokens sorted and joined by `_`, every value rendered by one
rule — integers and integer-valued reals without a decimal part, everything else as
printed — so one grid has one spelling.

# Examples

```jldoctest
julia> RUN_IDENTIFIER_ABBREVIATIONS["heavy_mass_max"]
"AHmax"

julia> RUN_IDENTIFIER_ABBREVIATIONS["charges_per_mass"]
"nZ"
```
"""
function run_identifier(configuration::Configuration)
    return savename(run_parameters(configuration); connector = "_", val_to_string = _token)
end

_token(x::Integer) = string(x)
_token(x::Bool) = string(x)
_token(x::AbstractString) = String(x)
# Integer-valued reals are written without a decimal part, so one grid has one spelling.
_token(x::Real) = isinteger(x) ? string(Int(round(x))) : string(x)

_model_token(::BackShiftedFermiGas) = "BSFG"
_model_token(::GilbertCameron) = "GC"
_cross_section_token(::VariableCrossSection) = "variable"
_cross_section_token(::ConstantCrossSection) = "constant"
_method_token(::TemperatureRatio) = "temperature_ratio"
_method_token(::ExcitationEnergyRatio) = "excitation_ratio"
_method_token(::ScissionModelling) = "scission"
_pair_weighting_token(::MultiplicityWeighting) = "multiplicity"
_pair_weighting_token(::ExcitationWeighting) = "excitation"
_average_token(::PerFissionAverage) = "per_fission"
_average_token(::PerNeutronAverage) = "per_neutron"
_sequence_weighting_token(::UniformSequenceWeighting) = "uniform"
_sequence_weighting_token(::EmissionProbabilityWeighting) = "emission_probability"
_product_count_token(::EmissionProbabilityCount) = "emission_probability"
_product_count_token(::SequenceCount) = "sequence_count"

"""
    write_results(result, directory) -> String

Write the tables of a run into `directory`: one two-column file per curve, the fragment yield
projections, the post-neutron product distributions with the independent yields
`Y_vs_Z_A_p.dat` in three columns, the charge distribution the run used beside the ones it
could have used, and a table of totals.

`directory` is created. A directory already holding a run is refused, not overwritten;
`scripts/run.jl` moves an existing run aside as `#1`, `#2`, …. The provenance record is written
by the script, from [`run_metadata`](@ref) and its own additions.
"""
function write_results(result::RunResult, directory::AbstractString)
    isfile(joinpath(directory, "totals.dat")) && throw(
        ArgumentError(
            "$(directory) already holds a run; write_results never overwrites one, move it \
             aside first",
        ),
    )
    configuration = result.configuration
    mkpath(directory)

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
        ("Y_vs_A_p.dat", "A_p Y", result.product_mass_yield),
        ("Y_vs_N_p.dat", "N_p Y", result.product_neutron_yield),
        ("E_K_p_vs_A_p.dat", "A_p E_K_p", result.product_kinetic_energy),
        ("E_K_p_vs_Z.dat", "Z E_K_p", result.product_kinetic_energy_by_charge),
        ("TKE_p_vs_A_p.dat", "A_p TKE_p", result.product_total_kinetic_energy),
        ("E_gamma_vs_A.dat", "A E_gamma", result.gamma_energy),
        ("E_gamma_tot_vs_A.dat", "A E_gamma_tot", result.gamma_energy_total_by_mass),
        ("E_gamma_tot_vs_TKE.dat", "TKE E_gamma_tot", result.gamma_energy_total_by_energy),
    )
    # A yield placed at the mean TKE of each mass has no TKE distribution, so the tables resolved
    # in TKE would describe one it does not have; they are not written. Nor is the product TKE,
    # which a yield without a TKE distribution does not have either.
    distributed = configuration.yield_width != "none"
    resolved_in_energy = (
        "nu_vs_TKE.dat",
        "E_star_vs_TKE.dat",
        "Y_vs_TKE.dat",
        "TKE_p_vs_A_p.dat",
        "E_gamma_tot_vs_TKE.dat",
    )
    for (name, header, values) in curves
        isempty(values) && continue
        distributed || !(name in resolved_in_energy) || continue
        open(joinpath(directory, name), "w") do io
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
            distributed || !(name in resolved_in_energy) || continue
            projection = yield_over(fragment_yield, dimension)
            open(joinpath(directory, name), "w") do io
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

    if !isempty(result.product_yield)
        open(joinpath(directory, "Y_vs_Z_A_p.dat"), "w") do io
            println(io, "Z A_p Y")
            for (Z, A_p, Y) in result.product_yield
                println(
                    io,
                    _significant(Float64(Z), digits),
                    " ",
                    _significant(Float64(A_p), digits),
                    " ",
                    _significant(Y, digits),
                )
            end
        end
    end

    # Where the yield is bounded at the physical limit, how each mass came out of it: the
    # measured and retained moments, the Gaussian cut, and the share removed.
    bounds = configuration.yield_kinetic_energy_bounds
    if !isempty(bounds)
        open(joinpath(directory, "TKE_bound_vs_A.dat"), "w") do io
            println(
                io,
                "A TKE sigma_TKE TKE_measured sigma_TKE_measured TKE_parent sigma_TKE_parent removed",
            )
            for bound in bounds
                println(
                    io,
                    bound.A,
                    " ",
                    join(
                        (
                            _significant(x, digits) for x in (
                                bound.mean,
                                bound.dispersion,
                                bound.measured_mean,
                                bound.measured_dispersion,
                                bound.parent_mean,
                                bound.parent_dispersion,
                                bound.removed,
                            )
                        ),
                        " ",
                    ),
                )
            end
        end
    end

    _write_charge_distributions(directory, configuration, digits)

    totals = [
        ("nu_per_fragment", result.total_multiplicity),
        ("nu_per_fission", 2 * result.total_multiplicity),
        ("mean_TKE", result.mean_kinetic_energy),
        ("mean_E", result.mean_neutron_lab_energy),
        ("E_gamma_per_fragment", result.total_gamma_energy),
        ("E_gamma_per_fission", 2 * result.total_gamma_energy),
        ("configurations_solved", Float64(result.configurations_solved)),
        ("configurations_attempted", Float64(result.configurations_attempted)),
        ("yield_reached", result.yield_reached),
        ("segment_coverage", result.segment_coverage),
    ]
    if configuration.products_enabled
        push!(totals, ("mean_product_mass_loss", result.mean_product_mass_loss))
        push!(totals, ("mean_TKE_p", result.mean_product_total_kinetic_energy))
    end
    open(joinpath(directory, "totals.dat"), "w") do io
        println(io, "quantity value")
        for (name, quantity) in totals
            println(io, name, " ", _significant(quantity, digits))
        end
    end
    return directory
end

"""
    _write_charge_distributions(destination, configuration, digits)

Write the isobaric charge distribution the run used, beside the alternatives, as up to four
tables in the `A dZ sigma_Z` form of an evaluation. A `Zₚ` model is reduced at each heavy mass to
the plain Gaussian of the same first two moments, the width being the rms on the charge lattice.

- `charge_distribution_vs_A.dat` — the distribution `charge_distribution_file` selected.
- `charge_distribution_vs_A_wahl_per_reaction.dat` — the `Zₚ` model of Wahl, At. Data Nucl. Data
  Tables 39, 1 (1988), doi:10.1016/0092-640X(88)90016-2, with this reaction's own parameters,
  where that evaluation covers it.
- `charge_distribution_vs_A_wahl_systematics.dat` — the same model with parameters from eq. (17)
  of LA-13928 (2002), doi:10.2172/809574, where the system lies inside its fitted range.
- `charge_distribution_vs_A_mean.dat` — the conventional means, without mass dependence.

These are the layers of `charge_model`. Where the run used one of the alternatives, the
corresponding files coincide.
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
    reaction = wahl_1988(configuration.system)
    reaction === nothing || _write_charge_distribution(
        joinpath(destination, "charge_distribution_vs_A_wahl_per_reaction.dat"),
        reaction,
        masses,
        digits,
    )
    if is_wahl_applicable(configuration.mass_excess_table, configuration.system)
        systematics = WahlSystematics(configuration.mass_excess_table, configuration.system)
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

_write_charge_distribution(path::AbstractString, model::ZpModel, masses, digits::Integer) =
    _write_charge_distribution(
        path,
        effective_charge_distribution(model, collect(masses)),
        masses,
        digits,
    )

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

_significant(x::Real, digits::Integer) =
    isfinite(x) ? string(round(Float64(x); sigdigits = digits)) : string(Float64(x))

"""
    input_mean_kinetic_energy(fragment_yield) -> Float64

The mean total kinetic energy of the input yield, `Σ Y TKE / Σ Y` over every fragment it holds,
in MeV. It is a property of the measurement the yield came from, compared in the run record with
[`recommended_mean_total_kinetic_energy`](@ref) where a standard exists.
"""
function input_mean_kinetic_energy(fragment_yield::FragmentYield)
    weighted = sum(Y * key[3] for (key, Y) in fragment_yield.values)
    return weighted / sum(values(fragment_yield.values))
end

"""
    run_metadata(result) -> Dict{String, Any}

What the library knows about a run, as nested tables ready for `TOML.print`: `"system"`, the
fissioning system; `"model"`, every modelling choice and input the configuration resolved;
`"result"`, the totals. The caller adds what only it knows — the run identifier, the paths,
the commit, the platform — and writes the whole; `scripts/run.jl` does.
"""
function run_metadata(result::RunResult)
    configuration = result.configuration
    grid = configuration.spectrum_grid
    fragment_yield = configuration.fragment_yield
    system = system_record(configuration.system)
    model = Dict{String, Any}(
        "heavy_masses" =>
            [first(configuration.heavy_masses), last(configuration.heavy_masses)],
        "charges_per_mass" => configuration.charges_per_mass,
        "kinetic_energies_MeV" => [
            first(configuration.kinetic_energies),
            step(configuration.kinetic_energies),
            last(configuration.kinetic_energies),
        ],
        "level_density" => sprint(show, configuration.model),
        "cross_section" => sprint(show, configuration.cross_section),
        "maximum_temperature_MeV" => configuration.maximum_temperature,
        "partition" => sprint(show, configuration.method),
        "partition_curve" => _curve_source(configuration.method),
        "charge_distribution" => charge_model_label(configuration.charge_distribution),
        "mass_excess" => configuration.mass_excess_table.source,
        "significant_digits" => configuration.significant_digits,
        "maxwellian_temperature_MeV" => configuration.maxwellian_temperature,
        "spectrum_pair_weighting" =>
            _pair_weighting_token(configuration.spectrum_pair_weighting),
        "spectrum_average" => _average_token(configuration.spectrum_average),
        "spectrum_sequence_weighting" =>
            _sequence_weighting_token(configuration.sequence_weighting),
        "spectrum_grid_MeV" => grid === nothing ? "none" : [first(grid), last(grid)],
        "spectrum_points" => grid === nothing ? 0 : length(grid),
        "products" => configuration.products_enabled,
        "product_count" => _product_count_token(configuration.product_count),
        "yield" => fragment_yield === nothing ? "none" : sprint(show, fragment_yield),
        "yield_symmetrized" => configuration.yield_symmetrized,
        "yield_kinetic_energy_width" => configuration.yield_width,
        "yield_inputs" => configuration.yield_inputs,
        "yield_kinetic_energy_limit" => configuration.yield_kinetic_energy_limit,
    )
    configuration.yield_kinetic_energy_limit == "grid" ||
        (model["yield_kinetic_energy_moments"] = configuration.yield_kinetic_energy_moments)
    bounds = configuration.yield_kinetic_energy_bounds
    if fragment_yield !== nothing && !isempty(bounds)
        model["kinetic_energy_bound"] = _bound_summary(bounds, fragment_yield)
    end
    # How far the input yield's mean TKE lies from the energy standard, where one exists: a
    # measured scale offset of 1 MeV moves ⟨ν⟩ by about 0.08.
    if fragment_yield !== nothing
        input = input_mean_kinetic_energy(fragment_yield)
        model["input_mean_TKE_MeV"] = input
        standard = recommended_mean_total_kinetic_energy(configuration.system)
        if standard !== nothing
            model["recommended_mean_TKE_MeV"] = [value(standard), uncertainty(standard)]
            model["input_mean_TKE_deviation_MeV"] = input - value(standard)
        end
    end
    totals = Dict{String, Any}(
        "nu_per_fission" => 2 * result.total_multiplicity,
        "mean_TKE_MeV" => result.mean_kinetic_energy,
        "mean_E_MeV" => result.mean_neutron_lab_energy,
        "E_gamma_per_fission_MeV" => 2 * result.total_gamma_energy,
        "configurations_solved" => result.configurations_solved,
        "configurations_attempted" => result.configurations_attempted,
        "yield_reached" => result.yield_reached,
        "segment_coverage" => result.segment_coverage,
        "solved_mean_TKE_MeV" => result.solved_mean_kinetic_energy,
        "solved_sigma_TKE_MeV" => result.solved_kinetic_energy_dispersion,
    )
    if configuration.products_enabled
        totals["mean_product_mass_loss"] = result.mean_product_mass_loss
        totals["mean_TKE_p_MeV"] = result.mean_product_total_kinetic_energy
    end
    return Dict{String, Any}("system" => system, "model" => model, "result" => totals)
end

# The bound of a yield summarised over masses with the yield each carries: the measured and the
# retained mean TKE, the pooled measured and retained widths (root mean square), the largest
# per-mass miss of the mean, the share removed, and the mean shift and width ratio of the
# Gaussians cut.
function _bound_summary(bounds::Vector{KineticEnergyBound}, fragment_yield::FragmentYield)
    by_mass = Dict(Int(A) => Y for (A, Y) in yield_over(fragment_yield, :mass))
    weighted(f, keep = _ -> true) = (
        numerator = sum((by_mass[b.A] * f(b) for b in bounds if keep(b)); init = 0.0);
        denominator = sum((by_mass[b.A] for b in bounds if keep(b)); init = 0.0);
        denominator > 0 ? numerator / denominator : NaN
    )
    widths(b) = isfinite(b.measured_dispersion)
    gaussian(b) = isfinite(b.parent_mean) && isfinite(b.parent_dispersion)
    return Dict{String, Any}(
        "measured_mean_TKE_MeV" => weighted(b -> b.measured_mean),
        "mean_TKE_MeV" => weighted(b -> b.mean),
        "measured_sigma_TKE_MeV" => sqrt(weighted(b -> b.measured_dispersion^2, widths)),
        "sigma_TKE_MeV" => sqrt(weighted(b -> b.dispersion^2, widths)),
        "max_mean_TKE_deviation_by_mass_MeV" =>
            maximum(abs(b.mean - b.measured_mean) for b in bounds),
        "removed_share" => weighted(b -> b.removed),
        "parent_mean_shift_MeV" => weighted(b -> b.parent_mean - b.measured_mean, gaussian),
        "parent_sigma_ratio" =>
            weighted(b -> b.parent_dispersion / b.measured_dispersion, gaussian),
        "masses" => length(bounds),
    )
end

_curve_source(method::TemperatureRatio) =
    isempty(method.curve.source) ? "supplied in memory" : method.curve.source
_curve_source(method::ExcitationEnergyRatio) =
    isempty(method.curve.source) ? "supplied in memory" : method.curve.source
_curve_source(::ScissionModelling) = "none; the partition is modelled at scission"
