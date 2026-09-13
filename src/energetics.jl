"""
    compound_nucleus_excitation(table, system) -> Union{Measurement{Float64},Nothing}

Excitation energy of the compound nucleus in MeV.

Zero for spontaneous fission. For neutron-induced fission it is the neutron separation energy
of the compound nucleus plus the incident kinetic energy, `E*_CN = Sₙ(A₀,Z₀) + Eₙ`. Returns
`nothing` when `Sₙ` cannot be formed from the mass table.
"""
function compound_nucleus_excitation(table::MassTable, system::FissioningSystem)
    system.spontaneous && return measurement(0.0, 0.0)
    Sₙ = neutron_separation_energy(table, system.compound)
    Sₙ === nothing && return nothing
    return Sₙ + system.incident_energy
end

"""
    total_excitation_energy(Q, E_compound, TKE) -> Union{Measurement{Float64},Nothing}

Total excitation energy shared by the complementary fragments at full acceleration, in MeV:

```
TXE = Q + E*_CN − TKE
```

Returns `nothing` when `TXE ≤ 0`, which is not a physical configuration and is the reason a
deterministic sweep over TKE produces fewer records than its declared grid.

`TKE` enters without an uncertainty because it is the swept independent variable of the model,
not a measurement.

**Uncertainty propagation stops here.** `TXE` carries the uncertainty of the mass excesses it
was built from, but the TXE partition methods that consume it — modelling at scission,
`E*_H/TXE (A_H)`, `R_T (A_H)` — are parameterizations with no quoted uncertainty, so `E*_L`
and `E*_H` inherit none. Carrying a `σ` past the partition would report the mass-excess
contribution while silently omitting the dominant model term.
"""
function total_excitation_energy(
    Q::Measurement{Float64},
    E_compound::Measurement{Float64},
    TKE::Real,
)
    TXE = Q + E_compound - TKE
    return value(TXE) > 0 ? TXE : nothing
end

"""
    kinetic_energy(system, fragment, TKE) -> Float64

Kinetic energy of one fragment in MeV, from conservation of momentum in a binary split:

```
KE(A) = TKE (A₀ − A) / A₀
```
"""
function kinetic_energy(system::FissioningSystem, fragment::Nuclide, TKE::Real)
    return TKE * (system.compound.A - fragment.A) / system.compound.A
end

"""
    Energetics

The energy balance of one fragmentation at one `TKE`, carrying the quantities that precede the
TXE partition.
"""
struct Energetics
    fragmentation::Fragmentation
    TKE::Float64
    Q::Measurement{Float64}
    TXE::Measurement{Float64}

    function Energetics(
        fragmentation::Fragmentation,
        TKE::Real,
        Q::Measurement{Float64},
        TXE::Measurement{Float64},
    )
        return new(fragmentation, Float64(TKE), Q, TXE)
    end
end

function Base.show(io::IO, energetics::Energetics)
    print(
        io,
        "Energetics(A_H=",
        energetics.fragmentation.heavy.A,
        ", Z_H=",
        energetics.fragmentation.heavy.Z,
        ", TKE=",
        energetics.TKE,
        ", TXE=",
        energetics.TXE,
        ")",
    )
    return nothing
end

"""
    energetics(table, system, fragmentation, TKE) -> Union{Energetics,Nothing}

The energy balance of one fragmentation at one kinetic energy, or `nothing` when the
configuration is unphysical — a missing mass excess, a non-positive `Q`, or `TXE ≤ 0`.
"""
function energetics(
    table::MassTable,
    system::FissioningSystem,
    fragmentation::Fragmentation,
    TKE::Real,
)
    Q = q_value(table, system, fragmentation.heavy)
    Q === nothing && return nothing
    E_compound = compound_nucleus_excitation(table, system)
    E_compound === nothing && return nothing
    TXE = total_excitation_energy(Q, E_compound, TKE)
    TXE === nothing && return nothing
    return Energetics(fragmentation, Float64(TKE), Q, TXE)
end

"""
    ExcludedConfiguration

Why one `(fragmentation, TKE)` configuration was left out of a sweep.

The prototype drops such configurations silently, so coverage cannot be audited from its
output: the archived ²⁵²Cf(SF) run writes 9819 records against a declared grid of 12 495.
"""
struct ExcludedConfiguration
    fragmentation::Fragmentation
    TKE::Float64
    reason::Symbol

    function ExcludedConfiguration(fragmentation::Fragmentation, TKE::Real, reason::Symbol)
        return new(fragmentation, Float64(TKE), reason)
    end
end

"""
    _validate_kinetic_energies(kinetic_energies)

Check the swept `TKE` grid the way [`fragmentation_domain`](@ref) checks the mass range.

It is the other swept axis of the model and it was reaching the sweep unvalidated, which is not
harmless: `TXE = Q + E*_CN − TKE`, so a negative `TKE` *raises* the excitation energy and the
configuration passes the `TXE > 0` test. A grid of `[-50.0, -10.0]` produced 490 accepted
configurations and no complaint. An empty grid produced a sweep of zero records, which looks
like a physical exclusion rather than a bad argument.

Ascending and distinct is required rather than merely tidy: the sweep emits in grid order, so a
permuted grid permutes the output for no stated reason.
"""
function _validate_kinetic_energies(kinetic_energies::AbstractVector{<:Real})
    isempty(kinetic_energies) &&
        throw(ArgumentError("the kinetic energy grid is empty; nothing would be swept"))
    # Ascending and distinct fall out of one adjacency test, which also keeps this off Base's
    # `allunique` specialisation for abstractly typed vectors — a path the static analysis
    # cannot follow through its reinterpret and reshape layers.
    start = firstindex(kinetic_energies)
    for index in eachindex(kinetic_energies)
        TKE = Float64(kinetic_energies[index])
        isfinite(TKE) ||
            throw(ArgumentError("kinetic energy number $index of the grid is not finite"))
        TKE > 0 || throw(
            ArgumentError(
                "total kinetic energy must be positive; entry $index is $TKE MeV",
            ),
        )
        index == start && continue
        previous = Float64(kinetic_energies[index - 1])
        previous < TKE || throw(
            ArgumentError(
                "the grid must ascend strictly; entry $index is $TKE after $previous",
            ),
        )
    end
    return nothing
end

"""
    sweep_energetics(table, system, domain, kinetic_energies)
        -> (Vector{Energetics}, Vector{ExcludedConfiguration})

Evaluate the energy balance over the whole domain and TKE grid, returning both what survived
and what did not, each exclusion carrying its reason: `:missing_mass`, `:nonpositive_q` or
`:nonpositive_txe`.
"""
function sweep_energetics(
    table::MassTable,
    system::FissioningSystem,
    domain::FragmentationDomain,
    kinetic_energies::AbstractVector{<:Real},
)
    _validate_kinetic_energies(kinetic_energies)

    accepted = Energetics[]
    excluded = ExcludedConfiguration[]
    sizehint!(accepted, length(domain) * length(kinetic_energies))

    E_compound = compound_nucleus_excitation(table, system)
    E_compound === nothing && throw(
        ArgumentError(
            "the neutron separation energy of the compound nucleus $(system.compound) is not " *
            "available from $(table.source); the excitation energy cannot be formed",
        ),
    )

    for fragmentation in domain
        Q = q_value(table, system, fragmentation.heavy)
        if Q === nothing
            reason = _q_failure_reason(table, system, fragmentation.heavy)
            for TKE in kinetic_energies
                push!(excluded, ExcludedConfiguration(fragmentation, Float64(TKE), reason))
            end
            continue
        end
        for TKE in kinetic_energies
            TXE = total_excitation_energy(Q, E_compound, TKE)
            if TXE === nothing
                push!(
                    excluded,
                    ExcludedConfiguration(fragmentation, Float64(TKE), :nonpositive_txe),
                )
            else
                push!(accepted, Energetics(fragmentation, Float64(TKE), Q, TXE))
            end
        end
    end
    return accepted, excluded
end

function _q_failure_reason(table::MassTable, system::FissioningSystem, heavy::Nuclide)
    light = complementary_fragment(system, heavy)
    if !haskey(table, heavy) || !haskey(table, light) || !haskey(table, system.compound)
        return :missing_mass
    end
    return :nonpositive_q
end

"""
    exclusion_summary(excluded) -> Dict{Symbol,Int}

Count the excluded configurations by reason, for the run's provenance record.
"""
function exclusion_summary(excluded::AbstractVector{ExcludedConfiguration})
    counts = Dict{Symbol, Int}()
    for exclusion in excluded
        counts[exclusion.reason] = get(counts, exclusion.reason, 0) + 1
    end
    return counts
end
