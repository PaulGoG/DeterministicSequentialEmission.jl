"""
    EmissionSequence

One emission sequence `k` of one fragment: the neutron separation energy of the emitting
nucleus, its level density parameter, the residual temperature, the residual excitation
energy and the mean centre-of-mass energy of the neutron emitted — all in MeV, except `a` in
MeV⁻¹.

`temperature` is `nothing` for a sequence that is energetically possible but for which the
equation admits no acceptable root. Those sequences still count towards `n`.
"""
struct EmissionSequence
    index::Int
    separation_energy::Float64
    level_density_parameter::Float64
    temperature::Union{Float64, Nothing}
    residual_energy::Union{Float64, Nothing}
    mean_neutron_energy::Union{Float64, Nothing}

    function EmissionSequence(
        index::Integer,
        separation_energy::Real,
        level_density_parameter::Real,
        temperature::Union{Real, Nothing},
        residual_energy::Union{Real, Nothing},
        mean_neutron_energy::Union{Real, Nothing},
    )
        return new(
            Int(index),
            Float64(separation_energy),
            Float64(level_density_parameter),
            temperature === nothing ? nothing : Float64(temperature),
            residual_energy === nothing ? nothing : Float64(residual_energy),
            mean_neutron_energy === nothing ? nothing : Float64(mean_neutron_energy),
        )
    end
end

is_solved(sequence::EmissionSequence) = sequence.temperature !== nothing

"""
    EmissionChain

The full sequential emission from one initial fragment at one `(A, Z, TKE)`.

Two counts characterise it and they are not interchangeable.

`possible` is `n`, the number of prompt neutrons emitted: the chain is continued while the
residual can still afford a separation energy once the kinetic energy carried away by the
neutrons already emitted is returned. Because `n` is an integer at every configuration, the
post-neutron fragment mass `A_p = A − n` is an integer without rounding — which is the reason
the model treats emission sequentially rather than through an integrated temperature
distribution.

`solved` is `k_max ≤ n`, the number of sequences whose equation admits an acceptable
temperature. Quantities defined only where a temperature exists — `T`, `Ē_r`, `⟨ε⟩` — are
averaged over `k_max`; padding the remainder with zeros would dilute a temperature with values
that are not temperatures.
"""
struct EmissionChain
    fragment::Nuclide
    sequences::Vector{EmissionSequence}
    possible::Int
    solved::Int

    function EmissionChain(fragment::Nuclide, sequences::Vector{EmissionSequence})
        possible = isempty(sequences) ? 0 : last(sequences).index
        solved = count(is_solved, sequences)
        return new(fragment, sequences, possible, solved)
    end
end

Base.length(chain::EmissionChain) = length(chain.sequences)
Base.iterate(chain::EmissionChain, state...) = iterate(chain.sequences, state...)
Base.eltype(::Type{EmissionChain}) = EmissionSequence

function Base.show(io::IO, chain::EmissionChain)
    print(
        io,
        "EmissionChain(",
        chain.fragment,
        ", n=",
        chain.possible,
        ", k_max=",
        chain.solved,
        ")",
    )
    return nothing
end

"""
    neutron_multiplicity(chain) -> Float64

The prompt neutron multiplicity `ν = (n+1)/2`, zero for a chain that emits nothing.

`ν` is a measurable quantity, averaged over the emission sequences of the configuration. It is
not the same as `n`, the number of sequences, and the two are not interchangeable even though
both describe emission from the same fragment.
"""
neutron_multiplicity(chain::EmissionChain) =
    chain.possible == 0 ? 0.0 : (chain.possible + 1) / 2

"""
    post_neutron_mass(chain) -> Int

The mass number of the last residual fragment, `A_p = A − n`. An integer by construction.
"""
post_neutron_mass(chain::EmissionChain) = chain.fragment.A - chain.possible

"""
    sequence_average(chain, quantity) -> Union{Float64,Nothing}

Average a per-sequence quantity over the sequences where it is defined, `k_max` of them.

`quantity` is one of `:temperature`, `:residual_energy` or `:mean_neutron_energy`.
"""
function sequence_average(chain::EmissionChain, quantity::Symbol)
    quantity in (:temperature, :residual_energy, :mean_neutron_energy) ||
        throw(ArgumentError("unknown per-sequence quantity $quantity"))
    total = 0.0
    counted = 0
    for sequence in chain.sequences
        value = getfield(sequence, quantity)
        value === nothing && continue
        total += value
        counted += 1
    end
    return counted == 0 ? nothing : total / counted
end

"""
    emission_chain(table, model, cross_section, fragment, excitation; kwargs...)
        -> EmissionChain

Solve the sequential equations for one fragment, starting from its excitation energy at full
acceleration.

Two stages, and the distinction between them is what separates `k_max` from `n`.

The first solves

```
Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²
```

while the residual can afford the separation energy, consuming `Sₙ + ⟨ε⟩` per step. It stops
at the first sequence whose equation admits no acceptable temperature.

The second returns the kinetic energy carried away by the neutrons emitted so far and
continues on `Sₙ` alone, counting the emissions that remain energetically possible. Those
sequences carry no temperature but do count towards `n`.
"""
function emission_chain(
    table::MassTable,
    model::LevelDensityModel,
    cross_section::EvaporationCrossSection,
    fragment::Nuclide,
    excitation::Real;
    minimum_temperature::Real = MINIMUM_TEMPERATURE,
    maximum_temperature::Union{Real, Nothing} = nothing,
)
    sequences = EmissionSequence[]
    residual = Float64(excitation)
    emitting = fragment
    carried_away = 0.0
    index = 1

    # First stage: sequences whose equation admits a temperature.
    while true
        Sₙ = neutron_separation_energy(table, emitting)
        Sₙ === nothing && break
        separation = value(Sₙ)
        residual >= separation || break

        daughter = Nuclide(emitting.Z, emitting.A - 1)
        a = level_density_parameter(model, daughter)
        a === nothing && break

        T = _solve_temperature(
            cross_section,
            a,
            residual - separation,
            daughter;
            minimum_temperature,
            maximum_temperature,
        )
        T === nothing && break

        ε = _mean_energy(cross_section, T, daughter)
        push!(
            sequences,
            EmissionSequence(index, separation, a, T, fermi_gas_energy(a, T), ε),
        )
        carried_away += ε
        residual = fermi_gas_energy(a, T)
        emitting = daughter
        index += 1
    end

    # Second stage: emissions that remain possible once the neutrons' kinetic energy is
    # returned. These have no temperature but are prompt neutrons all the same.
    residual += carried_away
    while true
        Sₙ = neutron_separation_energy(table, emitting)
        Sₙ === nothing && break
        separation = value(Sₙ)
        residual >= separation || break

        daughter = Nuclide(emitting.Z, emitting.A - 1)
        a = level_density_parameter(model, daughter)
        push!(
            sequences,
            EmissionSequence(
                index,
                separation,
                a === nothing ? NaN : a,
                nothing,
                nothing,
                nothing,
            ),
        )
        residual -= separation
        emitting = daughter
        index += 1
    end

    return EmissionChain(fragment, sequences)
end

function _solve_temperature(
    cross_section::ConstantCrossSection,
    a,
    available,
    ::Nuclide;
    minimum_temperature,
    maximum_temperature,
)
    return residual_temperature(
        cross_section,
        a,
        available;
        minimum_temperature,
        maximum_temperature,
    )
end

function _solve_temperature(
    cross_section::VariableCrossSection,
    a,
    available,
    daughter::Nuclide;
    minimum_temperature,
    maximum_temperature,
)
    return residual_temperature(
        cross_section,
        a,
        available;
        α = alpha_parameter(daughter.A),
        minimum_temperature,
        maximum_temperature,
    )
end

_mean_energy(cross_section::ConstantCrossSection, T, ::Nuclide) =
    mean_neutron_energy(cross_section, T)

_mean_energy(cross_section::VariableCrossSection, T, daughter::Nuclide) =
    mean_neutron_energy(cross_section, alpha_parameter(daughter.A), T)
