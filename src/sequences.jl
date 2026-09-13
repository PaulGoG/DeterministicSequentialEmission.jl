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

The prompt neutron multiplicity `ν`, zero for a chain that emits nothing.

`ν` is the sequence average of the emission index. The model reduces every per-sequence
quantity `q_k` the same way — Tudora, *EPJA* 54, 87 (2018) eq. (12); *EPJA* 56, 84 (2020)
eq. (4); the 2021 review eq. (20) —

```
⟨q⟩ = (1/n) Σ_{k=1}^{n} q_k
```

and `ν` is that operator evaluated at `q_k = k`, so `ν = (1 + 2 + … + n)/n = (n+1)/2`. It is
the multiplicity of the averaged emission sequence, and it is the quantity with an experimental
counterpart to compare against.

It is **deliberately not** `n`, and the two are not interchangeable. `n` counts the sequences
of one configuration; `ν` is a sequence average, obtained exactly as `⟨T⟩`, `⟨ε⟩` and `⟨Ē_r⟩`
are. [`post_neutron_mass`](@ref) takes `n`, not `ν`, and the reason is stated there.

Yield-weighted over the archived ²⁵²Cf(SF) run the two come out at `⟨ν⟩_pair = 3.48` against
`⟨n⟩_pair = 5.01`, with the measured multiplicity at 3.757. Both conventions are settled and
preserved; see the note in `post_neutron_mass` for what follows from holding both.
"""
neutron_multiplicity(chain::EmissionChain) =
    chain.possible == 0 ? 0.0 : (chain.possible + 1) / 2

"""
    post_neutron_mass(chain) -> Int

The mass number of the last residual fragment, `A_p = A − n`. An integer by construction.

`n`, not `ν`. This is a property of the chain rather than an average over it: the model steps
the cascade to exhaustion and `A − n` is the nucleus it ends on. Integer `A_p` without rounding
is the reason for treating emission sequentially at all — a global treatment yields only an
averaged `ν(A,Z,TKE)`, so its `A − ν` must be rounded, which obstructs independent fission
product yields (Tudora, *EPJA* 58, 126 (2022), §2.4).

The two conventions are held together deliberately and are not reconciled. Stated plainly so
nobody has to rediscover it: over the archived ²⁵²Cf(SF) run the chain removes 5.01 nucleons per
fission while `ν` reports 3.48, so a naive `⟨A⟩ − ⟨A_p⟩ = ⟨ν⟩` identity does not hold here and
should not be asserted anywhere downstream. `A_p` answers "which nucleus is left", `ν` answers
"how many neutrons does the averaged sequence carry"; they are different questions and the
model does not claim they have the same answer.
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
