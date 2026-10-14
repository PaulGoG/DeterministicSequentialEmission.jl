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

The full sequential emission from one initial fragment at one `(A, Z, TKE)`, characterised by
two distinct counts.

# Fields
- `possible`: `n`, the number of prompt neutrons emitted. The chain continues while the
  residual can afford a separation energy once the kinetic energy carried away by the neutrons
  already emitted is returned. `n` is an integer, so the post-neutron fragment mass
  `A_p = A − n` is an integer without rounding.
- `solved`: `k_max ≤ n`, the number of sequences whose equation admits an acceptable
  temperature. `T`, `Ē_r` and `⟨ε⟩` are averaged over these `k_max` sequences only.
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

`ν` is the sequence average of the emission index. Every per-sequence quantity `q_k` is reduced
by the same operator, eq. (12) of Tudora, *EPJA* 54, 87 (2018), doi:10.1140/epja/i2018-12521-7;
eq. (4) of Tudora, *EPJA* 56, 84 (2020), doi:10.1140/epja/s10050-020-00059-2; eq. (20) of Tudora, EPJ Web Conf. 256, 00017 (2021),
doi:10.1051/epjconf/202125600017,

```
⟨q⟩ = (1/n) Σ_{k=1}^{n} q_k
```

evaluated at `q_k = k`, so `ν = (1 + 2 + … + n)/n = (n+1)/2`. It is the multiplicity of the
averaged emission sequence and the quantity compared with measurement. It differs from `n`, the
number of sequences, which [`post_neutron_mass`](@ref) uses.
"""
neutron_multiplicity(chain::EmissionChain) =
    chain.possible == 0 ? 0.0 : (chain.possible + 1) / 2

"""
    post_neutron_mass(chain) -> Int

The mass number of the last residual fragment, `A_p = A − n`, an integer by construction.

`A_p` uses `n`, not `ν`: it is a property of the chain, the nucleus on which the cascade ends,
rather than a sequence average. A global treatment yields only an averaged `ν(A,Z,TKE)`, whose
`A − ν` must be rounded (Tudora, *EPJA* 58, 126 (2022), doi:10.1140/epja/s10050-022-00766-y, §2.4). Consequently
`⟨A⟩ − ⟨A_p⟩` equals the yield average of `n`, not `⟨ν⟩`. It is the product mass under
[`SequenceCount`](@ref); the default [`EmissionProbabilityCount`](@ref) of [`product_masses`](@ref)
spreads the product over the chain's `n` residuals, and the mass loss is then `ν`.
"""
post_neutron_mass(chain::EmissionChain) = chain.fragment.A - chain.possible

"""
    gamma_energy(chain, excitation) -> Float64

The excitation energy of a fragment that its neutrons leave, in MeV: what the model has left for
prompt γ rays, conversion electrons and isomers. Under the count of
[`neutron_multiplicity`](@ref), the chain ending equiprobably after `1 … n` emissions, the `k`-th
neutron is emitted with probability `P(N ≥ k) = (n+1−k)/n` and carries `Sₙ^(k) + ⟨ε⟩_k` away:

```
E_γ = E* − Σ_k P(N ≥ k) (Sₙ^(k) + ⟨ε⟩_k)
```

A sequence of the second stage carries no temperature and is charged its separation energy
alone, `⟨ε⟩_k = 0`, as the spectrum charges it, so the balance `E* = Σ_k P(N ≥ k)(Sₙ^(k) + ⟨ε⟩_k)
+ E_γ` closes on the energies the run's neutrons are given. The endings beyond `k_max`, taken
alone, would close below zero, since the second stage counts on `Sₙ` once the kinetic energy has
been returned; the mean over the endings is what the count predicts. The weights are those of
[`EmissionProbabilityWeighting`](@ref) whatever the spectrum's weighting, since they belong to the
count. A chain that emits nothing leaves its whole excitation.
"""
function gamma_energy(chain::EmissionChain, excitation::Real)
    chain.possible == 0 && return Float64(excitation)
    carried = 0.0
    for sequence in chain.sequences
        ε = something(sequence.mean_neutron_energy, 0.0)
        carried +=
            sequence_weight(EmissionProbabilityWeighting(), chain, sequence) *
            (sequence.separation_energy + ε)
    end
    return Float64(excitation) - carried
end

"""
    SequenceWeighting

How the emission sequences of one fragment are weighted in a quantity carried by its neutrons:
the centre-of-mass and laboratory spectra and the mean neutron energy.

[`neutron_multiplicity`](@ref) is `ν = (n+1)/2`, the mean of a number of emitted neutrons that
is equiprobable on `1 … n`. The `k`-th neutron of the chain is then emitted with probability
`P(N ≥ k) = (n+1−k)/n`, and a quantity averaged over the emitted neutrons carries that weight:

```
⟨q⟩ = Σ_k (n+1−k) q_k / Σ_k (n+1−k)
```

[`EmissionProbabilityWeighting`](@ref), the default, applies it. The weights fall linearly along
the chain, as the probability of reaching a step does in every treatment that carries a
distribution of excitation energy whose width is comparable with one neutron step: Terrell,
*Phys. Rev.* **113**, 527 (1959), doi:10.1103/PhysRev.113.527; Märten and Seeliger, *J. Phys. G*
**10**, 349 (1984), doi:10.1088/0305-4616/10/3/010; Kawano *et al.*, *Phys. Rev. C* **104**,
014611 (2021), doi:10.1103/PhysRevC.104.014611. [`UniformSequenceWeighting`](@ref) counts every
sequence once, `⟨q⟩ = (1/k_max) Σ_k q_k`, the operator of eq. (12) of Tudora, *EPJA* 54, 87
(2018), doi:10.1140/epja/i2018-12521-7, and of the thesis. It takes the last neutron of a chain
as emitted as often as the first, which over-represents the cold end of the chain and gives a
spectrum softer than the multiplicity it goes with implies.

The weighting applies to what the neutrons carry. The residual temperature and the residual
energy written against fragment mass belong to the nucleus the chain ends on; with the count
equiprobable on `1 … n` each residual is the last with the same probability, so their average
over the sequences is uniform under either weighting.
"""
abstract type SequenceWeighting end

"""
    UniformSequenceWeighting()

Every sequence that carries the quantity counts once.
"""
struct UniformSequenceWeighting <: SequenceWeighting end

"""
    EmissionProbabilityWeighting()

The `k`-th sequence of a chain of `n` counts with `(n+1−k)/n`.
"""
struct EmissionProbabilityWeighting <: SequenceWeighting end

"""
    sequence_weight(weighting, chain, sequence) -> Float64

The weight of one sequence of `chain` under `weighting`: `1` under
[`UniformSequenceWeighting`](@ref), `(n+1−k)/n` under [`EmissionProbabilityWeighting`](@ref),
with `n = chain.possible` and `k = sequence.index`.
"""
sequence_weight(::UniformSequenceWeighting, ::EmissionChain, ::EmissionSequence) = 1.0
sequence_weight(
    ::EmissionProbabilityWeighting,
    chain::EmissionChain,
    sequence::EmissionSequence,
) = (chain.possible + 1 - sequence.index) / chain.possible

"""
    sequence_average(chain, quantity[, weighting]) -> Union{Float64,Nothing}

Average a per-sequence quantity over the sequences where it is defined, `k_max` of them, with
the weights of `weighting`, a [`SequenceWeighting`](@ref); uniform where none is given.

`quantity` is one of `:temperature`, `:residual_energy` or `:mean_neutron_energy`.
"""
function sequence_average(
    chain::EmissionChain,
    quantity::Symbol,
    weighting::SequenceWeighting,
)
    quantity in (:temperature, :residual_energy, :mean_neutron_energy) ||
        throw(ArgumentError("unknown per-sequence quantity $quantity"))
    total = 0.0
    norm = 0.0
    for sequence in chain.sequences
        value = getfield(sequence, quantity)
        value === nothing && continue
        weight = sequence_weight(weighting, chain, sequence)
        total += weight * value
        norm += weight
    end
    return norm > 0 ? total / norm : nothing
end

sequence_average(chain::EmissionChain, quantity::Symbol) =
    sequence_average(chain, quantity, UniformSequenceWeighting())

"""
    emission_chain(table, model, cross_section, fragment, excitation; kwargs...)
        -> EmissionChain

Solve the sequential equations for one fragment, starting from its excitation energy at full
acceleration. Two stages separate `k_max` from `n`.

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
    table::MassExcessTable,
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
    # returned. These carry no temperature but count towards `n`.
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
