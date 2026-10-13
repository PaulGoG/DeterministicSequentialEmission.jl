"""
    ProductCount

How many neutrons a fragment's emission chain removes from its mass: the rule that turns a chain
of `n` possible sequences into the distribution of the post-neutron product mass `A_p`. The
charge is unchanged, since the products are taken before β decay.

[`EmissionProbabilityCount`](@ref), the default, ends the chain on each of its `n` residuals with
probability `1/n`, the reading that gives [`neutron_multiplicity`](@ref) `ν = (n+1)/2` and the
emission probabilities `(n+1−k)/n` of [`EmissionProbabilityWeighting`](@ref). The mass the
products lose is then the multiplicity, `⟨A − A_p⟩ = ν`, configuration by configuration.
[`SequenceCount`](@ref) subtracts all `n`, the prescription of Tudora, *EPJA* **58**, 126
(2022), doi:10.1140/epja/s10050-022-00766-y, under which the products lose `⟨n⟩ ≈ 2ν − 1`.
"""
abstract type ProductCount end

"""
    EmissionProbabilityCount()

A chain of `n` possible sequences ends on its `k`-th residual with probability `1/n`,

```
P(A_p = A − k) = 1/n,    k = 1 … n,
```

and a chain with `n = 0` leaves the fragment at `A_p = A`. The mean mass loss is
`(1/n) Σ_k k = (n+1)/2 = ν`.
"""
struct EmissionProbabilityCount <: ProductCount end

"""
    SequenceCount()

Every possible sequence is emitted, `A_p = A − n` with probability 1: eq. (8) and §2.4 of Tudora,
*EPJA* **58**, 126 (2022), doi:10.1140/epja/s10050-022-00766-y. It is [`post_neutron_mass`](@ref)
as a distribution.
"""
struct SequenceCount <: ProductCount end

"""
    product_masses(chain, count) -> Vector{Tuple{Int, Float64}}

The post-neutron product masses of the fragment `chain` starts from, with their probabilities,
under the [`ProductCount`](@ref) `count`. The probabilities sum to one. `n` is `chain.possible`,
the number of possible sequences, those without a temperature included.

- [`EmissionProbabilityCount`](@ref): `(A − k, 1/n)` for `k = 1 … n`; `(A, 1)` for `n = 0`.
- [`SequenceCount`](@ref): `(A − n, 1)`.

# Examples

```jldoctest
julia> chain = EmissionChain(Nuclide(50, 132), [
           EmissionSequence(1, 7.0, 13.0, 0.9, 10.5, 1.4),
           EmissionSequence(2, 5.5, 13.1, nothing, nothing, nothing),
       ]);

julia> product_masses(chain, EmissionProbabilityCount())
2-element Vector{Tuple{Int64, Float64}}:
 (131, 0.5)
 (130, 0.5)

julia> product_masses(chain, SequenceCount())
1-element Vector{Tuple{Int64, Float64}}:
 (130, 1.0)
```
"""
function product_masses(chain::EmissionChain, ::EmissionProbabilityCount)
    A = chain.fragment.A
    n = chain.possible
    n == 0 && return [(A, 1.0)]
    return [(A - k, 1 / n) for k in 1:n]
end

product_masses(chain::EmissionChain, ::SequenceCount) = [(post_neutron_mass(chain), 1.0)]

"""
    _product_distributions(chains, weights, compound_mass, count, total) -> NamedTuple

Contract the emission chains of a run into the post-neutron product distributions under the
[`ProductCount`](@ref) `count`. `chains` is the list [`run_pipeline`](@ref) builds, one
`(fragmentation, TKE, light, heavy, split)` per solved configuration; `weights` the yield
weights keyed by `(A, Z, TKE)`; `compound_mass` `A₀`.

Each fragment configuration `(A, Z, TKE)` with a chain and a weight `w` contributes once, through
the chain its multiplicity is taken from (a later occurrence replaces an earlier one, as in the
multiplicity store), with `w p` at each product mass `(A_p, p)` of [`product_masses`](@ref):

```
Y_p(A_p) = 𝒩 Σ w p,    Y_p(A_p, Z) = 𝒩 Σ_{Z fixed} w p,    Y_p(N_p = A_p − Z) = 𝒩 Σ w p
```

At `A₀/2` a configuration is reached from both labellings of a split, and twice from the
self-complementary one; it is counted once, as the input yield counts it. `𝒩 = total / Σ w`, the
sum running over the configurations solved, so that `Σ_{A_p} Y_p(A_p) = total`: the products of
the solved configurations stand for the whole distribution, as in the multiplicity contraction.

The kinetic energies are `KE = TKE (A₀ − A)/A₀` and `KE_p = KE A_p/A`, the velocity unchanged on
average; `KE_p(A_p)` and `KE_p(Z)` are means over `w p`. `TKE_p = KE_p,L + KE_p,H` is taken over
the pairs, the two products independent, with weight `w_pair p_L p_H`, `w_pair` the heavy
fragment's weight halved at mass symmetry (`_pair_weight`), and keyed by the heavy product
mass. Means with no positive weight are omitted.

Returns `mass_yield`, `yield` (`(A_p, Z, Y)`, sorted), `neutron_yield`, `kinetic_energy`,
`kinetic_energy_by_charge`, `total_kinetic_energy`, and the totals `mass_loss`, the mean
`⟨A − A_p⟩` per fragment over `w`, and `mean_total_kinetic_energy`, the mean `TKE_p` over the pairs.
"""
function _product_distributions(
    chains,
    weights::FragmentQuantity,
    compound_mass::Integer,
    count::ProductCount,
    total::Real,
)
    A₀ = Int(compound_mass)

    # One chain per fragment configuration, in the order first reached; a later occurrence
    # replaces the chain but keeps the place, so the sums run in a fixed order.
    position = Dict{Tuple{Int, Int, Float64}, Int}()
    fragments = Tuple{Tuple{Int, Int, Float64}, EmissionChain}[]
    for (entry, TKE, light, heavy, _) in chains
        for (fragment, chain) in ((entry.light, light), (entry.heavy, heavy))
            key = (fragment.A, fragment.Z, Float64(TKE))
            index = get(position, key, 0)
            if index == 0
                push!(fragments, (key, chain))
                position[key] = length(fragments)
            else
                fragments[index] = (key, chain)
            end
        end
    end

    by_mass = Dict{Int, Float64}()
    by_mass_charge = Dict{Tuple{Int, Int}, Float64}()
    by_neutron = Dict{Int, Float64}()
    energy_by_mass = Dict{Int, NTuple{2, Float64}}()
    energy_by_charge = Dict{Int, NTuple{2, Float64}}()
    weight_sum = 0.0
    loss = 0.0
    for (key, chain) in fragments
        raw = get(weights.values, key, nothing)
        raw === nothing && continue
        w = Float64(value(raw))
        A, Z, TKE = key
        weight_sum += w
        KE = TKE * (A₀ - A) / A₀
        for (A_p, p) in product_masses(chain, count)
            share = w * p
            by_mass[A_p] = get(by_mass, A_p, 0.0) + share
            by_mass_charge[(A_p, Z)] = get(by_mass_charge, (A_p, Z), 0.0) + share
            by_neutron[A_p - Z] = get(by_neutron, A_p - Z, 0.0) + share
            KE_p = KE * A_p / A
            _accumulate!(energy_by_mass, A_p, share, KE_p)
            _accumulate!(energy_by_charge, Z, share, KE_p)
            loss += share * (A - A_p)
        end
    end

    energy_by_heavy_mass = Dict{Int, NTuple{2, Float64}}()
    pair_sum = 0.0
    pair_energy = 0.0
    for (entry, TKE, light, heavy, _) in chains
        raw = get(weights.values, (entry.heavy.A, entry.heavy.Z, Float64(TKE)), nothing)
        raw === nothing && continue
        w = _pair_weight(entry, Float64(value(raw)))
        A_L = entry.light.A
        A_H = entry.heavy.A
        KE_L = TKE * (A₀ - A_L) / A₀
        KE_H = TKE * (A₀ - A_H) / A₀
        light_products = product_masses(light, count)
        for (A_pH, p_H) in product_masses(heavy, count), (A_pL, p_L) in light_products
            share = w * p_L * p_H
            TKE_p = KE_L * A_pL / A_L + KE_H * A_pH / A_H
            _accumulate!(energy_by_heavy_mass, A_pH, share, TKE_p)
            pair_sum += share
            pair_energy += share * TKE_p
        end
    end

    𝒩 = weight_sum > 0 ? Float64(total) / weight_sum : NaN
    yields = sort!([(A_p, Z, 𝒩 * Y) for ((A_p, Z), Y) in by_mass_charge])
    return (
        mass_yield = _scaled_curve(by_mass, 𝒩),
        yield = yields,
        neutron_yield = _scaled_curve(by_neutron, 𝒩),
        kinetic_energy = _mean_curve(energy_by_mass),
        kinetic_energy_by_charge = _mean_curve(energy_by_charge),
        total_kinetic_energy = _mean_curve(energy_by_heavy_mass),
        mass_loss = weight_sum > 0 ? loss / weight_sum : NaN,
        mean_total_kinetic_energy = pair_sum > 0 ? pair_energy / pair_sum : NaN,
    )
end

# Weighted sum and weight of a mean, accumulated in place.
function _accumulate!(store::Dict{Int, NTuple{2, Float64}}, key::Int, weight, quantity)
    weighted, norm = get(store, key, (0.0, 0.0))
    store[key] = (weighted + weight * quantity, norm + weight)
    return store
end

function _scaled_curve(store::Dict{Int, Float64}, scale::Float64)
    out = Vector{Tuple{Float64, Float64}}()
    for x in sort!(collect(keys(store)))
        push!(out, (Float64(x), scale * store[x]))
    end
    return out
end

function _mean_curve(store::Dict{Int, NTuple{2, Float64}})
    out = Vector{Tuple{Float64, Float64}}()
    for x in sort!(collect(keys(store)))
        weighted, norm = store[x]
        norm > 0 && push!(out, (Float64(x), weighted / norm))
    end
    return out
end
