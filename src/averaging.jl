"""
    FragmentQuantity

A model quantity tabulated over fragment configurations, `q(A, Z, TKE)`.

The rank-3 array the whole model is built on. Every observable is a contraction of one of these
against a weight — the isobaric charge distribution `p(Z,A)` when reducing over charge alone,
the yield `Y(A,Z,TKE)` otherwise — so the averaging routines take a quantity and a weight and
differ only in which indices they sum.

Values may be `Float64` or `Measurement{Float64}`. That choice is how a quantity says whether it
carries an uncertainty, and it is the whole of the decision: everything downstream of the TXE
partition has none to carry, and comes out bare without anyone writing a rule for it.
"""
struct FragmentQuantity{T <: Union{Float64, Measurement{Float64}}}
    values::Dict{Tuple{Int, Int, Float64}, T}
end

FragmentQuantity{T}() where {T} = FragmentQuantity(Dict{Tuple{Int, Int, Float64}, T}())

Base.length(q::FragmentQuantity) = length(q.values)
Base.getindex(q::FragmentQuantity, key) = q.values[key]
Base.get(q::FragmentQuantity, key, default) = get(q.values, key, default)
Base.haskey(q::FragmentQuantity, key) = haskey(q.values, key)
Base.setindex!(q::FragmentQuantity, value, key) = (q.values[key] = value)
Base.keys(q::FragmentQuantity) = keys(q.values)

function Base.show(io::IO, q::FragmentQuantity{T}) where {T}
    print(io, "FragmentQuantity{", T, "}(", length(q.values), " configurations)")
    return nothing
end

"""
    average(quantity, weights; over) -> Vector{Tuple{Float64,V}}

Contract `q(A,Z,TKE)` against a weight over the indices named in `over`, leaving the rest.

```
⟨q⟩(x) = Σ q(x,…) W(x,…) / Σ W(x,…)
```

`over` is any subset of `(:mass, :charge, :kinetic_energy)`; the result is indexed by whatever
is left, ascending. Reducing over all three gives a single number — use
[`total_average`](@ref) for that.

**A configuration enters only where both the quantity and the weight exist**, and it enters both
sums or neither. That is what keeps a coarse sweep from being diluted by a finer yield grid:
where the model solved nothing, nothing is averaged, rather than a zero being averaged in. The
prototype does the same with its `isassigned` guards, and it is the reason a grid mismatch
costs sampling rather than magnitude.

Uncertainties propagate by the standard expression for a weighted mean, with both terms,

```
δ²⟨q⟩ = Σᵢ [ (Wᵢ/ΣW) δqᵢ ]² + Σᵢ [ ((qᵢ − ⟨q⟩)/ΣW) δWᵢ ]²
```

the first from the quantity, the second from the weights. Neither is written out here: the
arithmetic is done in `Measurement`s and the derivatives fall out, which also keeps the terms
correlated where they share a source — `Y(A,Z,TKE)` cells at one `(A,TKE)` are `p(Z,A)` times a
single measurement, and treating them as independent understates the result by a third.
"""
function average(
    quantity::FragmentQuantity{Q},
    weights::FragmentQuantity{W};
    over::NTuple{N, Symbol},
) where {Q, W, N}
    _validate_indices(over)
    length(over) == 2 || throw(
        ArgumentError(
            "average reduces to a one-dimensional distribution, so it takes exactly two " *
            "indices to reduce over; got $over",
        ),
    )
    kept = only(_kept_indices(over))
    P = typeof(one(Q) * one(W))
    numerator = Dict{Float64, P}()
    denominator = Dict{Float64, W}()

    for (key, quantity_value) in quantity.values
        weight = get(weights, key, nothing)
        weight === nothing && continue
        label = Float64(key[kept])
        numerator[label] = get(numerator, label, zero(P)) + quantity_value * weight
        denominator[label] = get(denominator, label, zero(W)) + weight
    end

    out = Tuple{Float64, P}[]
    for (label, total) in denominator
        _positive(total) || continue
        push!(out, (label, numerator[label] / total))
    end
    sort!(out; by = first)
    return out
end

"""
    total_average(quantity, weights; masses = nothing) -> Union{Float64,Measurement{Float64}}

Contract over all three indices, optionally restricting to a set of fragment masses.

`masses` is how the light and heavy groups are separated: pass the light range for `⟨q⟩_L`, the
heavy range for `⟨q⟩_H`, and nothing for the whole distribution. The prototype takes the same
argument for the same reason.

Both ranges include `A₀/2` when the domain reaches symmetric fission, because both fragments of
a symmetric split carry that mass. The two groups therefore sum to slightly more than twice the
whole, which is a property of the grouping rather than an error.

`nothing` when no configuration carries both a value and a weight.
"""
function total_average(
    quantity::FragmentQuantity{Q},
    weights::FragmentQuantity{W};
    masses = nothing,
) where {Q, W}
    P = typeof(one(Q) * one(W))
    numerator = zero(P)
    denominator = zero(W)
    counted = 0
    for (key, quantity_value) in quantity.values
        masses === nothing || key[1] in masses || continue
        weight = get(weights, key, nothing)
        weight === nothing && continue
        numerator += quantity_value * weight
        denominator += weight
        counted += 1
    end
    counted == 0 && return nothing
    _positive(denominator) || return nothing
    return numerator / denominator
end

_positive(x::Real) = x > 0
_positive(x::Measurement) = value(x) > 0

const _INDEX_NAMES = (:mass, :charge, :kinetic_energy)

function _validate_indices(over::NTuple{N, Symbol}) where {N}
    for name in over
        name in _INDEX_NAMES ||
            throw(ArgumentError("unknown index $name; expected any of $(_INDEX_NAMES)"))
    end
    allunique(over) ||
        throw(ArgumentError("indices to reduce over must be distinct: $over"))
    length(over) < 3 || throw(
        ArgumentError("reducing over every index gives a number; use `total_average`"),
    )
    return nothing
end

function _kept_indices(over::NTuple{N, Symbol}) where {N}
    kept = Int[]
    for (position, name) in pairs(_INDEX_NAMES)
        name in over || push!(kept, position)
    end
    return kept
end

"""
    weights_from(yields) -> FragmentQuantity

Use a fragment yield distribution as the weight of a contraction.
"""
weights_from(yields::FragmentYield) = FragmentQuantity(copy(yields.values))

"""
    weights_from(domain, distribution, kinetic_energies) -> FragmentQuantity

Use the isobaric charge distribution as the weight, for a run with no experimental yield.

`p(Z,A)` alone, broadcast over every kinetic energy the quantity is tabulated at — which is what
the prototype falls back on when no `Y(A,TKE)` is supplied, and what its `main_DSE.out` writes in
place of a yield. Reducing over charge gives the same answer either way: `Y = p·Y(A,TKE)` and
`Y(A,TKE)` does not depend on `Z`, so it cancels between numerator and denominator.

Reducing over kinetic energy does **not**. A flat weight makes every `TKE` on the grid equally
likely, so the contraction sits at the midpoint of the grid rather than where the measurement
does, and every quantity built on `TXE` moves with it. Where no joint `Y(A,TKE)` was measured,
[`factorized_yield`](@ref) reconstructs one from the two marginals that were, which is a far
better weight than this.
"""
function weights_from(
    domain::FragmentationDomain,
    distribution::ChargeDistribution,
    kinetic_energies,
)
    probabilities = fragment_charge_probabilities(domain, distribution)
    weights = Dict{Tuple{Int, Int, Float64}, Float64}()
    for (nuclide, p) in probabilities, TKE in kinetic_energies
        weights[(nuclide.A, nuclide.Z, Float64(TKE))] = p
    end
    return FragmentQuantity(weights)
end
