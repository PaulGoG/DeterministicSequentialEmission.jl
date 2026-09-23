"""
    Fragmentation

One binary split of the compound nucleus into a complementary pair, with the probability of
that charge split.

The probability is a property of the pair, not of either fragment: with `Z_L = Z₀ − Z_H` and
the light fragment's centre at `Z₀ − Zₚ`, the Gaussian argument satisfies
`Z_L − (Z₀ − Zₚ) = −(Z_H − Zₚ)`, so both fragments carry the same `p(Z,A)`.
"""
struct Fragmentation
    heavy::Nuclide
    light::Nuclide
    probability::Float64

    function Fragmentation(heavy::Nuclide, light::Nuclide, probability::Real)
        probability >= 0 || throw(
            ArgumentError("a charge probability cannot be negative, got $probability"),
        )
        return new(heavy, light, Float64(probability))
    end
end

"""
    is_symmetric(fragmentation) -> Bool

Whether the split is mass-symmetric, `A_H == A_L`.
"""
is_symmetric(fragmentation::Fragmentation) = fragmentation.heavy.A == fragmentation.light.A

"""
    FragmentationDomain

The deterministic set of fragmentations the model sweeps: every heavy mass in
`heavy_masses`, each with `charges_per_mass` charge numbers around `Zₚ(A)`, paired with its
complement. Entries are sorted as whole rows by heavy mass and then by heavy charge.
"""
struct FragmentationDomain
    system::FissioningSystem
    entries::Vector{Fragmentation}
    heavy_masses::UnitRange{Int}
    charges_per_mass::Int

    function FragmentationDomain(
        system::FissioningSystem,
        entries::Vector{Fragmentation},
        heavy_masses::UnitRange{Int},
        charges_per_mass::Integer,
    )
        return new(system, entries, heavy_masses, Int(charges_per_mass))
    end
end

Base.length(domain::FragmentationDomain) = length(domain.entries)
Base.iterate(domain::FragmentationDomain, state...) = iterate(domain.entries, state...)
Base.eltype(::Type{FragmentationDomain}) = Fragmentation
Base.getindex(domain::FragmentationDomain, index) = domain.entries[index]

function Base.show(io::IO, domain::FragmentationDomain)
    print(
        io,
        "FragmentationDomain(",
        length(domain.entries),
        " fragmentations, A_H ∈ ",
        domain.heavy_masses,
        ", ",
        domain.charges_per_mass,
        " Z per A)",
    )
    return nothing
end

"""
    fragmentation_domain(system, distribution, heavy_masses; charges_per_mass = 5)
        -> FragmentationDomain

Build the fragmentation domain: heavy masses from symmetric fission to a maximum asymmetry in
steps of one mass unit, each carrying an odd number of charge numbers centred on
`Zₚ(A) = Z_UCD(A) + ΔZ(A)`, weighted by the isobaric charge distribution. `heavy_masses` is
validated against the compound nucleus, the nucleus that splits.

# Examples

```jldoctest
julia> system = spontaneous_fission(Nuclide(98, 252));

julia> domain = fragmentation_domain(system, mean_charge_distribution(), 126:174);

julia> length(domain)
245

julia> first(domain).heavy
Nuclide(Z=46, A=126)
```
"""
function fragmentation_domain(
    system::FissioningSystem,
    distribution::ChargeDistribution,
    heavy_masses::AbstractUnitRange{<:Integer};
    charges_per_mass::Integer = 5,
)
    isempty(heavy_masses) && throw(ArgumentError("the heavy fragment mass range is empty"))
    lower = first(heavy_masses)
    upper = last(heavy_masses)
    symmetric = symmetric_mass(system)
    lower >= symmetric || throw(
        ArgumentError(
            "the heavy fragment range must start at or above symmetric fission, " *
            "A₀/2 = $symmetric for $(system.compound); got A_H_min = $lower",
        ),
    )
    upper < system.compound.A || throw(
        ArgumentError(
            "the heavy fragment range must stay below the compound mass $(system.compound.A); " *
            "got A_H_max = $upper",
        ),
    )
    lower < upper || throw(
        ArgumentError(
            "the heavy fragment range must span more than one mass; got A_H_min = $lower and " *
            "A_H_max = $upper",
        ),
    )

    entries = Fragmentation[]
    sizehint!(entries, length(heavy_masses) * charges_per_mass)
    for A_heavy in heavy_masses
        ΔZ = charge_polarization(distribution, A_heavy)
        σ_Z = charge_dispersion(distribution, A_heavy)
        Zₚ = most_probable_charge(system, A_heavy, ΔZ)
        for Z_heavy in charge_numbers(Zₚ, charges_per_mass)
            # Both fragments must carry protons; a wide charge window or a Zₚ near the edge of
            # the range can leave one fragment with none.
            0 < Z_heavy < system.compound.Z || throw(
                ArgumentError(
                    "charge number $Z_heavy at A_H = $A_heavy leaves no protons for one " *
                    "fragment of $(system.compound); $charges_per_mass charge numbers about " *
                    "Zₚ = $(round(Zₚ; digits = 3)) is too many",
                ),
            )
            heavy = Nuclide(Z_heavy, A_heavy)
            light = complementary_fragment(system, heavy)
            push!(
                entries,
                Fragmentation(heavy, light, charge_probability(Z_heavy, Zₚ, σ_Z)),
            )
        end
    end

    sort!(entries; by = entry -> (entry.heavy.A, entry.heavy.Z))
    return FragmentationDomain(
        system,
        entries,
        UnitRange{Int}(lower, upper),
        Int(charges_per_mass),
    )
end

"""
    fragments(domain) -> Vector{Nuclide}

Every distinct nuclide the domain produces as a fragment, light and heavy, each once, in
ascending `(A, Z)`. Per-nuclide quantities such as `Sₙ(A,Z)` and `a(A,Z)` are evaluated once
over this set rather than at every configuration and emission sequence. The set carries no
probability, since a nuclide's weight depends on the quantity: `p(Z,A)` for building yields,
[`fragment_charge_probabilities`](@ref), or `Y(A,Z,TKE)` for averaging a model quantity. See
[`distinct_splits`](@ref) for the pair-level view, which is a different collection.
"""
function fragments(domain::FragmentationDomain)
    seen = Set{Nuclide}()
    for entry in domain.entries
        push!(seen, entry.heavy)
        push!(seen, entry.light)
    end
    return sort!(collect(seen); by = nuclide -> (nuclide.A, nuclide.Z))
end

"""
    distinct_splits(domain) -> Vector{Fragmentation}

Every physical split the domain contains, each once. Away from mass symmetry this is
`domain.entries`. At `A_H == A_L == A₀/2` the domain enumerates `(H=Z, L=Z₀−Z)` and
`(H=Z₀−Z, L=Z)`, which label the same split from both sides; `entries` keeps both, each carrying
its heavy fragment's own `p(Z,A)` as a fragment-indexed sweep requires, while this list keeps the
`Z_H >= Z_L` labelling, so the result does not depend on iteration order. A sum over pairs
iterates `entries` instead and halves the yield of every mass-symmetric one (`_pair_weight`),
which also covers the self-complementary split `(Z₀/2, Z₀/2)` that this list keeps at the yield
of two fragments.

# Examples

```jldoctest
julia> system = spontaneous_fission(Nuclide(98, 252));

julia> domain = fragmentation_domain(system, mean_charge_distribution(), 126:174);

julia> length(domain), length(distinct_splits(domain))
(245, 244)
```
"""
function distinct_splits(domain::FragmentationDomain)
    labelled = Set((entry.heavy.A, entry.heavy.Z) for entry in domain.entries)
    kept = Fragmentation[]
    for entry in domain.entries
        entry.heavy.A == entry.light.A || (push!(kept, entry); continue)
        # Mass-symmetric. Drop this labelling only if the mirror is itself in the domain —
        # the retained charge window is centred on Zₚ, which at A₀/2 need not be symmetric
        # about Z₀/2, so a split can have no mirror and must be kept whichever side labels it.
        mirror = (entry.heavy.A, entry.light.Z)
        if entry.heavy.Z >= entry.light.Z || !(mirror in labelled)
            push!(kept, entry)
        end
    end
    return kept
end

"""
    mass_range(domain) -> UnitRange{Int}

The full span of fragment masses the domain covers, from the lightest complement to the
heaviest fragment.
"""
function mass_range(domain::FragmentationDomain)
    lightest = domain.system.compound.A - last(domain.heavy_masses)
    return lightest:last(domain.heavy_masses)
end
