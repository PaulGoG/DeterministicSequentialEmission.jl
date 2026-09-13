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
complement.

Entries are sorted by heavy mass and then by heavy charge, as whole rows. The prototype sorted
the mass and charge columns independently and then looked each value up by `(A, Z)`, taking the
first match. Sorting whole rows is `O(N log N)` against that lookup's `O(N²)` and cannot depend
on a first match, but the stronger claim once made here — that the column-wise sort mispairs a
mass with a charge — was checked and could not be reproduced: over a full domain the two agree
on every row, because the charges are grouped by `unique(sort(A))` and so stay with their mass.
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
    fragmentation_domain(system, polarization, heavy_masses; charges_per_mass = 5)
        -> FragmentationDomain

Build the fragmentation domain: heavy masses from symmetric fission to a maximum asymmetry in
steps of one mass unit, each carrying an odd number of charge numbers centred on
`Zₚ(A) = Z_UCD(A) + ΔZ(A)`, weighted by the isobaric charge distribution.

`heavy_masses` is validated against the **compound** nucleus, which is the nucleus that
actually splits. Validating against the target instead — as the prototype does — happens to
agree only when the lower bound sits exactly at `A₀/2`.

# Examples

```jldoctest
julia> system = spontaneous_fission(Nuclide(98, 252));

julia> domain = fragmentation_domain(system, mean_charge_polarization(), 126:174);

julia> length(domain)
245

julia> first(domain).heavy
Nuclide(Z=46, A=126)
```
"""
function fragmentation_domain(
    system::FissioningSystem,
    polarization::ChargePolarization,
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
        ΔZ = polarization_at(polarization, A_heavy)
        rms = width_at(polarization, A_heavy)
        Zₚ = most_probable_charge(system, A_heavy, ΔZ)
        for Z_heavy in charge_numbers(Zₚ, charges_per_mass)
            # Both fragments have to be nuclei. Too many charge numbers per mass, or a Zₚ near
            # the edge of the range, walks the window off one end and asks for a fragment with
            # no protons or with all of them.
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
                Fragmentation(heavy, light, charge_probability(Z_heavy, Zₚ, rms)),
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
ascending `(A, Z)`.

The set, not a weighting. It exists so per-nuclide quantities — `Sₙ(A,Z)`, `a(A,Z)` — can be
evaluated once each rather than recomputed at every configuration and every emission sequence,
which is what makes the prototype's averaging routines quadratic.

**It carries no probability, deliberately.** A nuclide's weight depends on what is being asked:
the isobaric charge distribution `p(Z,A)` for building yields, which is
[`fragment_charge_probabilities`](@ref); the yield `Y(A,Z,TKE)` for averaging a model quantity;
neither of which is a property of the domain. An earlier version returned a summed probability
here, which was well defined and answered a question nothing asks.

See [`distinct_splits`](@ref) for the pair-level view, which is not the same collection.
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

Every physical split the domain contains, each once.

Away from mass symmetry this is just `domain.entries`. At `A_H == A_L == A₀/2` it is not: the
domain enumerates `(H=Z, L=Z₀−Z)` and `(H=Z₀−Z, L=Z)`, which label the **same** split from the
two sides. They are kept apart in `entries` because each labels a different heavy fragment and
carries that fragment's own `p(Z,A)` — correct for a fragment-indexed sweep, and what the
archived run writes. Counted as splits, they are one.

So anything summing a *pair* quantity — `ν_pair`, `TKE_p`, the pair spectrum — must iterate this
rather than `entries`, or it counts the symmetric splits twice.

The prototype guards the same case with the three-way branch in its partitioning routines: push
both fragments when `A_L != A_H`; push both once, skipping a repeat, when the masses match but
the charges do not; and push a single cell when the fragments are the same nuclide. The rule
here is the same, with the survivor chosen as the `Z_H >= Z_L` labelling so it does not depend
on iteration order.

# Examples

```jldoctest
julia> system = spontaneous_fission(Nuclide(98, 252));

julia> domain = fragmentation_domain(system, mean_charge_polarization(), 126:174);

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
