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

Entries are sorted by heavy mass and then by heavy charge, as whole rows. The prototype
sorted the mass and charge columns independently and then looked values up by `(A, Z)`, which
can pair a mass with a charge that never occurred together.
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

    entries = Fragmentation[]
    sizehint!(entries, length(heavy_masses) * charges_per_mass)
    for A_heavy in heavy_masses
        ΔZ = polarization_at(polarization, A_heavy)
        rms = width_at(polarization, A_heavy)
        Zₚ = most_probable_charge(system, A_heavy, ΔZ)
        for Z_heavy in charge_numbers(Zₚ, charges_per_mass)
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
    fragments(domain) -> Vector{Tuple{Nuclide,Float64}}

Every distinct nuclide the domain produces as a fragment, with the total weight it is produced
with — the sum of the probabilities of every fragmentation it appears in, on either side.

Summing rather than taking one value matters only at mass symmetry, and there it is the whole
question. Away from `A₀/2` a nuclide reaches the domain through exactly one fragmentation, and
the pair's probability `p(Z_H, A_H)` is also its own Gaussian weight: the complement's most
probable charge is `Z₀ − Zₚ(A_H)`, so `p(Z_L, A_L)` evaluates to the same number by mirror
symmetry. At `A_H == A_L == A₀/2` that stops holding. The nuclide `(Z, A₀/2)` is then the heavy
member of one fragmentation and the light member of another — two different splits, reached with
two different probabilities — and neither one alone is "its" weight. For ²⁵²Cf at `A = 126`,
`Z = 50`, the two are `0.4058` and `0.0143`, a factor of 28 apart.

Taking the first, as a `get!` does, silently assigned each symmetric nuclide the probability of
its charge mirror and lost 40 % of the weight of that mass row. The sum, `0.4201` here, is the
expected number of times the nuclide is produced per fission over the domain, which is what a
fragment-indexed sum wants. A nuclide that is its own complement, `Z == Z₀/2`, is counted twice
for one fragmentation, correctly: that split really does make two of it.

This is a fragment census, not the isobaric charge distribution. Anything needing `p(Z,A)`
itself — the yield construction `Y(A,Z,TKE) = p(Z,A) Y(A,TKE)` above all — must take it from
[`charge_probability`](@ref) against that mass's own `Zₚ`, not from here.
"""
function fragments(domain::FragmentationDomain)
    seen = Dict{Nuclide, Float64}()
    for entry in domain.entries
        seen[entry.heavy] = get(seen, entry.heavy, 0.0) + entry.probability
        seen[entry.light] = get(seen, entry.light, 0.0) + entry.probability
    end
    ordered = collect(seen)
    sort!(ordered; by = pair -> (pair.first.A, pair.first.Z))
    return [(nuclide, probability) for (nuclide, probability) in ordered]
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
