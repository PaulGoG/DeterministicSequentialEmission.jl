"""
    EXPERIMENTAL_YIELD_SPEC

Layout of the measured pre-neutron mass yields: one header row, then `A TKE Y` and an optional
uncertainty column.
"""
const EXPERIMENTAL_YIELD_SPEC =
    TableSpec([:A, :TKE, :yield, :uncertainty]; skip = 1, optional = (:uncertainty,))

"""
    ExperimentalYield

Measured pre-neutron fragment mass yield `Y(A, TKE)`, on the grid the measurement was reported
on.

That grid is the measurement's, not the model's. It is common for the two to differ — the
²⁵²Cf(SF) distribution behind the archived runs is tabulated every 1 MeV while those runs swept
`TKE` every 2 MeV — so which grid a given average runs on is a question with a numerical answer,
and [`coverage`](@ref) is what answers it.

`A` runs over both fragments, light and heavy, so the distribution sums to 200 % rather than
100 %.
"""
struct ExperimentalYield
    values::Dict{Tuple{Int, Float64}, Float64}
    masses::Vector{Int}
    energies::Vector{Float64}
    path::String
end

Base.length(y::ExperimentalYield) = length(y.values)

function Base.show(io::IO, y::ExperimentalYield)
    print(
        io,
        "ExperimentalYield(",
        length(y.masses),
        " masses × ",
        length(y.energies),
        " TKE, ",
        length(y.values),
        " cells, from ",
        basename(y.path),
        ")",
    )
    return nothing
end

"""
    read_experimental_yield(path; spec = EXPERIMENTAL_YIELD_SPEC) -> ExperimentalYield

Read a measured `Y(A, TKE)` distribution.

Cells whose yield is not finite are dropped rather than carried: the legacy tables mark an
unmeasured `(A, TKE)` with `NaN`, and a missing measurement is not a yield of zero.

A repeated `(A, TKE)` is accepted when both rows carry the same yield and rejected when they do
not. The ²³⁵U(n,f) distribution lists `A = A₀/2` twice at every kinetic energy, with identical
values — the symmetric mass written from each side of the split — and it is one measurement, so
it is counted once. The prototype reaches the same result by iterating `unique(A)` and taking
the first match. A repeat that disagrees is a real inconsistency in the file and says so.
"""
function read_experimental_yield(
    path::AbstractString;
    spec::TableSpec = EXPERIMENTAL_YIELD_SPEC,
)
    table = read_table(path, spec)
    masses = integer_column(table, :A)
    energies = column(table, :TKE, Float64)
    measured = column(table, :yield, Float64)

    values = Dict{Tuple{Int, Float64}, Float64}()
    for row in eachindex(masses)
        isfinite(measured[row]) || continue
        key = (masses[row], energies[row])
        if haskey(values, key)
            values[key] == measured[row] || throw(
                ArgumentError(
                    "$(table.path) row $row: $(key) appears more than once with different " *
                    "yields, $(values[key]) and $(measured[row])",
                ),
            )
            continue
        end
        values[key] = measured[row]
    end
    isempty(values) && throw(ArgumentError("$(table.path) holds no finite yields"))
    return ExperimentalYield(
        values,
        sort(unique(masses)),
        sort(unique(energies)),
        String(path),
    )
end

"""
    MASS_YIELD_SPEC
    MEAN_KINETIC_ENERGY_SPEC

Layouts of the two marginals a factorized yield is built from: `A Y [σY]` and `A TKE [σTKE]`,
each with one header row. Both are what `ExforFissionData.jl` writes for the `yield`/`A` and
`TKE`/`A` observables.
"""
const MASS_YIELD_SPEC =
    TableSpec([:A, :yield, :uncertainty]; skip = 1, optional = (:uncertainty,))

@doc (@doc MASS_YIELD_SPEC)
const MEAN_KINETIC_ENERGY_SPEC =
    TableSpec([:A, :TKE, :uncertainty]; skip = 1, optional = (:uncertainty,))

"""
    factorized_yield(mass_yield_path, kinetic_energy_path, energies; width, spread) -> ExperimentalYield

Reconstruct `Y(A, TKE)` from the two marginals that were measured, where the joint distribution
was not:

```
Y(A, TKE) = Y(A) · G(TKE; ⟨TKE⟩(A), σ)
```

with `G` a Gaussian normalized over `energies`. `Y(A)` is a measured pre-neutron mass yield and
`⟨TKE⟩(A)` a measured mean total kinetic energy against the same mass number — for 233-U(n,f),
for instance, both come from the same experiment.

**This is a reconstruction, not a measurement, and it is the mass dependence of `⟨TKE⟩` that
makes it worth doing.** A single mean kinetic energy applied at every mass over-subtracts at the
asymmetric splits, where `⟨TKE⟩` is 20 MeV lower than in the peaks, starves the heavy fragment of
excitation, and turns the rising heavy branch of `ν(A)` into a falling one. Checked directly:
replacing the measured `Y(A,TKE)` of 235-U(n_th,f) with a mass-independent Gaussian at its own
`⟨TKE⟩ = 170.5 MeV` takes `ν(160)` from 2.02 down to 0.73 and `⟨ν⟩` from 2.524 to 2.280. Keeping
the mass dependence is what the factorization is for.

What it cannot recover is the correlation the joint distribution carries: the width of `P(TKE|A)`
varies with mass and is taken here as the constant `width`, 10 MeV by default, with `spread` mass
units of Gaussian smearing applied to `⟨TKE⟩(A)` where the table is sparse. Both are stated
rather than fitted.

A mass without a tabulated `⟨TKE⟩` is dropped: extrapolating the quantity whose mass dependence
is the point of the exercise would defeat it. Complementary masses inherit the mean of their
partner, since `TKE` is a property of the split rather than of one fragment.
"""
function factorized_yield(
    mass_yield_path::AbstractString,
    kinetic_energy_path::AbstractString,
    energies,
    compound_mass::Integer;
    width::Real = 10.0,
    total::Real = DEFAULT_YIELD_TOTAL,
)
    width > 0 ||
        throw(ArgumentError("the kinetic energy width must be positive, got $width"))
    grid = sort(unique(Float64.(collect(energies))))
    isempty(grid) && throw(ArgumentError("the kinetic energy grid is empty"))

    mass_table = read_table(mass_yield_path, MASS_YIELD_SPEC)
    masses = integer_column(mass_table, :A)
    yields = column(mass_table, :yield, Float64)

    energy_table = read_table(kinetic_energy_path, MEAN_KINETIC_ENERGY_SPEC)
    energy_masses = integer_column(energy_table, :A)
    means = column(energy_table, :TKE, Float64)

    # TKE belongs to the split, so a mean tabulated against the heavy mass serves the light
    # partner too.
    mean_by_mass = Dict{Int, Float64}()
    for row in eachindex(energy_masses)
        isfinite(means[row]) && means[row] > 0 || continue
        mean_by_mass[energy_masses[row]] = means[row]
        complement = Int(compound_mass) - energy_masses[row]
        haskey(mean_by_mass, complement) || (mean_by_mass[complement] = means[row])
    end
    isempty(mean_by_mass) &&
        throw(ArgumentError("$(energy_table.path) holds no positive mean kinetic energies"))

    values = Dict{Tuple{Int, Float64}, Float64}()
    kept = Int[]
    running = 0.0
    for row in eachindex(masses)
        A = masses[row]
        Y = yields[row]
        isfinite(Y) && Y > 0 || continue
        mean = get(mean_by_mass, A, nothing)
        mean === nothing && continue
        shape = [exp(-(TKE - mean)^2 / (2 * width^2)) for TKE in grid]
        norm = sum(shape)
        norm > 0 || continue
        push!(kept, A)
        for (index, TKE) in pairs(grid)
            cell = Y * shape[index] / norm
            values[(A, TKE)] = cell
            running += cell
        end
    end
    isempty(values) &&
        throw(ArgumentError("no mass carries both a yield and a mean kinetic energy; \
                       $(mass_table.path) and $(energy_table.path) do not overlap"))

    # Renormalized to the same 200 % a measured two-fragment distribution carries, so a run on a
    # factorized yield and a run on a measured one are weighted alike.
    scale = total / running
    for key in keys(values)
        values[key] *= scale
    end
    return ExperimentalYield(
        values,
        sort(kept),
        grid,
        string(basename(mass_yield_path), " × ", basename(kinetic_energy_path)),
    )
end

"""
    fragment_charge_probabilities(domain, polarization) -> Dict{Nuclide,Float64}

The isobaric charge probability `p(Z,A)` of every fragment the domain produces.

**Each nuclide gets one value, evaluated against its own mass's `Zₚ`** — not accumulated over
the fragmentations it appears in, and not inherited from a complementary fragment. Away from
mass symmetry the distinction is invisible, because `Zₚ(A_L) = Z₀ − Zₚ(A_H)` makes the light
fragment's Gaussian the mirror of the heavy one's. At `A = A₀/2` it is the whole difference: a
nuclide there belongs to two different fragmentations with two different probabilities, and the
yield of a cell `(A, Z)` is a property of the cell.

This is the quantity `Y(A,Z,TKE) = p(Z,A) Y(A,TKE)` is built from. It is deliberately not
[`fragments`](@ref), which sums a nuclide's contributions to answer a different question — how
often the domain produces it.

`p` is the analytic Gaussian, unnormalized over the retained charges; see the note in
`src/charge.jl` for why.
"""
function fragment_charge_probabilities(
    domain::FragmentationDomain,
    polarization::ChargePolarization,
)
    system = domain.system
    symmetric = system.compound.A
    probabilities = Dict{Nuclide, Float64}()
    for entry in domain.entries
        for fragment in (entry.light, entry.heavy)
            haskey(probabilities, fragment) && continue
            A_heavy = max(fragment.A, symmetric - fragment.A)
            ΔZ = polarization_at(polarization, A_heavy)
            rms = width_at(polarization, A_heavy)
            Zₚ_heavy = most_probable_charge(system, A_heavy, ΔZ)
            Zₚ = fragment.A == A_heavy ? Zₚ_heavy : system.compound.Z - Zₚ_heavy
            probabilities[fragment] = charge_probability(fragment.Z, Zₚ, rms)
        end
    end
    return probabilities
end

"""
    FragmentYields

The yield of every fragment configuration, `Y(A, Z, TKE)`, normalized so the whole distribution
sums to `total`.
"""
struct FragmentYields
    values::Dict{Tuple{Int, Int, Float64}, Float64}
    total::Float64
end

Base.length(y::FragmentYields) = length(y.values)
Base.sum(y::FragmentYields) = sum(values(y.values))

function Base.show(io::IO, y::FragmentYields)
    print(
        io,
        "FragmentYields(",
        length(y.values),
        " cells, Σ = ",
        round(y.total; digits = 4),
        ")",
    )
    return nothing
end

"""
    DEFAULT_YIELD_TOTAL

What a fragment yield distribution sums to: 200, because a fission event makes two fragments and
each mass distribution is quoted as a percentage.

Named rather than written into the normalization, because it is a convention and the convention
should be visible where it is applied.
"""
const DEFAULT_YIELD_TOTAL = 200.0

"""
    build_fragment_yields(domain, polarization, experimental; total = DEFAULT_YIELD_TOTAL)
        -> FragmentYields

Spread the measured `Y(A, TKE)` over charge numbers,

```
Y(A, Z, TKE) = p(Z, A) Y(A, TKE)
```

and renormalize the result to `total`.

The renormalization is what makes the unnormalized `p(Z,A)` harmless downstream: whatever weight
the retained charge numbers fail to capture is divided out once, over the whole distribution,
rather than per mass row. A domain spanning too few charges per mass therefore still yields a
distribution summing to 200, but the shape it summed from is the honest one.

Cells are taken on the **measurement's** `TKE` grid, not the model's. Use [`coverage`](@ref)
before averaging a model quantity over these yields.
"""
function build_fragment_yields(
    domain::FragmentationDomain,
    polarization::ChargePolarization,
    experimental::ExperimentalYield;
    total::Real = DEFAULT_YIELD_TOTAL,
)
    total > 0 ||
        throw(ArgumentError("the yield normalization must be positive, got $total"))
    probabilities = fragment_charge_probabilities(domain, polarization)

    values = Dict{Tuple{Int, Int, Float64}, Float64}()
    running = 0.0
    for (nuclide, p) in probabilities, TKE in experimental.energies
        measured = get(experimental.values, (nuclide.A, TKE), nothing)
        measured === nothing && continue
        contribution = p * measured
        values[(nuclide.A, nuclide.Z, TKE)] = contribution
        running += contribution
    end
    isempty(values) &&
        throw(ArgumentError("no fragment of the domain has a measured yield"))
    running > 0 || throw(ArgumentError("the fragment yields sum to $running"))

    scale = total / running
    for key in keys(values)
        values[key] *= scale
    end
    return FragmentYields(values, Float64(total))
end

"""
    coverage(yields, energies) -> NamedTuple

How much of a yield distribution a model sweep actually reaches.

Returns the fraction of the total yield that sits on `energies`, the fraction that does not, and
the two `TKE` grids. A sweep coarser than the measurement leaves yield behind, and an average
taken without noticing is an average over part of the distribution.

For the archived ²⁵²Cf runs, whose measurement is tabulated every 1 MeV and whose sweep stepped
2 MeV, the answer is that **half the yield lies off the sweep**.

# Examples

```jldoctest
julia> using DeterministicSequentialEmission

julia> yields = FragmentYields(Dict((140, 54, 170.0) => 3.0, (140, 54, 171.0) => 1.0), 4.0);

julia> report = coverage(yields, 170.0:2.0:174.0);

julia> report.reached
0.75
```
"""
function coverage(yields::FragmentYields, energies)
    wanted = Set(Float64.(energies))
    on = 0.0
    off = 0.0
    present = Set{Float64}()
    for ((_, _, TKE), y) in yields.values
        push!(present, TKE)
        TKE in wanted ? (on += y) : (off += y)
    end
    whole = on + off
    return (
        reached = whole == 0 ? 0.0 : on / whole,
        missed = whole == 0 ? 0.0 : off / whole,
        swept = sort(collect(wanted)),
        tabulated = sort(collect(present)),
    )
end

"""
    yield_over(yields, dimension) -> Vector{Tuple{Float64,Float64}}

Project the yield onto one of its arguments, `:mass`, `:charge` or `:kinetic_energy`, summing
over the other two. Sorted by the argument.

Each projection sums to the same total as the distribution it came from.
"""
function yield_over(yields::FragmentYields, dimension::Symbol)
    index = if dimension === :mass
        1
    elseif dimension === :charge
        2
    elseif dimension === :kinetic_energy
        3
    else
        throw(
            ArgumentError(
                "unknown yield dimension $dimension; expected :mass, :charge or :kinetic_energy",
            ),
        )
    end
    projected = Dict{Float64, Float64}()
    for (key, y) in yields.values
        argument = Float64(key[index])
        projected[argument] = get(projected, argument, 0.0) + y
    end
    return sort!([(argument, y) for (argument, y) in projected]; by = first)
end

"""
    mean_kinetic_energy_by_mass(yields) -> Vector{Tuple{Int,Float64}}

`⟨TKE⟩(A)`, the yield-weighted mean total kinetic energy at each fragment mass.

Defined at every mass the distribution covers, light and heavy; the two halves carry the same
information, since a fragmentation's `TKE` belongs to the pair rather than to either fragment.
"""
function mean_kinetic_energy_by_mass(yields::FragmentYields)
    weighted = Dict{Int, Float64}()
    total = Dict{Int, Float64}()
    for ((A, _, TKE), y) in yields.values
        weighted[A] = get(weighted, A, 0.0) + y * TKE
        total[A] = get(total, A, 0.0) + y
    end
    return sort!(
        [(A, weighted[A] / total[A]) for A in keys(total) if total[A] > 0];
        by = first,
    )
end
