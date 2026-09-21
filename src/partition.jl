"""
    SegmentedCurve(masses, values, source = "")

A piecewise-linear function of heavy-fragment mass, given as segment endpoints `(A_H, y)` in
ascending `A_H`.

This is how `E*_H/TXE (A_H)` and `R_T (A_H)` are supplied: a handful of points fitted elsewhere,
joined by straight lines. The alternative — the raw point-by-point values — is unusable as an
input, because it is too scattered and too incompletely covered in `A_H` to fit.

`source` records the file the curve was read from, so that a run can say which tabulation it
partitioned with rather than leaving that to whoever staged it.

Outside the tabulated range the first segment is extended below the first point, and the
function is undefined above the last. That asymmetry reproduces the prototype, whose published
curves always end at `A_H_max`.
"""
struct SegmentedCurve
    masses::Vector{Int}
    values::Vector{Float64}
    source::String

    function SegmentedCurve(
        masses::Vector{Int},
        values::Vector{Float64},
        source::AbstractString = "",
    )
        length(masses) == length(values) ||
            throw(ArgumentError("segment masses and values must have equal length"))
        length(masses) >= 2 ||
            throw(ArgumentError("a segmented curve needs at least two points"))
        issorted(masses) || throw(
            ArgumentError("segment points must be given in ascending mass, got $masses"),
        )
        allunique(masses) ||
            throw(ArgumentError("segment points must have distinct masses, got $masses"))
        return new(copy(masses), copy(values), String(source))
    end
end

"""
    SegmentedCurve(points::AbstractVector{<:Tuple})

Build from `(A_H, value)` pairs.

# Examples

```jldoctest
julia> curve = SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]);

julia> curve(130)
0.14

julia> round(curve(140); digits = 4)
0.3735
```
"""
function SegmentedCurve(
    points::AbstractVector{<:Tuple{Integer, Real}},
    source::AbstractString = "",
)
    return SegmentedCurve(
        [Int(p[1]) for p in points],
        [Float64(p[2]) for p in points],
        source,
    )
end

Base.length(curve::SegmentedCurve) = length(curve.masses)

function Base.show(io::IO, curve::SegmentedCurve)
    print(
        io,
        "SegmentedCurve(",
        length(curve),
        " points, A_H ∈ ",
        first(curve.masses),
        ":",
        last(curve.masses),
        ")",
    )
    return nothing
end

"""
    (curve)(A) -> Union{Float64,Nothing}

Evaluate the curve at heavy-fragment mass `A`, or `nothing` above the last point or where the
interpolated value is non-positive.
"""
function (curve::SegmentedCurve)(A::Integer)
    for index in 2:length(curve)
        if curve.masses[index] >= A
            x₀ = curve.masses[index - 1]
            x₁ = curve.masses[index]
            y₀ = curve.values[index - 1]
            y₁ = curve.values[index]
            slope = (y₁ - y₀) / (x₁ - x₀)
            y = y₀ + slope * (A - x₀)
            return y > 0 ? y : nothing
        end
    end
    return nothing
end

"""
    PartitionMethod

How the total excitation energy is divided between the complementary fragments at full
acceleration.

The partition happens at scission, but not every method models it there: two of the three
divide `TXE` directly at full acceleration using a curve fitted elsewhere.

Whatever the method, the result carries no uncertainty. The curves quote none, so `E*_L` and
`E*_H` inherit none — which is why uncertainty propagation stops at `TXE`.
"""
abstract type PartitionMethod end

"""
    ExcitationEnergyRatio(curve)

Partition by `E*_H/TXE` as a function of `A_H`, supplied as a segmented curve.

The ratio is normally obtained from modelling at scission in a separate code, then handed to the
prompt-emission calculation as a tabulated curve. No manifest describes one, because no package
of this toolchain produces it — which is why this method takes a file path where
[`TemperatureRatio`](@ref) takes a run record.
"""
struct ExcitationEnergyRatio <: PartitionMethod
    curve::SegmentedCurve

    function ExcitationEnergyRatio(curve::SegmentedCurve)
        for (index, value) in pairs(curve.values)
            0 < value < 1 || throw(
                ArgumentError(
                    "E*_H/TXE is a fraction of the total excitation energy and must lie in " *
                    "(0, 1); point $index, at A_H = $(curve.masses[index]), is $value",
                ),
            )
        end
        return new(curve)
    end
end

"""
    TemperatureRatio(curve, model)

Partition by the ratio of fragment temperatures at full acceleration, `R_T = T_L/T_H`, as a
function of `A_H`.

With both fragments in the Fermi-gas regime, `E*_{L,H} = a_{L,H} T_{L,H}²` gives

```
E*_H/TXE = 1 / (1 + (a_L/a_H) R_T²)
```

so the temperature ratio and the level density parameters together fix the split. `R_T` is the
input parameter by which several prompt-emission codes partition `TXE` directly.
"""
struct TemperatureRatio{M <: LevelDensityModel} <: PartitionMethod
    curve::SegmentedCurve
    model::M

    function TemperatureRatio(
        curve::SegmentedCurve,
        model::M,
    ) where {M <: LevelDensityModel}
        for (index, value) in pairs(curve.values)
            value > 0 || throw(
                ArgumentError(
                    "R_T = T_L/T_H is a ratio of temperatures and must be positive; point " *
                    "$index, at A_H = $(curve.masses[index]), is $value",
                ),
            )
        end
        return new{M}(curve, model)
    end
end

"""
    ScissionModelling(deformation, model)

Partition by modelling at scission, using extra-deformation energies.

Two assumptions fix the split. The fragments are far more deformed at scission than at full
acceleration, and the difference

```
ΔE_def = E_def(β₂^sc) − E_def(β₂^a.t.)
```

reappears as excitation energy of the fully accelerated fragment. And the complementary
fragments are in thermal equilibrium at scission with Fermi-gas level densities, so
`T_L^sc = T_H^sc` gives `E_sc^L/E_sc^H = a_L/a_H`. Together:

```
E_sc^{L,H} = [TXE − (ΔE_def^L + ΔE_def^H)] / (1 + a_{H,L}/a_{L,H})
E*_{L,H}   = E_sc^{L,H} + ΔE_def^{L,H}
```
"""
struct ScissionModelling{M <: LevelDensityModel} <: PartitionMethod
    deformation::Dict{Nuclide, Float64}
    model::M
end

"""
    TEMPERATURE_RATIO_CURVE_SPEC

Layout of a tabulated `R_T(A_H)`: one header row, then `A_H,R_T` and an optional
`R_T_uncertainty`, comma-delimited. This is what a temperature-ratio run writes.

The uncertainty is read past and discarded: the partition carries none downstream, for the reason
[`PartitionMethod`](@ref) gives.
"""
const TEMPERATURE_RATIO_CURVE_SPEC = TableSpec(
    [:A_H, :R_T, :R_T_uncertainty];
    delimiter = ',',
    skip = 1,
    optional = (:R_T_uncertainty,),
)

"""
    EXCITATION_RATIO_CURVE_SPEC

Layout of a tabulated `E*_H/TXE (A_H)`: one header row, then `A_H,excitation_ratio` and an
optional `excitation_ratio_uncertainty`, comma-delimited.

A separate spec from [`TEMPERATURE_RATIO_CURVE_SPEC`](@ref) rather than one shared spec naming its
ordinate `value`, because the two tabulate different quantities and a column named for neither
would put the reader back to inferring from the file name what the file holds.
"""
const EXCITATION_RATIO_CURVE_SPEC = TableSpec(
    [:A_H, :excitation_ratio, :excitation_ratio_uncertainty];
    delimiter = ',',
    skip = 1,
    optional = (:excitation_ratio_uncertainty,),
)

"""
    read_segmented_curve(path; spec = TEMPERATURE_RATIO_CURVE_SPEC) -> SegmentedCurve

Read a curve tabulated against heavy-fragment mass.

**Columns are taken by position**, not by header text: the abscissa is the first column of `spec`
and the ordinate the second, whatever the file calls them. The spec names them only so that the
code reads like the table it is reading; nothing here may be changed to a lookup by header text,
because that would couple this package to the header a producer happens to write.

Written for the output of a temperature-ratio extraction, which takes `R_T = T_L/T_H` from
experimental `ν(A)` and writes it at **every integer `A_H`** of its fit range. Fed whole to
[`SegmentedCurve`](@ref) that is exact rather than an approximation: every abscissa the model
evaluates is itself a tabulated point, so no interpolation happens. Handing over only the
breakpoints of that extraction's own piecewise-linear fit would not be equivalent — it fits the
multiplicity ratio `r_ν`, and `R_T = √[(1−r_ν)/(R_a r_ν)]` is not piecewise linear where `r_ν`
is, because the level density parameter ratio `R_a(A_H)` carries the shell structure. That is also
why the pivots file is never the file to read; see [`temperature_ratio_path`](@ref).

Two things to expect rather than diagnose later. Each curve spans its own fit range, not the
configured `A_H` range, so a set that stops short leaves [`partition`](@ref) returning `nothing`
above its last point — correct, since the alternative is extrapolating a fit past its last
datum. And the ratio averaging must match: `R_T` produced for this consumer should come from a
run with `ratio_averaging = "mean_of_ratios"`, because the partition re-expands it with
per-`(A_H, Z_H)` level density parameters and reduces over `Z` afterwards, which inverts a
`Z`-averaged multiplicity ratio with `⟨a_L/a_H⟩` rather than `⟨a_L⟩/⟨a_H⟩`. That setting is
recorded in the run identifier the file name carries.
"""
function read_segmented_curve(
    path::AbstractString;
    spec::TableSpec = TEMPERATURE_RATIO_CURVE_SPEC,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, first(spec.columns))
    values = column(table, spec.columns[2], Float64)
    return SegmentedCurve(masses, values, path)
end

"""
    segment_coverage(method, heavy_masses) -> NamedTuple

Which heavy masses of a sweep the curve actually covers.

Returns the masses it does not reach and the fraction it does. A curve is a fit over the range
its data supported, and nothing ties that range to the domain being swept: the
prototype anchors its last segment point at `A_H_max` by construction, but a curve read from a
file carries its own range and may stop short.

Worth checking rather than discovering, because [`partition`](@ref) returns `nothing` above the
last point and that exclusion is indistinguishable from a physical one — a fragmentation the
model could not solve and a fragmentation no curve reaches look identical downstream.

# Examples

```jldoctest
julia> method = ExcitationEnergyRatio(SegmentedCurve([(126, 0.5), (150, 0.6)]));

julia> report = segment_coverage(method, 126:174);

julia> report.covered
0.5102040816326531

julia> first(report.uncovered), last(report.uncovered)
(151, 174)
```
"""
function segment_coverage(
    method::Union{ExcitationEnergyRatio, TemperatureRatio},
    heavy_masses,
)
    return _coverage(A -> method.curve(A) !== nothing, heavy_masses)
end

# No curve to cover the range here: a heavy mass is reached where the extra-deformation table
# holds a nuclide of that mass.
function segment_coverage(method::ScissionModelling, heavy_masses)
    tabulated = Set(nuclide.A for nuclide in keys(method.deformation))
    return _coverage(A -> A in tabulated, heavy_masses)
end

function _coverage(is_covered, heavy_masses)
    uncovered = Int[]
    total = 0
    for A in heavy_masses
        total += 1
        is_covered(A) || push!(uncovered, A)
    end
    return (
        covered = total == 0 ? 0.0 : (total - length(uncovered)) / total,
        uncovered = uncovered,
    )
end

"""
    EXTRA_DEFORMATION_SPEC

Layout of the extra-deformation energy files: one header row, then `A Z ΔE_def`.
"""
const EXTRA_DEFORMATION_SPEC = TableSpec([:A, :Z, :E_def]; skip = 1)

"""
    read_extra_deformation(path; spec = EXTRA_DEFORMATION_SPEC) -> Dict{Nuclide,Float64}

Read extra-deformation energies in MeV, indexed by nuclide.

**Columns are taken by position**, not by header text: mass, charge, `ΔE_def`.
"""
function read_extra_deformation(
    path::AbstractString;
    spec::TableSpec = EXTRA_DEFORMATION_SPEC,
)
    table = read_delimited_table(path, spec)
    masses = integer_column(table, :A)
    charges = integer_column(table, :Z)
    energies = column(table, :E_def, Float64)
    deformation = Dict{Nuclide, Float64}()
    for row in eachindex(masses)
        deformation[Nuclide(charges[row], masses[row])] = energies[row]
    end
    return deformation
end

"""
    ExcitationEnergies

The excitation energies of the two fully accelerated fragments, in MeV. By construction
`light + heavy == TXE`.
"""
struct ExcitationEnergies
    light::Float64
    heavy::Float64

    function ExcitationEnergies(light::Real, heavy::Real)
        return new(Float64(light), Float64(heavy))
    end
end

Base.sum(energies::ExcitationEnergies) = energies.light + energies.heavy

function Base.show(io::IO, energies::ExcitationEnergies)
    print(
        io,
        "ExcitationEnergies(L=",
        round(energies.light; digits = 4),
        ", H=",
        round(energies.heavy; digits = 4),
        " MeV)",
    )
    return nothing
end

"""
    partition(method, fragmentation, TXE) -> Union{ExcitationEnergies,Nothing}

Divide `TXE` between the complementary fragments, or `nothing` when the method cannot be
evaluated for this fragmentation — a mass outside the curve, a missing level
density parameter, an absent deformation energy, or a division that leaves either fragment
with negative excitation energy.

A split into two copies of the same nuclide divides `TXE` evenly, and does so before the
partition method is consulted — see [`_self_complementary`](@ref).
"""
function partition(method::ExcitationEnergyRatio, fragmentation::Fragmentation, TXE::Real)
    _self_complementary(fragmentation) && return _split(TXE, 0.5)
    ratio = method.curve(fragmentation.heavy.A)
    (ratio === nothing || ratio > 1) && return nothing
    return _split(TXE, ratio)
end

function partition(method::TemperatureRatio, fragmentation::Fragmentation, TXE::Real)
    _self_complementary(fragmentation) && return _split(TXE, 0.5)
    R_T = method.curve(fragmentation.heavy.A)
    R_T === nothing && return nothing
    a_light = level_density_parameter(method.model, fragmentation.light)
    a_heavy = level_density_parameter(method.model, fragmentation.heavy)
    (a_light === nothing || a_heavy === nothing) && return nothing
    ratio = 1 / (1 + (a_light / a_heavy) * R_T^2)
    return _split(TXE, ratio)
end

function partition(method::ScissionModelling, fragmentation::Fragmentation, TXE::Real)
    ΔE_light = get(method.deformation, fragmentation.light, nothing)
    ΔE_heavy = get(method.deformation, fragmentation.heavy, nothing)
    (ΔE_light === nothing || ΔE_heavy === nothing) && return nothing
    a_light = level_density_parameter(method.model, fragmentation.light)
    a_heavy = level_density_parameter(method.model, fragmentation.heavy)
    (a_light === nothing || a_heavy === nothing) && return nothing

    at_scission = TXE - (ΔE_light + ΔE_heavy)
    at_scission >= 0 || return nothing
    r = a_light / a_heavy
    heavy = at_scission / (1 + r) + ΔE_heavy
    light = r * at_scission / (1 + r) + ΔE_light
    (light >= 0 && heavy >= 0) || return nothing
    return ExcitationEnergies(light, heavy)
end

"""
    _self_complementary(fragmentation) -> Bool

Whether the two fragments are the same nuclide, `A_H == A_L` **and** `Z_H == Z_L`.

Such a split must divide `TXE` evenly, and the reason is a symmetry rather than a convention:
the fragments are indistinguishable, so `a_L == a_H` and `T_L == T_H` identically, hence
`E*_L == E*_H` and `R_T == 1`. No curve in `A_H` can know this — `R_T(A_H)` and
`E*_H/TXE (A_H)` are fitted, and nothing makes them pass through `1` and `0.5` at
`A₀/2`. Evaluated there with the ²⁵²Cf literature value `R_T = 1.2` the temperature-ratio
method returns `E*_H/TXE = 0.41`, dividing energy unequally between two copies of the
same nucleus.

Guarding on it restores the behaviour of the prototype, which special-cased
`A_L == A_H && Z_L == Z_H` in all three methods. `ScissionModelling` needs no guard and
does not carry one: identical fragments give `a_L/a_H == 1` and `ΔE_def^L == ΔE_def^H`, so it
divides evenly on its own arithmetic, and short-circuiting it would skip its check that the
energy remaining at scission is non-negative.

Both archived thesis curves happen to satisfy the constraint through their first
pivot — `R_T(126) = 1.0` and `E*_H/TXE (126) = 0.5` — which is why the reference runs are
unaffected either way.
"""
_self_complementary(fragmentation::Fragmentation) =
    fragmentation.heavy == fragmentation.light

function _split(TXE::Real, ratio::Real)
    heavy = TXE * ratio
    light = TXE - heavy
    (light >= 0 && heavy >= 0) || return nothing
    return ExcitationEnergies(light, heavy)
end
