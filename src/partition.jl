"""
    SegmentedCurve(masses, values, source = "")

A piecewise-linear function of heavy-fragment mass, given as segment endpoints `(A_H, y)` in
ascending `A_H`. `E*_H/TXE (A_H)` and `R_T (A_H)` are supplied in this form. `source` records
the file the curve was read from and enters the run metadata. Outside the tabulated range the
function evaluates to `nothing` on both sides, so a fit is never extrapolated.
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

Evaluate the curve at heavy-fragment mass `A`, or `nothing` outside the tabulated range or where
the interpolated value is non-positive.
"""
function (curve::SegmentedCurve)(A::Integer)
    A < first(curve.masses) && return nothing
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
acceleration. `ScissionModelling` models the division at scission; the two curve-based
methods divide `TXE` directly at full acceleration. The curves quote no uncertainty, so `E*_L`
and `E*_H` carry none and uncertainty propagation stops at `TXE`.
"""
abstract type PartitionMethod end

"""
    ExcitationEnergyRatio(curve)

Partition by `E*_H/TXE` as a function of `A_H`, supplied as a segmented curve. The ratio is
obtained from modelling at scission in a separate code and supplied as a tabulated file,
without a run manifest; [`TemperatureRatio`](@ref), by contrast, is resolved through a run
record. Every point must lie in `(0, 1)`.
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

so the temperature ratio and the level density parameters together fix the split. Every point
of the curve must be positive.
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
`R_T_uncertainty`, comma-delimited, as written by a temperature-ratio run. The uncertainty is
read and discarded, since the partition carries none downstream; see [`PartitionMethod`](@ref).
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
optional `excitation_ratio_uncertainty`, comma-delimited. Counterpart of
[`TEMPERATURE_RATIO_CURVE_SPEC`](@ref) for the excitation-energy ratio.
"""
const EXCITATION_RATIO_CURVE_SPEC = TableSpec(
    [:A_H, :excitation_ratio, :excitation_ratio_uncertainty];
    delimiter = ',',
    skip = 1,
    optional = (:excitation_ratio_uncertainty,),
)

"""
    read_segmented_curve(path; spec = TEMPERATURE_RATIO_CURVE_SPEC) -> SegmentedCurve

Read a curve tabulated against heavy-fragment mass. Columns are taken by position, not by
header text: the abscissa is the first column of `spec` and the ordinate the second.

The format is that of a temperature-ratio extraction, which tabulates `R_T = T_L/T_H` at every
integer `A_H` of its fit range, so [`SegmentedCurve`](@ref) evaluates only tabulated points and
no interpolation occurs. The breakpoints of the extraction's piecewise-linear fit of `r_ν` are
not equivalent, because `R_T = √[(1−r_ν)/(R_a r_ν)]` inherits the shell structure of
`R_a(A_H)`; see [`temperature_ratio_path`](@ref).

Each curve spans its own fit range, and [`partition`](@ref) returns `nothing` outside it. `R_T`
should come from a run with `ratio_averaging = "mean_of_ratios"`, since the partition applies
per-`(A_H, Z_H)` level density parameters and reduces over `Z` afterwards, which corresponds to
`⟨a_L/a_H⟩` rather than `⟨a_L⟩/⟨a_H⟩`; the setting is recorded in the run identifier the file
name carries.
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

The heavy masses of a sweep covered by the partition input: `covered`, the fraction reached,
and `uncovered`, the masses not reached. A curve carries its own fit range, independent of the
swept domain. [`partition`](@ref) returns `nothing` outside that range, which is
indistinguishable downstream from a physically excluded fragmentation.

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

# Scission modelling has no curve: a heavy mass is covered where the extra-deformation table
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

Read extra-deformation energies in MeV, indexed by nuclide. Columns are taken by position, not
by header text: mass, charge, `ΔE_def`.
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

Whether the two fragments are the same nuclide, `A_H == A_L` and `Z_H == Z_L`.

Such a split divides `TXE` evenly by symmetry: the fragments are indistinguishable, so
`a_L == a_H` and `T_L == T_H`, hence `E*_L == E*_H` and `R_T == 1`. Fitted curves `R_T(A_H)`
and `E*_H/TXE (A_H)` are not constrained to pass through `1` and `0.5` at `A₀/2`, so the two
curve-based methods apply this guard. `ScissionModelling` does not: identical fragments give
`a_L/a_H == 1` and `ΔE_def^L == ΔE_def^H`, so its arithmetic divides evenly and its check that
the energy remaining at scission is non-negative still applies.
"""
_self_complementary(fragmentation::Fragmentation) =
    fragmentation.heavy == fragmentation.light

function _split(TXE::Real, ratio::Real)
    heavy = TXE * ratio
    light = TXE - heavy
    (light >= 0 && heavy >= 0) || return nothing
    return ExcitationEnergies(light, heavy)
end
