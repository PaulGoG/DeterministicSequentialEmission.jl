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
partition method is consulted.
"""
function partition(method::ExcitationEnergyRatio, fragmentation::Fragmentation, TXE::Real)
    is_self_complementary(fragmentation) && return _split(TXE, 0.5)
    ratio = method.curve(fragmentation.heavy.A)
    (ratio === nothing || ratio > 1) && return nothing
    return _split(TXE, ratio)
end

function partition(method::TemperatureRatio, fragmentation::Fragmentation, TXE::Real)
    is_self_complementary(fragmentation) && return _split(TXE, 0.5)
    R_T = method.curve(fragmentation.heavy.A)
    R_T === nothing && return nothing
    a_light = level_density_parameter(method.model, fragmentation.light)
    a_heavy = level_density_parameter(method.model, fragmentation.heavy)
    (a_light === nothing || a_heavy === nothing) && return nothing
    ratio = heavy_excitation_fraction(R_T, a_light / a_heavy)
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

# A split into two copies of one nuclide divides TXE evenly by symmetry: `a_L == a_H` and
# `T_L == T_H`, hence `E*_L == E*_H` and `R_T == 1`. Fitted curves `R_T(A_H)` and `E*_H/TXE (A_H)`
# are not constrained to pass through 1 and 0.5 at `A₀/2`, so the two curve-based methods apply
# this guard. `ScissionModelling` does not need it: identical fragments give `a_L/a_H == 1` and
# `ΔE_def^L == ΔE_def^H`, so its arithmetic divides evenly.

function _split(TXE::Real, ratio::Real)
    heavy = TXE * ratio
    light = TXE - heavy
    (light >= 0 && heavy >= 0) || return nothing
    return ExcitationEnergies(light, heavy)
end
