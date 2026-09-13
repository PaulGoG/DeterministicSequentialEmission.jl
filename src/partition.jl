"""
    JointedSegments(points)

A piecewise-linear function of heavy-fragment mass, given as jointed segment endpoints
`(A_H, y)` in ascending `A_H`.

This is how `E*_H/TXE (A_H)` and `R_T (A_H)` are supplied: a handful of points fitted
elsewhere, joined by straight lines. The alternative — the raw point-by-point values — is
unusable as an input, because it is too scattered and too incompletely covered in `A_H` to
parameterize.

Outside the tabulated range the first segment is extended below the first point, and the
function is undefined above the last. That asymmetry reproduces the prototype, whose
published parameterizations always end at `A_H_max`.
"""
struct JointedSegments
    masses::Vector{Int}
    values::Vector{Float64}

    function JointedSegments(masses::Vector{Int}, values::Vector{Float64})
        length(masses) == length(values) ||
            throw(ArgumentError("segment masses and values must have equal length"))
        length(masses) >= 2 || throw(
            ArgumentError("a jointed-segment parameterization needs at least two points"),
        )
        issorted(masses) || throw(
            ArgumentError("segment points must be given in ascending mass, got $masses"),
        )
        allunique(masses) ||
            throw(ArgumentError("segment points must have distinct masses, got $masses"))
        return new(copy(masses), copy(values))
    end
end

"""
    JointedSegments(points::AbstractVector{<:Tuple})

Build from `(A_H, value)` pairs.

# Examples

```jldoctest
julia> segments = JointedSegments([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]);

julia> segments(130)
0.14

julia> round(segments(140); digits = 4)
0.3735
```
"""
function JointedSegments(points::AbstractVector{<:Tuple{Integer, Real}})
    return JointedSegments([Int(p[1]) for p in points], [Float64(p[2]) for p in points])
end

Base.length(segments::JointedSegments) = length(segments.masses)

function Base.show(io::IO, segments::JointedSegments)
    print(
        io,
        "JointedSegments(",
        length(segments),
        " points, A_H ∈ ",
        first(segments.masses),
        ":",
        last(segments.masses),
        ")",
    )
    return nothing
end

"""
    (segments)(A) -> Union{Float64,Nothing}

Evaluate the parameterization at heavy-fragment mass `A`, or `nothing` above the last point
or where the interpolated value is non-positive.
"""
function (segments::JointedSegments)(A::Integer)
    for index in 2:length(segments)
        if segments.masses[index] >= A
            x₀ = segments.masses[index - 1]
            x₁ = segments.masses[index]
            y₀ = segments.values[index - 1]
            y₁ = segments.values[index]
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
divide `TXE` directly at full acceleration using a parameterization obtained elsewhere.

Whatever the method, the result carries no uncertainty. The parameterizations quote none, so
`E*_L` and `E*_H` inherit none — which is why uncertainty propagation stops at `TXE`.
"""
abstract type PartitionMethod end

"""
    ExcitationEnergyRatio(segments)

Partition by `E*_H/TXE` as a function of `A_H`, supplied as jointed segments.

The ratio is normally obtained from modelling at scission in a separate code, then handed to
the prompt-emission calculation as a parameterization.
"""
struct ExcitationEnergyRatio <: PartitionMethod
    segments::JointedSegments
end

"""
    TemperatureRatio(segments, model)

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
    segments::JointedSegments
    model::M
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
    EXTRA_DEFORMATION_SPEC

Layout of the extra-deformation energy files: one header row, then `A Z ΔE_def`.
"""
const EXTRA_DEFORMATION_SPEC = TableSpec([:A, :Z, :energy]; skip = 1)

"""
    read_extra_deformation(path; spec = EXTRA_DEFORMATION_SPEC) -> Dict{Nuclide,Float64}

Read extra-deformation energies in MeV, indexed by nuclide.
"""
function read_extra_deformation(
    path::AbstractString;
    spec::TableSpec = EXTRA_DEFORMATION_SPEC,
)
    table = read_table(path, spec)
    masses = integer_column(table, :A)
    charges = integer_column(table, :Z)
    energies = column(table, :energy, Float64)
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
evaluated for this fragmentation — a mass outside the parameterization, a missing level
density parameter, an absent deformation energy, or a division that leaves either fragment
with negative excitation energy.

A split into two copies of the same nuclide divides `TXE` evenly, and does so before the
prescription is consulted — see [`_self_complementary`](@ref).
"""
function partition(method::ExcitationEnergyRatio, fragmentation::Fragmentation, TXE::Real)
    _self_complementary(fragmentation) && return _split(TXE, 0.5)
    ratio = method.segments(fragmentation.heavy.A)
    (ratio === nothing || ratio > 1) && return nothing
    return _split(TXE, ratio)
end

function partition(method::TemperatureRatio, fragmentation::Fragmentation, TXE::Real)
    _self_complementary(fragmentation) && return _split(TXE, 0.5)
    R_T = method.segments(fragmentation.heavy.A)
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
`E*_L == E*_H` and `R_T == 1`. No parameterization of `A_H` can know this — `R_T(A_H)` and
`E*_H/TXE (A_H)` are fitted curves, and nothing makes them pass through `1` and `0.5` at
`A₀/2`. Evaluated there with the ²⁵²Cf literature value `R_T = 1.2` the temperature-ratio
prescription returns `E*_H/TXE = 0.41`, dividing energy unequally between two copies of the
same nucleus.

Guarding on it restores the behaviour of the prototype, which special-cased
`A_L == A_H && Z_L == Z_H` in all three prescriptions. `ScissionModelling` needs no guard and
does not carry one: identical fragments give `a_L/a_H == 1` and `ΔE_def^L == ΔE_def^H`, so it
divides evenly on its own arithmetic, and short-circuiting it would skip its check that the
energy remaining at scission is non-negative.

Both archived thesis parameterizations happen to satisfy the constraint through their first
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
