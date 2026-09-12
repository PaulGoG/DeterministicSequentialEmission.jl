"""
    ChargePolarization

Charge polarization `ΔZ(A)` and width `rms(A)` of the isobaric charge distribution, indexed by
fragment mass number, with the values used where a mass is not tabulated.

The defaults `⟨ΔZ⟩ = -0.5` and `⟨rms⟩ = 0.6` are the conventional means applied when no `Zₚ`
systematics is available for the fission case under study. `ΔZ` is quoted for the **heavy**
fragment and carries the opposite sign for the light one.
"""
struct ChargePolarization
    polarization::Dict{Int, Float64}
    width::Dict{Int, Float64}
    default_polarization::Float64
    default_width::Float64
    source::String

    function ChargePolarization(
        polarization::Dict{Int, Float64},
        width::Dict{Int, Float64},
        default_polarization::Real,
        default_width::Real,
        source::AbstractString,
    )
        default_width > 0 ||
            throw(ArgumentError("the default rms must be positive, got $default_width"))
        return new(
            polarization,
            width,
            Float64(default_polarization),
            Float64(default_width),
            String(source),
        )
    end
end

"""
    CHARGE_POLARIZATION_SPEC

Layout of the `ΔZ(A)`, `rms(A)` files: one header row, then `A ΔZ rms`.
"""
const CHARGE_POLARIZATION_SPEC = TableSpec([:A, :polarization, :width]; skip = 1)

const DEFAULT_POLARIZATION = -0.5
const DEFAULT_WIDTH = 0.6

"""
    read_charge_polarization(path; spec, default_polarization, default_width) -> ChargePolarization

Read tabulated `ΔZ(A)` and `rms(A)`.
"""
function read_charge_polarization(
    path::AbstractString;
    spec::TableSpec = CHARGE_POLARIZATION_SPEC,
    default_polarization::Real = DEFAULT_POLARIZATION,
    default_width::Real = DEFAULT_WIDTH,
)
    table = read_table(path, spec)
    masses = column(table, :A, Int)
    polarization = column(table, :polarization, Float64)
    width = column(table, :width, Float64)

    polarizations = Dict{Int, Float64}()
    widths = Dict{Int, Float64}()
    for row in eachindex(masses)
        width[row] > 0 || throw(
            ArgumentError(
                "$(table.path) row $row: rms must be positive, got $(width[row]) for A = $(masses[row])",
            ),
        )
        polarizations[masses[row]] = polarization[row]
        widths[masses[row]] = width[row]
    end
    return ChargePolarization(
        polarizations,
        widths,
        Float64(default_polarization),
        Float64(default_width),
        String(path),
    )
end

"""
    mean_charge_polarization(; default_polarization, default_width) -> ChargePolarization

The fallback systematics with no tabulated values: every mass number takes the mean
`ΔZ` and `rms`.
"""
function mean_charge_polarization(;
    default_polarization::Real = DEFAULT_POLARIZATION,
    default_width::Real = DEFAULT_WIDTH,
)
    return ChargePolarization(
        Dict{Int, Float64}(),
        Dict{Int, Float64}(),
        Float64(default_polarization),
        Float64(default_width),
        "mean values",
    )
end

function Base.show(io::IO, polarization::ChargePolarization)
    print(
        io,
        "ChargePolarization(",
        length(polarization.polarization),
        " masses from ",
        basename(polarization.source),
        ")",
    )
    return nothing
end

"""
    polarization_at(polarization, A) -> Float64

Charge polarization `ΔZ` for heavy-fragment mass `A`, falling back to the mean value.
"""
polarization_at(polarization::ChargePolarization, A::Integer) =
    get(polarization.polarization, Int(A), polarization.default_polarization)

"""
    width_at(polarization, A) -> Float64

Width `rms` of the isobaric charge distribution at mass `A`, falling back to the mean value.
"""
width_at(polarization::ChargePolarization, A::Integer) =
    get(polarization.width, Int(A), polarization.default_width)

"""
    unchanged_charge_distribution(system, A) -> Float64

`Z_UCD(A) = A Z₀ / A₀`, the charge a fragment would carry if the compound nucleus split
without charge rearrangement.
"""
function unchanged_charge_distribution(system::FissioningSystem, A::Integer)
    return A * system.compound.Z / system.compound.A
end

"""
    most_probable_charge(system, A, ΔZ) -> Float64

`Zₚ(A) = Z_UCD(A) + ΔZ(A)`, the centre of the isobaric charge distribution.
"""
function most_probable_charge(system::FissioningSystem, A::Integer, ΔZ::Real)
    return unchanged_charge_distribution(system, A) + ΔZ
end

"""
    charge_probability(Z, Zₚ, rms) -> Float64

The isobaric charge distribution, a Gaussian of width `rms` centred on `Zₚ`:

```
p(Z, A) = exp(−(Z − Zₚ)² / 2 rms²) / (√(2π) rms)
```
"""
function charge_probability(Z::Real, Zₚ::Real, rms::Real)
    rms > 0 || throw(ArgumentError("rms must be positive, got $rms"))
    return exp(-(Z - Zₚ)^2 / (2 * rms^2)) / (sqrt(2π) * rms)
end

"""
    charge_numbers(Zₚ, count) -> UnitRange{Int}

The `count` integer charge numbers nearest `Zₚ`, centred on it.

`count` must be odd so that the set is symmetric about the rounded most probable charge; an
even count would have to break the tie in one direction and silently bias the domain.

Rounding follows Julia's default, ties to even, so a `Zₚ` landing exactly halfway between two
integers — which happens at symmetric fission under the mean polarization, where
`Zₚ = Z₀/2 − 0.5` — resolves downwards for an even neighbour and upwards for an odd one.
"""
function charge_numbers(Zₚ::Real, count::Integer)
    isodd(count) ||
        throw(ArgumentError("the number of charges per mass must be odd, got $count"))
    count > 0 || throw(ArgumentError("the number of charges per mass must be positive"))
    centre = round(Int, Zₚ)
    half = (count - 1) ÷ 2
    return (centre - half):(centre + half)
end
