"""
    WAHL_LOW_ENERGY

Coefficients of the `Zₚ` model of A. C. Wahl, *Systematics of Fission-Product Yields*, LA-13928
(2002), Table 2, for the low-energy branch `PE ≤ 8 MeV`. Each entry is `(P₁ … P₅)` of eq. (17),

```
Par = P₁ + P₂[Z_F − 92] + P₃[A_F − 236] + P₄[PE − 6.551] + P₅[A_F − 236]²
```

The even-odd factors `F_Z` and `F_N` of that model are not carried: they modulate the charge
distribution about `Zₚ` and this package treats the distribution as a plain Gaussian, so only the
centroid `ΔZ` and the dispersion `σ_Z` enter. Adding them would change `p(Z,A)` at every mass,
not only where they are tabulated, and is a separate decision from supplying a `Zₚ` systematics.
"""
const WAHL_LOW_ENERGY = (
    σZ140 = (0.566, 0.0, 0.0064, 0.0109, 0.0),
    ΔZ140 = (-0.487, 0.0, 0.0180, 0.0, -0.00203),
    σZSL = (-0.0038, 0.0, 0.0, 0.0, 0.0),
    ΔZSL = (-0.0080, 0.0, 0.0, 0.0, 0.0),
    SL50 = (0.191, 0.0, -0.0076, 0.0, 0.0),
    σZ50 = (0.356, 0.060, 0.0, 0.0, 0.0),
    ΔZmax = (0.699, 0.0, 0.0, 0.0, 0.0),
    σZSLW = (-0.045, 0.0094, 0.0, 0.0, 0.0),
    ΔZSLW = (0.0, -0.0045, 0.0, 0.0, 0.0),
)

"""
    WAHL_VALIDITY

The range of fissioning nuclides the systematics was fitted over: `Z_F ∈ [90, 98]`,
`A_F ∈ [230, 252]`, and, for the low-energy branch used here, `PE ≤ 8 MeV`.
"""
const WAHL_VALIDITY = (charge = 90:98, mass = 230:252, excitation = 8.0)

"""
    WahlSystematics

The `Zₚ` model evaluated for one fissioning nucleus: the nine parameters of eq. (17) and the
mass-number boundaries that divide `A'` into wing, peak and near-symmetry regions.

`A'` is the pre-neutron fragment mass, which is the mass number this package indexes fragments
by, so `ΔZ` and `σ_Z` transfer without a shift. Both are quoted for the heavy fragment, the
convention of [`ChargePolarization`](@ref) and of eqs. (9c) and (9d) of the report.

The boundaries follow eq. (10). `Ba` and `Bb` are the ends of the steep interval where `Zₚ`
crosses the `Z = 50` shell closure; the report prints `Ba = A'_max` and `Bb = A_F − A'_max`,
which puts the heavy-side crossing at the lower end. Figs. 18b, 19b and 20b place the labels the
other way round, and only that reading is continuous with the peak regions on either side, so it
is the one used. `Bb` then lands within a mass unit of where the evaluated tables put the maximum
of `ΔZ`, which is what it is defined to be; the constructor records the measurement.

**A tabulated `ΔZ(A)`, `rms(A)` evaluated for the reaction at hand is the better object and takes
precedence wherever one exists.** Fig. 18 quotes reduced `χ²` of 7.9 for these systematics against
2.9 for a per-reaction least-squares fit of the same model, and the evaluated tables are the
latter. This is the fallback for a system that has none, in place of a single mean `ΔZ` and `rms`
applied at every mass.
"""
struct WahlSystematics
    Z_F::Int
    A_F::Int
    PE::Float64
    σZ140::Float64
    ΔZ140::Float64
    σZSL::Float64
    ΔZSL::Float64
    SL50::Float64
    σZ50::Float64
    ΔZmax::Float64
    σZSLW::Float64
    ΔZSLW::Float64
    B2::Float64
    B4::Float64
    B5::Float64
    B6::Float64
    Ba::Float64
    Bb::Float64
end

_wahl_parameter(coefficients, Z_F::Integer, A_F::Integer, PE::Real) =
    coefficients[1] +
    coefficients[2] * (Z_F - 92) +
    coefficients[3] * (A_F - 236) +
    coefficients[4] * (PE - 6.551) +
    coefficients[5] * (A_F - 236)^2

"""
    WahlSystematics(Z_F, A_F, PE) -> WahlSystematics

Evaluate the `Zₚ` model for a fissioning nucleus `(Z_F, A_F)` at precursor excitation energy
`PE` in MeV — zero for spontaneous fission, `Sₙ + Eₙ` for neutron-induced fission.

Throws for a nucleus or an excitation energy outside the range the systematics was fitted over,
[`WAHL_VALIDITY`](@ref): extrapolating a set of empirical polynomials past their data is not a
degradation that announces itself in the numbers.
"""
function WahlSystematics(Z_F::Integer, A_F::Integer, PE::Real)
    Z_F in WAHL_VALIDITY.charge || throw(
        ArgumentError("the Zₚ systematics is fitted for Z_F in $(WAHL_VALIDITY.charge), \
                       got $Z_F"),
    )
    A_F in WAHL_VALIDITY.mass || throw(
        ArgumentError("the Zₚ systematics is fitted for A_F in $(WAHL_VALIDITY.mass), \
                       got $A_F"),
    )
    (isfinite(PE) && PE >= 0) || throw(
        ArgumentError("the precursor excitation energy must be finite and non-negative, \
                             got $PE MeV"),
    )
    PE <= WAHL_VALIDITY.excitation ||
        throw(ArgumentError("the low-energy branch of the Zₚ systematics holds to \
                       PE = $(WAHL_VALIDITY.excitation) MeV, got $PE MeV"))

    evaluate(key) = _wahl_parameter(getproperty(WAHL_LOW_ENERGY, key), Z_F, A_F, PE)
    σZ140 = evaluate(:σZ140)
    ΔZ140 = evaluate(:ΔZ140)
    σZSL = evaluate(:σZSL)
    ΔZSL = evaluate(:ΔZSL)
    SL50 = evaluate(:SL50)
    σZ50 = evaluate(:σZ50)
    ΔZmax = evaluate(:ΔZmax)
    σZSLW = evaluate(:σZSLW)
    ΔZSLW = evaluate(:ΔZSLW)

    # Eq. (12): the A' at which Zₚ reaches 50, blended from two estimates.
    #
    # Only the second, AK2, satisfies the condition the quantity is named for — substituting it
    # into Zₚ(A') = A'(Z_F/A_F) + ΔZ_max returns exactly 50, while AK1 returns 49.27 for 236-U
    # and 48.99 for 240-Pu — which makes the blend look wrong on its face. It is not. `Bb` is by
    # construction the mass at which ΔZ attains ΔZ_max, and the evaluated tables say where that
    # is: 125 for 235-U, 124 for 239-Pu, 126 for 252-Cf. The blend gives 124.60, 124.04 and
    # 126.77, within 0.04 to 0.8 of each; AK2 alone gives 126.47, 125.87 and 126.77 and is worse
    # in the two cases where the two differ. Dropping AK1 also widens the near-symmetry region
    # away from the tabulated one and takes the largest ΔZ difference over the heavy branch from
    # 0.227 to 0.583 for 235-U. The published formula is right and is kept.
    F₁ = clamp((250 - A_F) / 14, 0.0, 1.0)
    A_max = F₁ * (50 * (A_F / Z_F) - ΔZmax / SL50) + (1 - F₁) * ((50 - ΔZmax) * (A_F / Z_F))

    B2 = 77 + 0.036 * (A_F - 236)
    B4 = (ΔZmax - ΔZ140 + A_max * SL50 + 140 * ΔZSL) / (SL50 + ΔZSL)
    return WahlSystematics(
        Int(Z_F),
        Int(A_F),
        Float64(PE),
        σZ140,
        ΔZ140,
        σZSL,
        ΔZSL,
        SL50,
        σZ50,
        ΔZmax,
        σZSLW,
        ΔZSLW,
        B2,
        B4,
        A_F - B2,
        A_F - 70.0,
        A_F - A_max,
        A_max,
    )
end

"""
    WahlSystematics(table, system) -> Union{WahlSystematics,Nothing}

The systematics for a fissioning system, taking the precursor excitation energy from the
compound nucleus: zero for spontaneous fission, `Sₙ + Eₙ` otherwise. Returns `nothing` when the
mass table cannot form `Sₙ`.
"""
function WahlSystematics(table::MassTable, system::FissioningSystem)
    excitation = compound_nucleus_excitation(table, system)
    excitation === nothing && return nothing
    return WahlSystematics(system.compound.Z, system.compound.A, value(excitation))
end

_wahl_peak_polarization(w::WahlSystematics, A) = w.ΔZ140 + w.ΔZSL * (A - 140)
_wahl_peak_width(w::WahlSystematics, A) = w.σZ140 + w.σZSL * (A - 140)

"""
    wahl_polarization(systematics, A) -> Float64

Charge polarization `ΔZ` of the heavy fragment at pre-neutron mass `A`, eqs. (11a), (12b),
(12d), (12f), (13b) and (14b).

Defined for `A ≥ A_F/2`. The light partner takes `−ΔZ` at its complement, which is how
[`most_probable_charge`](@ref) is applied here and what eq. (9d) states.
"""
function wahl_polarization(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return _wahl_peak_polarization(systematics, systematics.B5)
    A >= systematics.B5 && return _wahl_peak_polarization(systematics, systematics.B5) -
           systematics.ΔZSLW * (A - systematics.B5)
    A >= systematics.B4 && return _wahl_peak_polarization(systematics, A)
    edge =
        _wahl_peak_polarization(systematics, systematics.B4) +
        systematics.SL50 * (systematics.B4 - systematics.Bb)
    A >= systematics.Bb && return _wahl_peak_polarization(systematics, systematics.B4) +
           systematics.SL50 * (systematics.B4 - A)
    # Ba to Bb: the Z = 50 crossing, a straight line between values equal and opposite by
    # complementarity, so ΔZ vanishes at symmetry as charge conservation requires.
    return -edge + 2 * edge * (A - systematics.Ba) / (systematics.Bb - systematics.Ba)
end

"""
    wahl_width(systematics, A) -> Float64

Dispersion `σ_Z` of the isobaric charge distribution at pre-neutron mass `A`, eqs. (11b), (12c),
(12e), (12g), (13d) and (14c). Defined for `A ≥ A_F/2`.

Unlike `ΔZ` this is a step function: eqs. (12c) and (12g) hold `σ_Z(50)` flat over the two legs
adjoining the `Z = 50` crossing while the crossing itself keeps the peak value, and eq. (14c)
returns the far wing to `σ_Z(B5)` rather than continuing the wing slope. Fig. 18a, p. 31, draws
both discontinuities.
"""
function wahl_width(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return _wahl_peak_width(systematics, systematics.B5)
    A >= systematics.B5 && return _wahl_peak_width(systematics, systematics.B5) +
           systematics.σZSLW * (A - systematics.B5)
    A >= systematics.B4 && return _wahl_peak_width(systematics, A)
    A >= systematics.Bb && return systematics.σZ50
    return _wahl_peak_width(systematics, systematics.Bb)
end

function _wahl_check_heavy(systematics::WahlSystematics, A::Real)
    2 * A >= systematics.A_F ||
        throw(ArgumentError("the Zₚ systematics is evaluated for the heavy fragment, \
                       A ≥ A_F/2 = $(systematics.A_F / 2), got $A"))
    return nothing
end

"""
    wahl_charge_polarization(systematics, masses) -> ChargePolarization
    wahl_charge_polarization(table, system, masses) -> Union{ChargePolarization,Nothing}

Tabulate the `Zₚ` systematics over heavy-fragment masses `masses`, as the same
[`ChargePolarization`](@ref) a measured table would produce. This is the route for a fissioning
system with no evaluated `ΔZ(A)`, `rms(A)` of its own: it is a fitted systematics rather than a
measurement, but it carries the mass dependence a single mean value discards.

```julia
polarization = wahl_charge_polarization(masses, spontaneous_fission(Nuclide(98, 252)), 126:174)
```
"""
function wahl_charge_polarization(
    systematics::WahlSystematics,
    masses::AbstractVector{<:Integer},
)
    isempty(masses) && throw(ArgumentError("no masses to tabulate the systematics over"))
    polarization = Dict{Int, Float64}()
    width = Dict{Int, Float64}()
    for A in masses
        polarization[Int(A)] = wahl_polarization(systematics, A)
        width[Int(A)] = wahl_width(systematics, A)
    end
    return ChargePolarization(
        polarization,
        width,
        DEFAULT_POLARIZATION,
        DEFAULT_WIDTH,
        "Wahl Zₚ systematics, LA-13928 (2002), " *
        "Z_F = $(systematics.Z_F), A_F = $(systematics.A_F), " *
        "PE = $(round(systematics.PE; digits = 3)) MeV",
    )
end

wahl_charge_polarization(systematics::WahlSystematics, masses::AbstractRange{<:Integer}) =
    wahl_charge_polarization(systematics, collect(masses))

function wahl_charge_polarization(
    table::MassTable,
    system::FissioningSystem,
    masses::AbstractVector{<:Integer},
)
    systematics = WahlSystematics(table, system)
    systematics === nothing && return nothing
    return wahl_charge_polarization(systematics, masses)
end

wahl_charge_polarization(
    table::MassTable,
    system::FissioningSystem,
    masses::AbstractRange{<:Integer},
) = wahl_charge_polarization(table, system, collect(masses))

function Base.show(io::IO, systematics::WahlSystematics)
    print(
        io,
        "WahlSystematics(Z_F = ",
        systematics.Z_F,
        ", A_F = ",
        systematics.A_F,
        ", PE = ",
        round(systematics.PE; digits = 3),
        " MeV)",
    )
    return nothing
end
