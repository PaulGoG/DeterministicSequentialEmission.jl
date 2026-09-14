"""
    WAHL_LOW_ENERGY

Coefficients of the `Zₚ` model of A. C. Wahl, *Systematics of Fission-Product Yields*, LA-13928
(2002), Table 2, for the low-energy branch `PE ≤ 8 MeV`. Each entry is `(P₁ … P₅)` of eq. (17),

```
Par = P₁ + P₂[Z_F − 92] + P₃[A_F − 236] + P₄[PE − 6.551] + P₅[A_F − 236]²
```

The even-odd factors `F_Z` and `F_N` are carried but not applied by default; see
[`wahl_charge_distribution`](@ref) for what turning them on does and does not buy.

Table 2 also lists `FZ SL = 0.0030`, a slope of `F_Z` through the peak regions. Eq. (11d) as the
report prints it makes `F_Z(A')` constant there, so no slope is applied, and the parameter is not
carried. Footnote b gives the `F_N` slope as `−0.0006(10)`, consistent with zero.
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
    FZ140 = (1.207, 0.0, -0.0420, 0.0, 0.0022),
    FN140 = (1.076, 0.0, 0.0, 0.0, 0.0),
    FZSLW = (0.159, -0.028, 0.0, 0.0, 0.0),
    FNSLW = (0.039, 0.0, 0.0, 0.0, 0.0),
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
convention of [`ChargeDistribution`](@ref) and of eqs. (9c) and (9d) of the report.

The boundaries follow eq. (10). `Ba` and `Bb` are the ends of the steep interval where `Zₚ`
crosses the `Z = 50` shell closure; the report prints `Ba = A'_max` and `Bb = A_F − A'_max`,
which puts the heavy-side crossing at the lower end. Figs. 18b, 19b and 20b place the labels the
other way round, and only that reading is continuous with the peak regions on either side, so it
is the one used. `Bb` then lands within a mass unit of where the evaluated tables put the maximum
of `ΔZ`, which is what it is defined to be; the constructor records the measurement.

**A tabulated `ΔZ(A)`, `σ_Z(A)` evaluated for the reaction at hand is the better object and takes
precedence wherever one exists.** Fig. 18 quotes reduced `χ²` of 7.9 for these systematics against
2.9 for a per-reaction least-squares fit of the same model, and the evaluated tables are the
latter. This is the fallback for a system that has none, in place of a single mean `ΔZ` and `σ_Z`
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
    FZ140::Float64
    FN140::Float64
    FZSLW::Float64
    FNSLW::Float64
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
    FZ140 = evaluate(:FZ140)
    FN140 = evaluate(:FN140)
    FZSLW = evaluate(:FZSLW)
    FNSLW = evaluate(:FNSLW)

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
        FZ140,
        FN140,
        FZSLW,
        FNSLW,
        B2,
        B4,
        A_F - B2,
        A_F - 70.0,
        A_F - A_max,
        A_max,
    )
end

"""
    is_wahl_applicable(table, system) -> Bool

Whether the `Zₚ` systematics can be evaluated for a fissioning system at all: the compound
nucleus inside the fitted `Z_F` and `A_F` ranges of [`WAHL_VALIDITY`](@ref), and a precursor
excitation energy the mass table can form and the low-energy branch covers.

Asked before constructing, so that a caller offering the systematics as one option among several
can leave it out where it does not apply instead of catching the constructor's refusal.
"""
function is_wahl_applicable(table::MassExcessTable, system::FissioningSystem)
    system.compound.Z in WAHL_VALIDITY.charge || return false
    system.compound.A in WAHL_VALIDITY.mass || return false
    excitation = compound_nucleus_excitation(table, system)
    excitation === nothing && return false
    PE = value(excitation)
    return isfinite(PE) && 0 <= PE <= WAHL_VALIDITY.excitation
end

"""
    WahlSystematics(table, system) -> Union{WahlSystematics,Nothing}

The systematics for a fissioning system, taking the precursor excitation energy from the
compound nucleus: zero for spontaneous fission, `Sₙ + Eₙ` otherwise. Returns `nothing` when the
mass table cannot form `Sₙ`.
"""
function WahlSystematics(table::MassExcessTable, system::FissioningSystem)
    excitation = compound_nucleus_excitation(table, system)
    excitation === nothing && return nothing
    return WahlSystematics(system.compound.Z, system.compound.A, value(excitation))
end

_wahl_peak_polarization(w::WahlSystematics, A) = w.ΔZ140 + w.ΔZSL * (A - 140)
_wahl_peak_dispersion(w::WahlSystematics, A) = w.σZ140 + w.σZSL * (A - 140)

"""
    wahl_polarization(systematics, A) -> Float64

Charge distribution `ΔZ` of the heavy fragment at pre-neutron mass `A`, eqs. (11a), (12b),
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
    wahl_dispersion(systematics, A) -> Float64

Dispersion `σ_Z` of the isobaric charge distribution at pre-neutron mass `A`, eqs. (11b), (12c),
(12e), (12g), (13d) and (14c). Defined for `A ≥ A_F/2`.

Unlike `ΔZ` this is a step function: eqs. (12c) and (12g) hold `σ_Z(50)` flat over the two legs
adjoining the `Z = 50` crossing while the crossing itself keeps the peak value, and eq. (14c)
returns the far wing to `σ_Z(B5)` rather than continuing the wing slope. Fig. 18a, p. 31, draws
both discontinuities.
"""
function wahl_dispersion(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return _wahl_peak_dispersion(systematics, systematics.B5)
    A >= systematics.B5 && return _wahl_peak_dispersion(systematics, systematics.B5) +
           systematics.σZSLW * (A - systematics.B5)
    A >= systematics.B4 && return _wahl_peak_dispersion(systematics, A)
    A >= systematics.Bb && return systematics.σZ50
    return _wahl_peak_dispersion(systematics, systematics.Bb)
end

"""
    wahl_even_odd_factors(systematics, A) -> (F_Z, F_N)

The even-odd proton and neutron factors at heavy-fragment mass `A`, eqs. (11d), (11e), (12a),
(13f), (13h), (14d) and (14e).

Both are `1` through the near-symmetry region, where eq. (12a) sets `F(A) = 1`, and both take
their peak values through the peak and far-wing regions. Only the wing regions carry a slope.

These are the factors the plain-Gaussian charge distribution of this package does **not** apply.
[`wahl_charge_distribution`](@ref) can fold them in; what that is worth is set out there.
"""
function wahl_even_odd_factors(systematics::WahlSystematics, A::Real)
    _wahl_check_heavy(systematics, A)
    A >= systematics.B6 && return (systematics.FZ140, systematics.FN140)
    A >= systematics.B5 && return (
        systematics.FZ140 + systematics.FZSLW * (A - systematics.B5),
        systematics.FN140 + systematics.FNSLW * (A - systematics.B5),
    )
    A >= systematics.B4 && return (systematics.FZ140, systematics.FN140)
    return (1.0, 1.0)
end

"""
    wahl_charge_yield(systematics, A, charges) -> (Vector{Int}, Vector{Float64})

The fractional independent yields `FI(Z, A)` of the `charges` charge numbers nearest `Zₚ(A)`,
Gaussian modulated by the even-odd factors and renormalized to sum to one.

The renormalization is not optional and is not a truncation correction. LA-13928 p. 21: the factor
"applied to achieve Σ(FI) = 1.00 for each A, is required because the even-odd factors, F(A),
destroy the inherent normalization properties of Gaussian distributions". A plain Gaussian keeps
its own normalization and needs none, which is why nothing else in this package renormalizes.

`F(A)` is composed from the parities of `Z` and `N = A − Z`: `F_Z F_N` for even-even,
`F_Z/F_N` for even-odd, `F_N/F_Z` for odd-even and `1/(F_Z F_N)` for odd-odd.
"""
function wahl_charge_yield(systematics::WahlSystematics, A::Integer, charges::Integer = 11)
    charges > 0 || throw(ArgumentError("need at least one charge number, got $charges"))
    F_Z, F_N = wahl_even_odd_factors(systematics, A)
    ΔZ = wahl_polarization(systematics, A)
    σ_Z = wahl_dispersion(systematics, A)
    Zₚ = A * systematics.Z_F / systematics.A_F + ΔZ
    half = charges ÷ 2
    centre = round(Int, Zₚ)
    numbers = (centre - half):(centre + half)
    yields = Float64[]
    for Z in numbers
        N = A - Z
        factor = if iseven(Z) && iseven(N)
            F_Z * F_N
        elseif iseven(Z)
            F_Z / F_N
        elseif iseven(N)
            F_N / F_Z
        else
            1 / (F_Z * F_N)
        end
        push!(yields, charge_probability(Z, Zₚ, σ_Z) * factor)
    end
    total = sum(yields)
    total > 0 || throw(ArgumentError("the charge yields at A = $A sum to $total"))
    yields ./= total
    return collect(numbers), yields
end

function _wahl_check_heavy(systematics::WahlSystematics, A::Real)
    2 * A >= systematics.A_F ||
        throw(ArgumentError("the Zₚ systematics is evaluated for the heavy fragment, \
                       A ≥ A_F/2 = $(systematics.A_F / 2), got $A"))
    return nothing
end

"""
    wahl_charge_distribution(systematics, masses) -> ChargeDistribution
    wahl_charge_distribution(table, system, masses) -> Union{ChargeDistribution,Nothing}

Tabulate the `Zₚ` systematics over heavy-fragment masses `masses`, as the same
[`ChargeDistribution`](@ref) a measured table would produce. This is the route for a fissioning
system with no evaluated `ΔZ(A)`, `σ_Z(A)` of its own: it is a fitted systematics rather than a
measurement, but it carries the mass dependence a single mean value discards.

```julia
distribution = wahl_charge_distribution(masses, spontaneous_fission(Nuclide(98, 252)), 126:174)
```

# The even-odd option

`even_odd = true` folds in the factors `F_Z` and `F_N` the model carries and this package
otherwise omits: `p(Z,A)` is built by [`wahl_charge_yield`](@ref), which is a Gaussian modulated
by `F(A)` and renormalized, and a plain Gaussian is then fitted to it at each mass by its first
two moments. `ΔZ(A)` is that Gaussian's centre less `Z_UCD(A)`, and `σ_Z(A)` its width. `charges`
is how many charge numbers the fit runs over.

**This is the construction of the evaluated tables**, which are effective Gaussian parameters of
exactly this kind rather than the model's own `σ_Z`. It is why those tables oscillate with a
period of about five masses — `Zₚ` advances by `Z_F/A_F` per mass unit, so even-charge dominance
recurs every `2A_F/Z_F ≈ 5.1` — and a smooth systematics cannot.

**What it buys, measured against the three evaluated tables.** It reproduces the ripple: the
modulation it introduces is `0.36` in `ΔZ` and `0.30` in `σ_Z`, against `0.32` and `0.30` actually
present, and the two correlate at `0.62` for ²³⁵U, `0.36` for ²³⁹Pu and `0.07` for ²⁵²Cf —
following the even-odd effect's own decline with the mass of the fissioning system.

**What it does not buy.** The mean absolute difference from an evaluated table barely moves,
because that difference is dominated by the systematics' own parameter error and not by the
missing even-odd layer. Fig. 18 of LA-13928 quantifies it: reduced `χ²` of 7.9 for the systematics
against 2.9 for a per-reaction least-squares fit. Eq. (17) predicts the parameters *on average
across reactions*; an evaluated table was built from that reaction's own fitted parameters. So
this is the better of two estimates where no table exists, and not a route to reproducing one.

Off by default, because the tabulated distributions this package normally reads already carry the
even-odd effect inside their fitted widths, and folding the factors in on top of one of those
would count it twice.
"""
function wahl_charge_distribution(
    systematics::WahlSystematics,
    masses::AbstractVector{<:Integer};
    even_odd::Bool = false,
    charges::Integer = 11,
)
    isempty(masses) && throw(ArgumentError("no masses to tabulate the systematics over"))
    distribution = Dict{Int, Float64}()
    σ_Z = Dict{Int, Float64}()
    for A in masses
        if even_odd
            distribution[Int(A)], σ_Z[Int(A)] =
                _wahl_effective_gaussian(systematics, Int(A), charges)
        else
            distribution[Int(A)] = wahl_polarization(systematics, A)
            σ_Z[Int(A)] = wahl_dispersion(systematics, A)
        end
    end
    return ChargeDistribution(
        distribution,
        σ_Z,
        DEFAULT_CHARGE_POLARIZATION,
        DEFAULT_CHARGE_DISPERSION,
        "Wahl Zₚ systematics, LA-13928 (2002), " *
        (even_odd ? "even-odd factors folded in, " : "") *
        "Z_F = $(systematics.Z_F), A_F = $(systematics.A_F), " *
        "PE = $(round(systematics.PE; digits = 3)) MeV",
    )
end

# The plain Gaussian that best represents the modulated distribution, by its first two moments.
# Returns (ΔZ, σ_Z), the centre expressed against unchanged charge division.
function _wahl_effective_gaussian(
    systematics::WahlSystematics,
    A::Integer,
    charges::Integer,
)
    numbers, yields = wahl_charge_yield(systematics, A, charges)
    centre = sum(Z * y for (Z, y) in zip(numbers, yields))
    variance = sum((Z - centre)^2 * y for (Z, y) in zip(numbers, yields))
    return centre - A * systematics.Z_F / systematics.A_F, sqrt(variance)
end

wahl_charge_distribution(
    systematics::WahlSystematics,
    masses::AbstractRange{<:Integer};
    kwargs...,
) = wahl_charge_distribution(systematics, collect(masses); kwargs...)

function wahl_charge_distribution(
    table::MassExcessTable,
    system::FissioningSystem,
    masses::AbstractVector{<:Integer};
    kwargs...,
)
    systematics = WahlSystematics(table, system)
    systematics === nothing && return nothing
    return wahl_charge_distribution(systematics, masses; kwargs...)
end

wahl_charge_distribution(
    table::MassExcessTable,
    system::FissioningSystem,
    masses::AbstractRange{<:Integer};
    kwargs...,
) = wahl_charge_distribution(table, system, collect(masses); kwargs...)

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
