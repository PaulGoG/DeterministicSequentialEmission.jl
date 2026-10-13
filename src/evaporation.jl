"""
    EvaporationCrossSection

Analytical form of the cross section for compound nucleus formation in the process inverse to
evaporation, `σ_c(ε)`.

The choice determines how the sequential equations are solved: a constant cross section
reduces them to a quadratic with a closed-form root, while an energy-dependent one makes them
transcendental.
"""
abstract type EvaporationCrossSection end

"""
    ConstantCrossSection()

Energy-independent `σ_c(ε) = σ₀`.

The evaporation spectrum and its first moment become

```
φ(ε; T) = (ε/T²) exp(−ε/T)
⟨ε⟩(T)  = 2T
```

so the sequential equation `a T² + 2T − (Ē_r − Sₙ) = 0` is quadratic with the positive root

```
T = [√(1 + a(Ē_r − Sₙ)) − 1] / a
```
"""
struct ConstantCrossSection <: EvaporationCrossSection end

"""
    VariableCrossSection()

Iwamoto's energy-dependent form, a constant term plus one carrying the s-wave neutron strength
function,

```
σ_c(ε) = σ₀ (1 + α/√ε),   α = π ħ² S₀ / (m r₀² A^(2/3))
```

which gives

```
φ(ε; T) = (ε + α√ε) exp(−ε/T) / [T^(3/2)(√T + α√π/2)]
⟨ε⟩(T)  = T (2√T + 3α√π/4) / (√T + α√π/2)
```

The sequential equation is then transcendental and is solved numerically. It has exactly one
positive root: `f(0) = −(Ē_r − Sₙ) < 0`, `f → +∞`, and `f′ > 0` throughout, with
`f′(0) = 3/2`.
"""
struct VariableCrossSection <: EvaporationCrossSection end

# Constants of the α term. `α` carries units of MeV^(1/2): S₀ is defined through the reduced
# width Γ_n⁰ = Γ_n √(E₀/E), E₀ = 1 eV, so
#
#     α = π ħ² S₀ / (m r₀² A^(2/3) √E₀),   1/√E₀ = 10³ MeV^(-1/2)
#
# The `0.1` on the radius and the `10` in `alpha_parameter` give (1/0.1²)·10 = 10³ = 1/√E₀;
# they are not an area-unit conversion. Closed form: α(A) = 9.041 (S₀/10⁻⁴) / A^(2/3) MeV^(1/2).
const REDUCED_PLANCK_TIMES_C = 197.3268601   # MeV fm
const ATOMIC_MASS_UNIT = 931.50176           # MeV
const NEUTRON_MASS_NUMBER = 1.008665         # u
const RADIUS_PARAMETER = 1.2                 # fm
const ALPHA_COEFFICIENT =
    (π * REDUCED_PLANCK_TIMES_C)^2 / (NEUTRON_MASS_NUMBER * ATOMIC_MASS_UNIT)

"""
    neutron_strength_function(A) -> Float64

The s-wave neutron strength function `S₀(A)`: a piecewise-linear interpolation, with plateaux,
through selected nuclei of the RIPL-3 s-wave resonance table, Capote et al., Nucl. Data Sheets
110, 3107 (2009), doi:10.1016/j.nds.2009.10.004 — ²³Na, ⁵⁵Mn, ¹²⁰Sn, ¹⁴⁶Nd, ¹⁵²Gd among the
nodes — not a least-squares fit. Returned bare, not in the customary units of `10⁻⁴`: `S₀(140)`
is `2.33e-4`, as [`alpha_parameter`](@ref) expects. In cross-validation against the table it
describes the stable nuclei as well as any function of `A` alone; as such a function it cannot
carry the `Z`-dependence of `S₀`, largest at `A = 111`–`124` and across the `N = 82` closure at
`A = 136`–`146`, and nothing is measured for the neutron-rich fragments 4–8 neutrons beyond the
table. The resulting `α` runs from 0.037 at `A = 120` to 1.24 at `A = 146`. The closed-form
`⟨ε⟩(T)` of the variable cross section was compared with an optical-model (Becchetti–Greenlees)
calculation by Tudora, Hambsch and Tobosaru, Eur. Phys. J. A 54, 87 (2018),
doi:10.1140/epja/i2018-12521-7, §2.1, over `A = 76`–`160` and `S₀ = 0.7`–`1.1 × 10⁻⁴`, i.e.
`α = 0.21`–`0.55`, agreeing within 4 % for `T = 0.2`–`2` MeV; about half of the fragment masses
fall outside that envelope here. The observables are weakly sensitive to the choice: across the
alternatives tried, `ν(A)` is unchanged, `⟨ε⟩(A)` moves by at most 0.13 MeV at one mass,
`⟨E⟩` by under 1.5 % and the ²⁵²Cf spectrum ratio at 10 MeV by under 2 %.
"""
function neutron_strength_function(A::Integer)
    A <= 23 && return 2.5e-5
    A <= 55 && return 1.29695e-5 * A - 2.7328e-4
    A <= 69 && return -2.1429e-5 * A + 1.6186e-3
    A <= 79 && return 1.4e-4
    A <= 90 && return -7.8182e-6 * A + 7.5764e-4
    A <= 110 && return 5.4e-5
    A <= 120 && return -4.4e-6 * A + 5.38e-4
    A <= 135 && return 6.6667e-6 * A - 7.9e-4
    A <= 146 && return 2.4545e-5 * A - 3.2036e-3
    A <= 152 && return 3.8e-4
    A <= 158 && return -3.8333e-5 * A + 6.2067e-3
    A <= 173 && return 1.5e-4
    A <= 186 && return 1.0769e-5 * A - 1.7131e-3
    A <= 208 && return -8.6364e-6 * A + 1.8964e-3
    return 1e-4
end

"""
    alpha_parameter(A) -> Float64

The coefficient `α` of the energy-dependent term of `σ_c`, for a residual nucleus of mass `A`.
"""
function alpha_parameter(A::Integer)
    σ₀ = π * (RADIUS_PARAMETER * 0.1)^2 * A^(2 / 3)
    return 10 * ALPHA_COEFFICIENT * neutron_strength_function(A) / σ₀
end

"""
    mean_neutron_energy(::ConstantCrossSection, T) -> Float64
    mean_neutron_energy(::VariableCrossSection, α, T) -> Float64

First moment `⟨ε⟩` of the evaporation spectrum in the centre-of-mass system, in MeV.
"""
mean_neutron_energy(::ConstantCrossSection, T::Real) = 2 * T

function mean_neutron_energy(::VariableCrossSection, α::Real, T::Real)
    return T * (2 * sqrt(T) + 3 * α * sqrt(π) / 4) / (sqrt(T) + α * sqrt(π) / 2)
end

"""
    fermi_gas_energy(a, T) -> Float64

Residual excitation energy in the Fermi-gas regime, `Ē_r = a T²`, in MeV.
"""
fermi_gas_energy(a::Real, T::Real) = a * T^2

"""
    residual_temperature(cross_section, a, available; kwargs...) -> Union{Float64,Nothing}

Solve one sequential equation for the residual temperature `T`, where
`available = Ē_r^(k−1) − Sₙ^(k−1)` is the energy left after removing the neutron.

Returns `nothing` when the equation admits no acceptable root, in particular below
`minimum_temperature`, the numerical floor of [`MINIMUM_TEMPERATURE`](@ref).
"""
function residual_temperature(
    ::ConstantCrossSection,
    a::Real,
    available::Real;
    minimum_temperature::Real = MINIMUM_TEMPERATURE,
    maximum_temperature::Union{Real, Nothing} = nothing,
)
    available >= 0 || return nothing
    T = (sqrt(1 + a * available) - 1) / a
    return _accept_temperature(T, minimum_temperature, maximum_temperature)
end

function residual_temperature(
    cross_section::VariableCrossSection,
    a::Real,
    available::Real;
    α::Real,
    minimum_temperature::Real = MINIMUM_TEMPERATURE,
    maximum_temperature::Union{Real, Nothing} = nothing,
)
    available >= 0 || return nothing
    equation(T) = a * T^2 + mean_neutron_energy(cross_section, α, T) - available
    T = _bisect(equation, 0.0, _upper_bracket(equation, a, available))
    T === nothing && return nothing
    return _accept_temperature(T, minimum_temperature, maximum_temperature)
end

"""
    MINIMUM_TEMPERATURE

Lower bound on an acceptable residual temperature, in MeV; the convention of the thesis.

A numerical floor, not the crossover from the Fermi-gas to the constant-temperature regime,
which lies near `T ≈ 0.5`–`0.8` MeV for fission fragments. It rejects roots only where
`Ē_r − Sₙ` has fallen below about 17 keV, truncating a chain that is already ending. Applied
unconditionally, unlike the optional upper bound.
"""
const MINIMUM_TEMPERATURE = 1.0e-2

function _accept_temperature(T, minimum_temperature, maximum_temperature)
    (isfinite(T) && T >= minimum_temperature) || return nothing
    maximum_temperature === nothing && return T
    return T <= maximum_temperature ? T : nothing
end

# Expand upwards until the equation changes sign, so the bracket is found rather than assumed.
function _upper_bracket(equation, a, available)
    upper = max(sqrt(max(available, 0.0) / max(a, eps())), 1.0)
    for _ in 1:64
        equation(upper) > 0 && return upper
        upper *= 2
    end
    return upper
end

function _bisect(equation, lower, upper; tolerance = 1.0e-12, iterations = 200)
    low = equation(lower)
    high = equation(upper)
    (isfinite(low) && isfinite(high)) || return nothing
    sign(low) == sign(high) && return nothing
    for _ in 1:iterations
        middle = (lower + upper) / 2
        value = equation(middle)
        (upper - lower) <= tolerance && return middle
        if sign(value) == sign(low)
            lower = middle
            low = value
        else
            upper = middle
        end
    end
    return (lower + upper) / 2
end
