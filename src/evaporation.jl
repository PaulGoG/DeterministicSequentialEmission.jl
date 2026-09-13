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

# Constants of the α term.
#
# Read the scaling before changing anything here. `α` multiplies `1/√ε`, so it carries units of
# √energy, but `π ħ² S₀ / (m r₀² A^(2/3))` as the papers print it evaluates to MeV: ħ²/m is
# 41.44 MeV fm², S₀ is dimensionless, r₀²A^(2/3) is fm². The missing half power is the
# strength function's own reference energy. S₀ is defined through the reduced neutron width
# Γ_n⁰ = Γ_n √(E₀/E) with E₀ = 1 eV, so the transmission coefficient is 2πS₀√(ε/E₀) rather than
# 2πS₀√ε, and the honest expression is
#
#     α = π ħ² S₀ / (m r₀² A^(2/3) √E₀),   1/√E₀ = 1/√(10⁻⁶ MeV) = 10³ MeV^(-1/2)
#
# The `0.1` on the radius and the `10` in `alpha_parameter` are how that 10³ is spelled here:
# (1/0.1²)·10 = 10³ = √(1 MeV / 1 eV) exactly, for every A. They are NOT a fm²-to-barn
# conversion — a genuine round trip would be a net factor of 1, and S₀ has no area units to
# convert. Closed form: α(A) = 9.041 (S₀/10⁻⁴) / A^(2/3) MeV^(1/2).
#
# Making the units "consistent" as a cleanup would divide α by a thousand and collapse the
# variable cross section onto the constant one, silently. The magnitude is checked by the
# paper's own figure: the deviation of 2T from ⟨ε⟩ must run 13-14 % at T ≈ 0.2 MeV and 7-8 %
# near T = 2 MeV, which α ~ 0.2-0.8 reproduces and α ~ 4×10⁻⁴ would miss by three orders.
const REDUCED_PLANCK_TIMES_C = 197.3268601   # MeV fm
const ATOMIC_MASS_UNIT = 931.50176           # MeV
const NEUTRON_MASS_NUMBER = 1.008665         # u
const RADIUS_PARAMETER = 1.2                 # fm
const ALPHA_COEFFICIENT =
    (π * REDUCED_PLANCK_TIMES_C)^2 / (NEUTRON_MASS_NUMBER * ATOMIC_MASS_UNIT)

"""
    neutron_strength_function(A) -> Float64

The s-wave neutron strength function `S₀(A)`, a jointed-segment parameterization of the RIPL-3
resonance data.

Returned bare, not in the customary units of `10⁻⁴`: `S₀(140)` is `2.33e-4`, not `2.33`. That
is what [`alpha_parameter`](@ref) expects.

Piecewise linear in `A` with plateaux where the data support one. Distinct from the
parameterization used in earlier work, and one of the two elements the thesis presents as new.

Verified against RIPL-3. The resonance table it was fitted to is byte-identical to
`resonances0.dat` as the IAEA distributes it, and over the 266 nuclides carrying a measured
`S₀` this reproduces them with a median ratio of 1.02 and 76 % of nuclides inside 2σ. It is a
sixfold improvement on the parameterization it replaced, which is essentially flat at `1e-4`
and biases `S₀` low by about 29 %: `χ²/N` of 5.4 against 31.3.

Two caveats, both real and neither a defect of the fit.

The residual scatter is Z-dependence that no function of `A` alone can carry. Worst in
`A = 111`–`124`, where the data span a factor of twelve at fixed mass — `A = 122` has Sn at
`0.13e-4` against Te at `1.22e-4` — and this follows the low tin branch, understating the
tellurium and cadmium fragments severalfold. The single worst segment is `A = 136`–`146`, the
rising edge to the 4s maximum, which runs straight through the `N = 82` closure.

And the published validation of the variable cross section does not cover the range this
implies. `⟨ε⟩(T)` was shown to be within 4 % of an optical-model calculation for `S₀` between
`0.7e-4` and `1.1e-4`, i.e. `α` between 0.21 and 0.55 — a band containing only 14 % of the
RIPL-3 nuclides in the fragment range. This parameterization runs `α` from 0.037 at `A = 120`
to 1.24 at `A = 147` because *the data run that far*, not because the fit strays. But the 4 %
figure does not transfer, and it is weakest where it matters most: the heavy fragment group of
²⁵²Cf peaks near `A ≈ 140`–`145`, on that worst segment, where `α` swings from 0.31 at
`A = 132` to 1.23 at `A = 148`.
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

Returns `nothing` when the equation admits no acceptable root — in particular below
`minimum_temperature`, where the level density has left the Fermi-gas regime for the
constant-temperature one and a Fermi-gas temperature has no meaning.
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

Lower bound on an acceptable residual temperature, in MeV.

A numerical floor, not a physical threshold, and the distinction is worth keeping straight. The
Fermi gas really does give way to the constant-temperature regime as a fragment cools, but that
crossover sits near `T ≈ 0.5`–`0.8` MeV for these nuclei, fifty times higher than this bound —
so this value cannot be justified as marking it, and a chain that stops here has already spent
several sequences outside the regime the closure assumes. What it does do is reject roots
returned where there is essentially nothing left to share: it bites only when `Ē_r − Sₙ` has
fallen below about 17 keV, so it truncates a chain that was ending anyway.

Unconditional, unlike the optional upper bound. Inherited from the prototype, where it is the
`ΔT` of the input options; no paper in the reference set states either bound.
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
