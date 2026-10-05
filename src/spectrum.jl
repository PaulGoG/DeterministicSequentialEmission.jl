"""
    centre_of_mass_spectrum(cross_section, α, T, ε) -> Float64

The Weisskopf-Ewing evaporation spectrum of one emission sequence, in the centre-of-mass frame
of the emitting fragment, normalized to unit area.

```
φ(ε; T) = K(T) ε σ_c(ε) exp(−ε/T),   K(T) = [∫₀^∞ ε σ_c(ε) exp(−ε/T) dε]⁻¹
```

`σ₀` cancels in `K`, so only the shape of `σ_c` survives. With a constant cross section this is
`ε exp(−ε/T)/T²`; with Iwamoto's `σ₀(1 + α/√ε)` it is
`(ε + α√ε) exp(−ε/T) / [T^(3/2)(√T + α√π/2)]`. Both integrate to 1, and their first moments are
`2T` and [`mean_neutron_energy`](@ref) respectively.

`α` is ignored by the constant cross section and is the residual nucleus's
[`alpha_parameter`](@ref) for the variable one, the value the sequential equation was solved
with, so the spectrum is consistent with the energy balance that produced `T`.
"""
function centre_of_mass_spectrum(::ConstantCrossSection, ::Real, T::Real, ε::Real)
    (ε <= 0 || T <= 0) && return 0.0
    return ε * exp(-ε / T) / T^2
end

function centre_of_mass_spectrum(::VariableCrossSection, α::Real, T::Real, ε::Real)
    (ε <= 0 || T <= 0) && return 0.0
    return (ε + α * sqrt(ε)) * exp(-ε / T) / (T^(3 / 2) * (sqrt(T) + α * sqrt(π) / 2))
end

"""
    laboratory_spectrum(cross_section, α, T, E_f, E) -> Float64

One emission sequence's spectrum transformed into the laboratory frame.

A neutron leaving a fragment that moves with kinetic energy per nucleon `E_f` arrives with
`E = ε + E_f + 2√(ε E_f) cos θ`. Isotropy in the centre-of-mass makes `cos θ` uniform, so at
fixed `ε` the lab energy is uniform over a band of width `4√(ε E_f)`, and folding over `φ`:

```
N(E) = 1/(4√E_f) ∫_{u₁}^{u₂} φ(ε) ε^(-1/2) dε,   u_{1,2} = (√E ∓ √E_f)²
```

The factor `4√(ε E_f)` is the width of that band, by the identity
`(√ε+√E_f)² − (√ε−√E_f)² = 4√(ε E_f)`, so the kernel is measure-preserving and `N` inherits
unit area from `φ` exactly. With `ε = s²` the kernel is uniform in `s` over an interval of fixed
length `2√E_f`,

```
N(E) = 1/(2√E_f) ∫_{|√E−√E_f|}^{√E+√E_f} φ(s²) ds
```

Evaluated in closed form rather than by quadrature. With `a_i = u_i/T`,

```
constant σ_c:   N = [Γ(3/2,a₁) − Γ(3/2,a₂)] / (4√(E_f T))
variable σ_c:   N = { [Γ(3/2,a₁) − Γ(3/2,a₂)] + (α/√T)[exp(−a₁) − exp(−a₂)] }
                    / [4√E_f (√T + α√π/2)]
```

The second reduces to the first as `α → 0`. The differences are taken over the upper incomplete
gamma, since the equivalent lower-gamma difference cancels catastrophically in the tail and is
in error by a factor of four beyond about 20 MeV.

Source: eq. (12) of Madland & Nix, *Nucl. Sci. Eng.* **81**, 213 (1982), doi:10.13182/NSE82-5.
"""
function laboratory_spectrum(
    cross_section::EvaporationCrossSection,
    α::Real,
    T::Real,
    E_f::Real,
    E::Real,
)
    (E <= 0 || T <= 0 || E_f <= 0) && return 0.0
    root = sqrt(E)
    speed = sqrt(E_f)
    u₁ = (root - speed)^2
    u₂ = (root + speed)^2
    return _laboratory(cross_section, α, T, E_f, u₁, u₂)
end

# Γ(3/2, x), the upper incomplete gamma, normalized out of `gamma_inc`'s regularized pair.
_upper_gamma(x::Real) = gamma_inc(1.5, x)[2] * sqrt(π) / 2

function _laboratory(::ConstantCrossSection, ::Real, T::Real, E_f::Real, u₁::Real, u₂::Real)
    band = _upper_gamma(u₁ / T) - _upper_gamma(u₂ / T)
    return band / (4 * sqrt(E_f * T))
end

function _laboratory(
    ::VariableCrossSection,
    α::Real,
    T::Real,
    E_f::Real,
    u₁::Real,
    u₂::Real,
)
    band = _upper_gamma(u₁ / T) - _upper_gamma(u₂ / T)
    tail = α / sqrt(T) * (exp(-u₁ / T) - exp(-u₂ / T))
    return (band + tail) / (4 * sqrt(E_f) * (sqrt(T) + α * sqrt(π) / 2))
end

"""
    chain_spectrum(chain, cross_section, E_f, E[, weighting]) -> Float64

The laboratory spectrum of one fragment, averaged over the emission sequences that carry a
temperature.

```
N_f(E) = Σ_{k solved} w_k N(E; α_k, T_k, E_f) / Σ_{k solved} w_k
```

`w_k` is the [`sequence_weight`](@ref) of `weighting`, uniform (`w_k = 1`, denominator
`k_max`) where none is given. The uniform case is the average of
[`sequence_average`](@ref) and of eq. (12) of Tudora, *EPJA* 54, 87 (2018),
doi:10.1140/epja/i2018-12521-7. Sequences without a temperature carry no spectrum and are
skipped, although they count towards `n` and hence towards the multiplicity weighting of the
pair spectrum.

Returns `0` for a chain that solved nothing.
"""
function chain_spectrum(
    chain::EmissionChain,
    cross_section::EvaporationCrossSection,
    E_f::Real,
    E::Real,
    weighting::SequenceWeighting,
)
    total = 0.0
    norm = 0.0
    for sequence in chain.sequences
        is_solved(sequence) || continue
        weight = sequence_weight(weighting, chain, sequence)
        α = alpha_parameter(chain.fragment.A - sequence.index)
        total +=
            weight * laboratory_spectrum(cross_section, α, sequence.temperature, E_f, E)
        norm += weight
    end
    return norm > 0 ? total / norm : 0.0
end

chain_spectrum(
    chain::EmissionChain,
    cross_section::EvaporationCrossSection,
    E_f::Real,
    E::Real,
) = chain_spectrum(chain, cross_section, E_f, E, UniformSequenceWeighting())

"""
    chain_centre_of_mass_spectrum(chain, cross_section, ε[, weighting]) -> Float64

The sequence-averaged centre-of-mass spectrum of one fragment,
`Φ(ε) = Σ_{k solved} w_k φ_k(ε) / Σ_{k solved} w_k`, with `w_k` the
[`sequence_weight`](@ref) of `weighting`, uniform (`w_k = 1`, denominator `k_max`) where none is
given.

Same weights and the same skip rule as [`chain_spectrum`](@ref). This is the quantity
measured centre-of-mass spectra are compared against, and its first moment is
`sequence_average(chain, :mean_neutron_energy, weighting)` by construction.
"""
function chain_centre_of_mass_spectrum(
    chain::EmissionChain,
    cross_section::EvaporationCrossSection,
    ε::Real,
    weighting::SequenceWeighting,
)
    total = 0.0
    norm = 0.0
    for sequence in chain.sequences
        is_solved(sequence) || continue
        weight = sequence_weight(weighting, chain, sequence)
        α = alpha_parameter(chain.fragment.A - sequence.index)
        total += weight * centre_of_mass_spectrum(cross_section, α, sequence.temperature, ε)
        norm += weight
    end
    return norm > 0 ? total / norm : 0.0
end

chain_centre_of_mass_spectrum(
    chain::EmissionChain,
    cross_section::EvaporationCrossSection,
    ε::Real,
) = chain_centre_of_mass_spectrum(chain, cross_section, ε, UniformSequenceWeighting())

"""
    SpectrumWeighting

How the light and heavy fragment spectra are combined into the pair spectrum.

`MultiplicityWeighting` uses `ν_{L,H}/ν_pair`, eq. (16) of Tudora, *EPJA* 54, 87 (2018),
doi:10.1140/epja/i2018-12521-7, and eq. (33) of Capote *et al.*, *Nucl. Data Sheets* **131**, 1
(2016), doi:10.1016/j.nds.2015.12.002; it is the convention of the thesis and the one the
pipeline applies. `ExcitationWeighting` uses `E*_{L,H}/TXE`, eq. (31) of Tudora, EPJ Web Conf.
256, 00017 (2021), doi:10.1051/epjconf/202125600017. The two are close by construction, since
`E*_H/TXE` behaves like `ν_H/(ν_L + ν_H)`, but not identical. Only `MultiplicityWeighting`
combined with [`PerNeutronAverage`](@ref) makes the weight of a fragment spectrum in the total
its neutron count, `Y ν_{L,H}`; `ExcitationWeighting` with the same average is a hybrid.
"""
abstract type SpectrumWeighting end

"""
    MultiplicityWeighting()

Combine by `ν_{L,H}/ν_pair`.
"""
struct MultiplicityWeighting <: SpectrumWeighting end

"""
    ExcitationWeighting()

Combine by `E*_{L,H}/TXE`.
"""
struct ExcitationWeighting <: SpectrumWeighting end

"""
    pair_weights(weighting, light, heavy, energies) -> Union{Tuple{Float64,Float64},Nothing}

The light and heavy weights of one fragmentation, summing to 1, or `nothing` where the
fragmentation emits nothing and has no pair spectrum to weight.
"""
function pair_weights(
    ::MultiplicityWeighting,
    light::EmissionChain,
    heavy::EmissionChain,
    ::ExcitationEnergies,
)
    ν_light = neutron_multiplicity(light)
    ν_heavy = neutron_multiplicity(heavy)
    total = ν_light + ν_heavy
    total > 0 || return nothing
    return (ν_light / total, ν_heavy / total)
end

function pair_weights(
    ::ExcitationWeighting,
    ::EmissionChain,
    ::EmissionChain,
    energies::ExcitationEnergies,
)
    total = sum(energies)
    total > 0 || return nothing
    return (energies.light / total, energies.heavy / total)
end

"""
    SpectrumAverage

How the pair spectra of the fragmentations combine into the total laboratory spectrum.

[`PerNeutronAverage`](@ref), the default, weights each unit-normalised pair spectrum by the
yield times the pair multiplicity, `N = Σ Y ν_pair N_pair / Σ Y ν_pair`. Every emitted neutron
counts once, so this is the unit-area shape `χ/ν̄` of the neutron yield per fission that
evaluated spectra represent: eq. (99) of Capote *et al.*, *Nucl. Data Sheets* **131**, 1 (2016),
doi:10.1016/j.nds.2015.12.002. It is also the weight of eq. (53) of Madland and Nix, *Nucl. Sci.
Eng.* **81**, 213 (1982), doi:10.13182/NSE82-5, across fission chances, and of eqs. (8)–(9) of
Kawano *et al.*, *Phys. Rev. C* **104**, 014611 (2021), doi:10.1103/PhysRevC.104.014611, across
fragment pairs. [`PerFissionAverage`](@ref) weights it by the yield alone, `N = Σ Y N_pair / Σ Y`,
the form of the thesis and of Tudora, Hambsch, Tobosaru, *EPJA* **54**, 87 (2018),
doi:10.1140/epja/i2018-12521-7, text after eq. (16). The two differ where `ν_pair` correlates
with the hardness of `N_pair`, and the per-fission form under-weights the hot, high-multiplicity
fragmentations.
"""
abstract type SpectrumAverage end

"""
    PerFissionAverage()

Combine pair spectra by `Y`.
"""
struct PerFissionAverage <: SpectrumAverage end

"""
    PerNeutronAverage()

Combine pair spectra by `Y ν_pair`.
"""
struct PerNeutronAverage <: SpectrumAverage end

"""
    total_weight(average, yield, light, heavy) -> Float64

The weight of one fragmentation's pair spectrum in the total: `yield` under
[`PerFissionAverage`](@ref), `yield` times the pair multiplicity under
[`PerNeutronAverage`](@ref).
"""
total_weight(::PerFissionAverage, yield::Float64, ::EmissionChain, ::EmissionChain) = yield

function total_weight(
    ::PerNeutronAverage,
    yield::Float64,
    light::EmissionChain,
    heavy::EmissionChain,
)
    return yield * (neutron_multiplicity(light) + neutron_multiplicity(heavy))
end
