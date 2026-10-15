```@meta
CurrentModule = DeterministicSequentialEmission
```

# DeterministicSequentialEmission

Prompt neutron emission from fission fragments, treated sequence by sequence. The literature behind
the formulation and the input data is listed with DOIs under [References](references.md).

The model solves the energy conservation equation of every emission sequence separately,

```
Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … k_max
```

with Weisskopf-Ewing evaporation spectra in the centre-of-mass system and fragment level
densities in the Fermi-gas regime. `Ē_r^(k) = a_k T_k²` is the residual excitation energy after
the `k`-th emission, with `Ē_r^(0) = E*`, the excitation energy of the fragment; `Sₙ^(k-1)` is
the neutron separation energy of the nucleus that emits the `k`-th neutron; `⟨ε⟩_k` is the mean
centre-of-mass energy of that neutron; `a_k` and `T_k` are the level density parameter and the
temperature of the `k`-th residual nucleus. The equation is solved for `k = 1 … k_max`, while the
residual nucleus can afford `Sₙ` and the equation admits a temperature. The emissions that remain
possible once the neutrons' kinetic energy is returned to the residual are then counted on `Sₙ`
alone, which gives the number of sequences `n ≥ k_max` of the configuration `(A, Z, TKE)` and
its multiplicity `ν = (n+1)/2`.

Because `n` is an integer at every configuration, every product mass `A − 1 … A − n` is an
integer, so the products are distributed without rounding. By default a chain ends on each of
its `n` residuals with probability `1/n`, so the mass its products lose is its multiplicity;
`[products] count = "sequence_count"` takes the published `A_p = A − n` [Tudora2022](@cite)
instead.

The physics is the literature's; the implementation is mine. The model is that of Tudora, Hambsch
and Tobosaru [TudoraHambschTobosaru2018](@cite), with the systematics of its residual quantities
[Tudora2020](@cite), its validation against measurement [Tudora2019](@cite), its post-neutron
fragment yields [Tudora2022](@cite) and its place beside the point-by-point model
[Tudora2021](@cite). The partition of the total excitation energy by the temperature ratio and
the evaporation framework are those of the point-by-point model [TudoraHambsch2017](@cite); the
extraction of the temperature ratio from measured multiplicities, which
[FissionTemperatureRatio.jl](https://github.com/PaulGoG/FissionTemperatureRatio.jl) implements,
is [TudoraGogita2024](@cite). The evaporation spectrum in closed form and its
transformation to the laboratory follow Madland and Nix [MadlandNix1982](@cite) and Terrell
[Terrell1959](@cite); the conventions of the evaluated spectra are those of Capote et al.
[Capote2016](@cite). This package is my Julia implementation of that physics.

A run is driven by a TOML configuration and writes its tables; see [Running a case](running.md),
which also documents the configuration and the pipeline. The vocabulary every name in this
package is drawn from is set out in [Naming](naming.md). What the model takes from the fragments
before emission begins — nuclides and systems, masses, the charge distribution and the
fragmentation domain, level densities, the energy balance, yields and the temperature-ratio run
record — lives in the domain package, summarised under
[The fragmentation domain](domain.md). The reference of this package follows the order a
configuration is consumed in: [the excitation-energy partition](energetics.md), then
[the evaporation, the emission chain and the spectrum](emission.md).

```@index
```

```@docs
DeterministicSequentialEmission
PACKAGE_ROOT
```

