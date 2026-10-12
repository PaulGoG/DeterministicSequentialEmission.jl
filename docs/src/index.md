```@meta
CurrentModule = DeterministicSequentialEmission
```

# DeterministicSequentialEmission

Prompt neutron emission from fission fragments, treated sequence by sequence. The literature behind
the formulation and the input data is listed with DOIs under [References](references.md).

The model solves the energy conservation equation of every emission sequence separately,

```
Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … n(A,Z,TKE)
```

with Weisskopf-Ewing evaporation spectra in the centre-of-mass system and fragment level
densities in the Fermi-gas regime. Because the number of sequences `n` is an integer at every
configuration, the post-neutron fragment mass `A_p = A − n` is an integer without rounding.

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

