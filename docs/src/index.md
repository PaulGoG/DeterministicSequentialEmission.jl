```@meta
CurrentModule = DeterministicSequentialEmission
```

# DeterministicSequentialEmission

Prompt neutron emission from fission fragments, treated sequence by sequence.

The model solves the energy conservation equation of every emission sequence separately,

```
Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … n(A,Z,TKE)
```

with Weisskopf-Ewing evaporation spectra in the centre-of-mass system and fragment level
densities in the Fermi-gas regime. Because the number of sequences `n` is an integer at every
configuration, the post-neutron fragment mass `A_p = A − n` is an integer without rounding.

A run is driven by a TOML configuration and writes its tables; see [Running a case](running.md).
The vocabulary every name in this package is drawn from is set out in [Naming](naming.md).
The reference below follows the order a configuration is consumed in: the fissioning system and
its nuclides, the tabular input contract, the mass evaluation, the charge distribution and the
fragmentation domain it produces, the level density, the energy balance, the TXE partition, the
evaporation and the emission chain it drives, the yields, the contractions that turn a
configuration-resolved quantity into an observable, and the spectrum.

```@index
```

```@docs
DeterministicSequentialEmission
PACKAGE_ROOT
```

