```@meta
CurrentModule = DeterministicSequentialEmission
```

# DeterministicSequentialEmission

DeterministicSequentialEmission.jl calculates prompt neutron emission from fission fragments in low-energy
and thermal fission, tracking de-excitation cascades sequence by sequence across the complete pre-neutron
fragmentation space. Primary literature citations and data evaluations are catalogued with DOIs under
[References](references.md).

## Physical framework

In contrast to global statistical treatments—such as the Los Alamos model of Madland and Nix
[MadlandNix1982](@cite), which integrates evaporation spectra over an assumed triangular residual
temperature distribution—the Deterministic Sequential Emission (DSE) framework tracks individual cascade
branches step by step.

### Energy conservation and sequential cascade

For a given pre-neutron fragment `(A, Z)` at initial excitation energy `E*`, de-excitation proceeds
through discrete neutron evaporation stages `k = 1 … k_max`. At each stage `k`, conservation
of energy is evaluated explicitly:

```
Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … k_max
```

where:
- `Ē_r^(k-1)` is the residual excitation energy prior to the `k`-th emission, with boundary condition
  `Ē_r^(0) = E*`;
- `Sₙ^(k-1)` is the neutron separation energy of the emitting precursor nucleus `(A − (k − 1), Z)`;
- `⟨ε⟩_k` is the mean centre-of-mass kinetic energy of the `k`-th evaporated neutron;
- `a_k` and `T_k` are the nuclear level density parameter and residual nuclear temperature of the daughter
  nucleus `(A − k, Z)`, connected via the Fermi-gas equation of state `Ē_r^(k) = a_k T_k²`.

The cascade step equation is solved successively while the available excitation satisfies
`Ē_r^(k-1) ≥ Sₙ^(k-1)` and admits a physical temperature solution (`T_k > 0`), terminating
at `k_max`. Subsequent emissions that remain energetically possible when the emitted neutrons' kinetic
energy is accounted for are evaluated against separation thresholds on `Sₙ` alone. This defines the
maximum allowable sequence count `n ≥ k_max` for configuration `(A, Z, TKE)`, yielding an average
prompt neutron multiplicity:

```
ν(A, Z, TKE) = (n + 1) / 2
```

The residual excitation energy remaining after neutron cascade cessation,

```
E_γ = Ē_r^(n) = E* − Σ_{k=1}^n (Sₙ^(k-1) + ⟨ε⟩_k)
```

is retained in the residual nucleus and dissipated via prompt γ-ray cascades.

### Excitation energy partition at scission

The total excitation energy (TXE) available to a complementary fragment split `(A_L, Z_L) + (A_H, Z_H)` at
total kinetic energy `TKE` is governed by the reaction balance:

```
TXE(A, Z, TKE) = Q(A, Z) + E*_CN − TKE
```

where `Q(A, Z)` is the ground-state reaction Q-value and `E*_CN` is the compound nucleus excitation
energy.

This energy is partitioned between light and heavy fragments using the temperature ratio
`R_T(A_H) = T_L/T_H`, following the Point-by-Point (PbP) formulation of Tudora and Hambsch
[TudoraHambsch2017](@cite):

```
E*_H / TXE = 1 / (1 + R_a R_T²(A_H))
```

where `R_a = a_L / a_H` is the ratio of fragment level density parameters. In this package, `R_T(A_H)` is
obtained from the model-independent inversion of experimental prompt neutron multiplicities `ν(A)`
implemented in [FissionTemperatureRatio.jl](https://github.com/PaulGoG/FissionTemperatureRatio.jl)
[TudoraGogita2024](@cite).

### Evaporation spectra and laboratory kinematics

Neutron evaporation from fully accelerated fragments is modeled using Weisskopf-Ewing statistical theory.
Closed-form centre-of-mass evaporation spectra and their laboratory transformations follow Terrell
[Terrell1959](@cite), Madland and Nix [MadlandNix1982](@cite), Märten and Seeliger [Maerten1984](@cite),
and Kawano et al. [Kawano2021](@cite). Multiplicity weighting and integration across fragmentations follow
the evaluation standards of Capote et al. [Capote2016](@cite).

### Post-neutron product distributions

Because the allowable sequence count `n` is an integer for every configuration `(A, Z, TKE)`,
post-neutron product masses `A_p ∈ [A − n, A]` remain discrete integers without artificial rounding.
- **Emission-probability count (default)**: The cascade terminates on each of its `n` residual stages
  with equal probability `1/n`, yielding `P(A_p = A − k) = 1/n` (`k = 1 … n`). This formulation
  ensures exact mass conservation, so that the mean mass loss per fragment equals the calculated multiplicity:
  `⟨A − A_p⟩ = ν`.
- **Sequence count (`[products] count = "sequence_count"`)**: Imposes `A_p = A − n` with probability one,
  reproducing the convention of Tudora [Tudora2022](@cite), applied to ²³⁹Pu [TudoraFujio2024](@cite) and
  ²³³U [Tudora2025](@cite). While representing maximum potential depletion, it removes `⟨n⟩ ≈ 2ν − 1`
  mass units on average.
- Post-neutron product kinetic energies `E_K_p = E_K A_p/A` and total product kinetic energy
  `TKE_p = E_K_p,L + E_K_p,H` follow the kinematic relations of Tudora [Tudora2022](@cite),
  assuming isotropic emission leaves average fragment velocities unchanged.

### Provenance and literature lineage

The DSE framework was introduced by Tudora, Hambsch, and Tobosaru [TudoraHambschTobosaru2018](@cite),
with systematic residual properties established in [Tudora2020](@cite), experimental validations in
[Tudora2019](@cite), and comparative analysis with Point-by-Point modeling in [Tudora2021](@cite). Nuclear
masses, reaction Q-values, and separation energies are taken from the AME2020 evaluation [Wang2021AME](@cite),
[Huang2021AME](@cite). Isobaric charge distributions follow Wahl's `Zₚ` model [Wahl1988](@cite),
[Wahl2002](@cite). Nuclear level densities are calculated using the Back-Shifted Fermi Gas model of von Egidy
and Bucurescu [EgidyBucurescu2005](@cite), [EgidyBucurescu2009](@cite) and the composite formula of Gilbert
and Cameron [GilbertCameron1965](@cite).

This package is a modern, modular Julia re-implementation of the computational code originally developed
for my 2023 MSc thesis (Faculty of Physics, University of Bucharest, supervised by Prof. Anabella Tudora).
Where current defaults refine earlier formulations (such as mass-conserving emission-probability counts and
per-neutron spectral weighting), options are provided to reproduce my original thesis results identically.

## Package organization

A run is driven by a TOML configuration and writes its tables; see [Running a case](running.md),
which also documents the configuration schema and the simulation pipeline. Standard nomenclature is
defined in [Naming](naming.md). Auxiliary fragmentation properties (masses, charge distributions,
level densities, energy balance, and yields) are provided by the companion domain package, summarised in
[The fragmentation domain](domain.md). Detailed module references cover
[the excitation-energy partition](energetics.md) and
[the evaporation, the emission chain and the spectrum](emission.md).

```@index
```

```@docs
DeterministicSequentialEmission
PACKAGE_ROOT
```

