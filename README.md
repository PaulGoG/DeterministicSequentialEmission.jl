# DeterministicSequentialEmission.jl

Prompt neutron emission from fission fragments, treated sequence by sequence.

Where global treatments integrate an evaporation spectrum over a residual temperature
distribution and return averages, the DSE model solves the energy conservation equation of every
emission sequence separately:

```
Ē_r^(k-1) − S_n^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … n(A,Z,TKE)
```

with Weisskopf-Ewing evaporation spectra in the centre-of-mass system and fragment level
densities in the Fermi-gas regime. Because the number of sequences `n` is an integer at every
configuration, the post-neutron fragment mass `A_p = A − n` is an integer too, without rounding —
which is what makes the model suitable for independent fission product yields.

**This branch is a rewrite in progress.** The code as submitted with the 2023 MSc thesis, and
the thesis itself, are on the [`original`](../../tree/original) branch. `main` reproduces its
physics behind a tested, configured, documented interface; where it deviates, the deviation is
deliberate and recorded.

## File structure

```
DeterministicSequentialEmission.jl/
├── src/                  library code behind the public interface
├── test/                 test suite, own environment (Aqua, JET, ExplicitImports)
├── docs/                 Documenter site, own environment
├── data/                 input data — see data/README.md; only the mass evaluation ships
├── activate.jl           activates and instantiates the root environment
└── Project.toml
```

## Input data

Almost none of it is version-controlled: it is third-party scientific data that the licence on
this code does not extend to. `data/README.md` records every file, its format and its source,
and points at the two companion packages that write this layout directly. The 2020 atomic mass
evaluation is the exception and ships with the package — it is openly published, it is small,
and keeping it means the tests pinned against the archived reference run execute on a bare
clone. Tests needing anything else announce themselves as skipped rather than passing vacuously.

Files still at the repository root — `main_DSE.jl`, `aux_func.jl`, `dse_eq_solvers.jl` and the
rest, plus `Data conversion/` — are the prototype, retained for the duration of the port and
removed as their content moves under `src/`.

## Getting started

```
julia activate.jl                      # root environment
julia test/activate.jl                 # test environment
julia docs/activate.jl                 # documentation environment
```

Run the test suite:

```
julia --project -e 'using Pkg; Pkg.test()'
```

`Pkg.test()` resolves `test/Project.toml` itself, so it is run from the root environment, not
from the test one.

Build the documentation:

```
julia --project=docs docs/make.jl
```

Running the pipeline is not yet wired on this branch; use the `original` branch for that, or
wait for the port. The entry point will be a thin script over a validated TOML configuration.

## Status

| Component | State |
|---|---|
| Package environment, test and docs environments, CI | done |
| Tabular input contract, `ExforFissionData.jl`-compatible | done |
| Mass data, `Q`, `S_n`, `TXE` with uncertainties | done |
| Fragmentation domain, `Z_p(A)`, `p(Z,A)` | done |
| Level density parameter — Gilbert-Cameron, Egidy-Bucurescu BSFG | done |
| TXE partition — scission modelling, `E*_H/TXE (A_H)`, `R_T (A_H)` | done |
| Sequential equations — constant and variable evaporation cross section | done |
| Fragment yields `Y(A,Z,TKE)` and its projections | done |
| Prompt fission neutron spectrum — centre-of-mass and laboratory, in closed form | done |
| Post-emission distributions | to port |
| Figures (CairoMakie) | to port |
| Configuration, provenance, entry point | to port |

Everything marked done is pinned against the archived thesis run: `Q(A_H, Z_H)` over all 245
fragmentations, the level density parameter over all 244 where it is defined, `E*_L` and `E*_H`
over 1200 configurations under each of the two archived partitions, 662 and 683 emission
sequences spanning the domain under each, and the fragment yield projections over 97 masses,
41 charges and 101 kinetic energies.

Beyond the pinned sample, **both archived runs have been reproduced in their entirety**: 9819
and 9788 configurations, and every one of the 44 679 and 45 070 emission sequences they contain,
each quantity to 5 × 10⁻⁹ MeV — exact agreement at the eight decimals the archive was written
with. Both sequence counts match at every configuration. The fragment yield matrix is reproduced
cell for cell, 48 985 of them, and its projections to double precision.

The prompt neutron spectrum has no archived counterpart to pin against, so it is validated
instead against the identities its derivation implies — unit area, `⟨E⟩ = E_f + ⟨ε⟩`, and the
second moment — over the whole domain.

The rewrite is pinned against the thesis run archived outside this repository: ²⁵²Cf(SF),
`A_H ∈ [126, 174]`, `TKE = 130:2:230` MeV, 5 charge numbers per `A`, AME2020, variable
evaporation cross section, under both the `E*_H/TXE (A_H)` and `R_T (A_H)` partitions.

## Companion packages

`R_T (A_H)` need not be supplied by hand. [`FissionTemperatureRatio.jl`][ftr] extracts it from
experimental prompt neutron multiplicity data and tabulates it at every integer `A_H`;
[`ExforFissionData.jl`][exfor] retrieves that data from the IAEA EXFOR archive. `read_segments`
reads what the first of those writes, so a temperature ratio can come from measurement rather
than from a fitted curve typed into a configuration:

```julia
segments = read_segments("R_T_parameterized_<set>_<run>.csv")
method = TemperatureRatio(segments, BackShiftedFermiGas(masses))
```

Produce that file from a run with `ratio_averaging = "mean_of_ratios"`: this partition re-expands
`R_T` with per-`(A_H, Z_H)` level density parameters and reduces over `Z` afterwards, so it
inverts a `Z`-averaged multiplicity ratio with `⟨a_L/a_H⟩` rather than `⟨a_L⟩/⟨a_H⟩`.

[ftr]: https://github.com/PaulGoG/FissionTemperatureRatio.jl
[exfor]: https://github.com/PaulGoG/ExforFissionData.jl

## Licence

MIT, see `LICENSE`.
