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
├── input_data/           mass excesses, charge polarization, level density and yield data
├── activate.jl           activates and instantiates the root environment
└── Project.toml
```

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
julia --project=test -e 'using Pkg; Pkg.test()'
```

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
| Yield averaging, post-emission distributions | to port |
| Prompt fission neutron spectrum | to port |
| Figures (CairoMakie) | to port |
| Configuration, provenance, entry point | to port |

Everything marked done is pinned against the archived thesis run: `Q(A_H, Z_H)` over all 245
fragmentations, the level density parameter over all 244 where it is defined, `E*_L` and `E*_H`
over 1200 configurations, and 233 emission sequences with their separation energies,
temperatures, residual energies and mean neutron energies.

The rewrite is pinned against the thesis run archived outside this repository: ²⁵²Cf(SF),
`A_H ∈ [126, 174]`, `TKE = 130:2:230` MeV, 5 charge numbers per `A`, AME2020, variable
evaporation cross section, under both the `E*_H/TXE (A_H)` and `R_T (A_H)` partitions.

## Licence

MIT, see `LICENSE`.
