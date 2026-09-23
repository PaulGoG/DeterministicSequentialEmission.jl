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

![Prompt neutron multiplicity against fragment mass, four fissioning systems](plots/multiplicity_vs_experiment.png)

The sawtooth is the signature observable, and it is not fitted: the total excitation energy is
shared between the fragments by a temperature ratio `R_T(A_H)` extracted from measured `ν(A)` by
a companion package, the emission chain is then solved configuration by configuration, and the
result is contracted against the measured fragment yield `Y(A, TKE)`. Total multiplicities come
out within −7 % to +5 % of the evaluated values, with no parameter adjusted per system.

![Prompt fission neutron spectrum of 252-Cf as a ratio to a Maxwellian](plots/prompt_neutron_spectrum.png)

The spectrum is the same statement in energy. Model and measurement are both absolute here —
the model spectrum integrates to one neutron by construction, and every retrieved set is
reported against the same 1.42 MeV Maxwellian — so nothing is scaled to anything. The model is
about 8 % soft in `⟨E⟩`, 1.96 MeV against 2.13, and the deficit is in the tail above 3 MeV.

## File structure

```
DeterministicSequentialEmission.jl/
├── src/                  library code behind the public interface
├── config/               one TOML per fissioning system, named for it, validated on load
├── scripts/run.jl        entry point: solve a system and write its tables
├── plotting/             CairoMakie environment; redraws figures from those tables
├── test/                 test suite, own environment (Aqua, JET, ExplicitImports)
├── docs/                 Documenter site, own environment
├── data/                 input data — see data/README.md; only the mass evaluation ships
│   └── sims/             written by a run, data/sims/<system>/<run> — not tracked
├── plots/                the two figures above, and the diagnostic panels, which are not tracked
├── formatter/            pinned JuliaFormatter, so the gate cannot drift
├── check.jl              the gate CI applies: format, then test
├── activate.jl           activates and instantiates the root environment
└── Project.toml
```

## Getting started

```
julia activate.jl                      # root environment
julia test/activate.jl                 # test environment
julia docs/activate.jl                 # documentation environment
julia plotting/activate.jl             # plotting environment
```

Solve one fissioning system:

```
julia scripts/run.jl config/Cf252_sf.toml
```

That writes `data/sims/Cf252_sf/<run>/`: `nu_vs_A.dat`, `T_vs_A.dat`, `E_star_vs_A.dat`, the
residual quantities, the projections against total kinetic energy, the fragment yield
projections, the laboratory spectrum, a table of totals, and the provenance — `metadata.toml`
recording the configuration, the package version, the commit and the platform, beside copies of
the configuration and the resolved `Manifest.toml`. `<run>` is the run identifier,
built from the configuration keys that change the result, so two runs of one system that differ
in any of them land side by side rather than one displacing the other. Nothing is overwritten — a
directory that already holds a run is moved aside with a numbered suffix first.

Redraw the figures from what those runs wrote, without solving anything again:

```
julia plotting/figures.jl
```

That draws the two figures above from every system a run has written, a panel of the quantities
internal to the emission chain — `E*(A)`, `Ē_r(A)`, `⟨T⟩(A)`, `⟨ε⟩(A)` — and a diagnostic of the
charge distribution each run was weighted by, against the two it could have taken instead: the
`Zₚ` model built without a table, and the conventional means. Those have no experimental counterpart, so they
are diagnostics and are not tracked.

The charge-distribution diagnostic is worth reading whenever a system is investigated. It is what
shows where an evaluated table and the systematics part company — which is the uncertainty a
system without an evaluation inherits — and where a run has silently fallen off the end of its
table onto the flat fallback. For ²⁵²Cf it does: the evaluation stops at `A_H = 169` while the
sweep runs to 174, so five masses carry `ΔZ = −0.5`, `σ_Z = 0.6`. They hold 0.09 % of the yield,
so nothing moves, but it was not visible before.

Run the test suite, and build the documentation:

```
julia check.jl --check
julia docs/make.jl
```

`check.jl --check` verifies the formatting without rewriting anything, then runs the test suite
with its static analysis.

Before committing, apply the gate CI applies — format against `.JuliaFormatter.toml` with the
pinned formatter, then run the tests:

```
julia check.jl            # format in place, then test
julia check.jl --check    # report differences without rewriting
```

## Configuration

Everything a run needs is in one file. `config/Cf252_sf.toml` is the annotated reference; the
other four are the same file with the system changed. Each is named for the system it declares,
and the test suite asserts that, so a file cannot drift from what it runs. Each key carries its allowed values or its
bounds as a comment, and the loader enforces exactly what the comment says — a configuration
that cannot be honoured fails on load, naming the key.

| System | Configuration | Note |
|---|---|---|
| ²⁵²Cf(sf) | `Cf252_sf.toml` | the archived thesis run's system |
| ²³⁵U(nth,f) | `U235_nth.toml` | |
| ²³⁹Pu(nth,f) | `Pu239_nth.toml` | |
| ²³³U(nth,f) | `U233_nth.toml` | no joint `Y(A,TKE)` exists, so it is reconstructed from two marginals |
| ²³⁵U(nres,f) at 580 eV | `U235_nres.toml` | the resonance-region beam the GELINA measurements use |

A system is named `<ElementSymbol><A>_<entrance channel>` — the channel as the field spells it, not
as EXFOR codes it — and that one token names the configuration, the input tree under `data/` and
the output tree under `data/sims/`. It is also what separates the thermal run of ²³⁵U from the
resonance one, which the reaction code `n,f` alone cannot.

The last two are where the data runs out rather than where the model does. No joint `Y(A, TKE)`
for ²³³U exists anywhere in EXFOR, so the yield is reconstructed from the two marginals that were
measured, `Y(A)` and `⟨TKE⟩(A)`, from the same experiment, with a measured `σ_TKE(A)`. The
substitution is named in the configuration rather than assumed.

The resonance case exists to establish a null result. At 580 eV the compound excitation is
5.8 × 10⁻⁴ MeV above thermal out of 6.5 MeV, and `ν(A)` moves by at most 5 × 10⁻⁴ neutrons, the
total by 10⁻⁷. Below the second-chance fission threshold the model is flat in incident energy, so
a measurement made on a resonance beam is directly comparable with a thermal one. Above that
threshold multiple-chance fission opens, which this model does not treat.

## Input data

Almost none of it is version-controlled: it is third-party scientific data that the licence on
this code does not extend to. `data/README.md` records every file, its format and its source,
and points at the two companion packages that write this layout directly. The 2020 atomic mass
evaluation is the exception and ships with the package — it is openly published, it is small,
and keeping it means the tests pinned against the archived reference run execute on a bare clone.
Tests needing anything else announce themselves as skipped rather than passing vacuously.

**Retrieved data is preferred wherever it can be had.** A measurement comes from EXFOR through
`ExforFissionData.jl` if the archive holds it; where the retrieval has no joint `Y(A, TKE)`, the
joint distribution is reconstructed by `factorized_yield` from the two marginals it does hold,
`Y(A)` and `⟨TKE⟩(A)`, which is the route ²³³U takes; a staged table is the last resort.
`data/README.md` says which systems sit where.

The isobaric charge distribution is the one input that cannot be retrieved. EXFOR holds the
independent yields it is fitted to — 608 points in 30 subentries for ²³³U — but only four of the
43 heavy masses carry three or more distinct `Z`, which is too sparse to fit a Gaussian per mass
and is precisely why the `Zₚ` model exists. So the model stands in for the measurement, with its
parameters taken from the best-founded source the system reaches:

1. **this reaction's own least-squares parameters**, for the four reactions Table A of the 1988
   evaluation covers — ²³⁵U(n_th,f), ²³³U(n_th,f), ²³⁹Pu(n_th,f), ²⁵²Cf(sf);
2. **the systematics**, eq. (17) of LA-13928, which estimates the same parameters from `Z_F`,
   `A_F` and the excitation energy, for any nucleus in its fitted range;
3. **the conventional means**, `ΔZ = −0.5` and `σ_Z = 0.6` at every mass.

`charge_distribution_file = "wahl"` takes the first of those three that answers, and the run
metadata records which one did. An evaluated table named in the configuration outranks all three.
Everything here is fitted to spontaneous, thermal and resonance fission; above roughly 8 MeV
multiple-chance fission opens and the model's separate high-energy branches, which this package
does not implement, would be needed.

## Companion packages

`R_T(A_H)` need not be supplied by hand. [`FissionTemperatureRatio.jl`][ftr] extracts it from
experimental prompt neutron multiplicity data and tabulates it at every integer `A_H`;
[`ExforFissionData.jl`][exfor] retrieves that data, and the yields and spectra, from the IAEA
EXFOR archive. Clone either from the URL below and run it; neither is a dependency of this
package, and this package reads only what they write.

**A temperature-ratio curve is taken through the producing run's own record, never through a path
typed into a configuration.** Copy that run's results directory — its `manifest_<run>.toml` and
the tables beside it — under `data/<system>/`, and name it:

```toml
[partition]
method = "temperature_ratio"
subdirectory = "U235_nth"
curve_label = "systematic_trend"
```

The manifest declares what its curves tabulate, and this package refuses one declaring anything
but `R_T` against `A_H`. That is worth the indirection: the run writes `R_T(A_H)` and the
multiplicity ratio `r_ν(A_H)` side by side, in the same column of two files, and `R_T ≈ 1.2` and
`r_ν ≈ 0.5` are both plausible numbers. A curve picked by hand can be the wrong quantity, and the
model will run to completion and report a result nobody can tell is wrong. Reading the record
makes the quantity a matter of record too.

Produce that run with `ratio_averaging = "mean_of_ratios"`: this partition re-expands `R_T` with
per-`(A_H, Z_H)` level density parameters and reduces over `Z` afterwards, so it inverts a
`Z`-averaged multiplicity ratio with `⟨a_L/a_H⟩` rather than `⟨a_L⟩/⟨a_H⟩`.

The retrieval records which datasets it kept and which it rejected, with the reason, and that
record is worth reading before using what it wrote. The archive is not uniformly reliable: sets
carrying `DERIV` are derived rather than measured, `SPA` means averaged over a fission spectrum
rather than a thermal one, and a handful of subentries are miscoded and have to be established
as per fragment from the data itself. The comparisons above use only sets that survive those
checks.

[ftr]: https://github.com/PaulGoG/FissionTemperatureRatio.jl
[exfor]: https://github.com/PaulGoG/ExforFissionData.jl

## Provenance

This branch is a rewrite. The code as submitted with the 2023 MSc thesis, and the thesis itself,
are on the [`original`](../../tree/original) branch. `main` reproduces its physics behind a
tested, configured, documented interface; where it deviates, the deviation is deliberate and
recorded.

The rewrite is pinned against the thesis run archived outside this repository: ²⁵²Cf(SF),
`A_H ∈ [126, 174]`, `TKE = 130:2:230` MeV, 5 charge numbers per `A`, AME2020, variable
evaporation cross section, under both the `E*_H/TXE (A_H)` and `R_T (A_H)` partitions.
`Q(A_H, Z_H)` over all 245 fragmentations, the level density parameter over all 244 where it is
defined, `E*_L` and `E*_H` over 1200 configurations under each partition, 662 and 683 emission
sequences, and the fragment yield projections over 97 masses, 41 charges and 101 kinetic
energies are asserted cell by cell in the test suite.

Beyond the pinned sample, **both archived runs have been reproduced in their entirety**: 9819
and 9788 configurations, and every one of the 44 679 and 45 070 emission sequences they contain,
each quantity to 5 × 10⁻⁹ MeV — exact agreement at the eight decimals the archive was written
with. Both sequence counts match at every configuration. The fragment yield matrix is reproduced
cell for cell, 48 985 of them, and its projections to double precision. The averaged observables
close the chain: `ν(A)`, `T(A)`, `E*(A)` and `E*(TKE)` are reproduced over every mass and kinetic
energy to between 4 × 10⁻¹⁵ and 7 × 10⁻¹⁴, and the total averages to every digit the archive's
summary file prints.

The prompt neutron spectrum has no archived counterpart to pin against, so it is validated
instead against the identities its derivation implies — unit area, `⟨E⟩ = E_f + ⟨ε⟩`, and the
second moment — over the whole domain.

## Status

| Component | State |
|---|---|
| Package environment, test and docs environments, CI | done |
| Tabular input contract, `ExforFissionData.jl`-compatible | done |
| Mass data, `Q`, `S_n`, `TXE` with uncertainties | done |
| Fragmentation domain, `Z_p(A)`, `p(Z,A)` | done |
| Charge distribution — evaluated tables, then `Zₚ` per-reaction, systematics, means | done |
| One name per quantity, from configuration key to column header (`docs/src/naming.md`) | done |
| Level density parameter — Gilbert-Cameron, Egidy-Bucurescu BSFG | done |
| TXE partition — scission modelling, `E*_H/TXE (A_H)`, `R_T (A_H)` | done |
| Sequential equations — constant and variable evaporation cross section | done |
| Fragment yields `Y(A,Z,TKE)`, its projections, and its reconstruction from marginals | done |
| Yield averaging — contractions of `q(A,Z,TKE)`, with uncertainty | done |
| Prompt fission neutron spectrum — centre-of-mass and laboratory, in closed form | done |
| Configuration, provenance, entry point | done |
| Figures (CairoMakie), redrawn from written tables | done |
| Post-emission distributions | to port |

## Licence

MIT, see `LICENSE`.
