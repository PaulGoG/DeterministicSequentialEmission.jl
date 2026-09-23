# Changelog

Notable changes to DeterministicSequentialEmission.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

The code as submitted with the 2023 MSc thesis is on the `original` branch and is not versioned
here. This changelog begins with the rewrite.

## [Unreleased]

### Added

- `charge_distribution_file = "wahl"` resolves the `Zₚ` model in three layers rather than always
  taking the systematics, through the new `build_charge_distribution`. `WAHL_PER_REACTION`
  transcribes Table A of the 1988 evaluation — the parameters fitted by least squares to the
  measured fractional independent yields of ²³⁵U(n_th,f), ²³³U(n_th,f), ²³⁹Pu(n_th,f) and
  ²⁵²Cf(sf) — and `wahl_reaction_parameters` matches a system to them on the fissioning nucleus
  and the entrance channel, so the resonance channel of ²³⁵U finds the thermal fit while
  ²⁴⁰Pu(sf) does not take ²³⁹Pu(n_th,f)'s. A reaction Table A does not cover falls through to the
  systematics of eq. (17), and a nucleus outside its fitted range to the conventional means. The
  layer that answered is written into the distribution's `source` and so into the run metadata.

  The channel is part of the match because the nucleus alone is not enough. ²⁴⁰Pu(sf) and
  ²³⁹Pu(n_th,f) reach the same compound nucleus 6.53 MeV apart, and over that gap eq. (17) moves
  `σ_Z` by 0.071 — twelve times the uncertainty Table A quotes on the fitted value. `nres` joins
  `nth` on the same test and passes it: at the 580 eV resonance ²³⁵U's compound excitation is
  6.546097 MeV against 6.545517 MeV at thermal.

  The motivation is that the two layers are the same model differing only in where its parameters
  came from, and Fig. 18 of LA-13928 puts that difference at a reduced `χ²` of 2.9 against 7.9.
  Measured as the mean absolute difference from the four evaluated tables, the fitted parameters
  give

  | | `ΔZ` systematics → per-reaction | `σ_Z` systematics → per-reaction |
  |---|---|---|
  | ²³⁵U(n_th,f) | 0.0617 → 0.0591 | 0.0645 → 0.0483 |
  | ²³³U(n_th,f) | 0.0540 → 0.0442 | 0.0581 → 0.0444 |
  | ²³⁹Pu(n_th,f) | 0.1122 → 0.0443 | 0.0494 → 0.0281 |
  | ²⁵²Cf(sf) | 0.2618 → 0.1193 | 0.1114 → 0.0173 |

  No shipped configuration changes behaviour: the four reactions Table A covers are the four that
  already read an evaluated table, and a table named in a configuration still outranks the model.
  What changes is what a system without one gets.

  The 1988 evaluation carries no slope on `σ_Z` — one average width outside the `Z = 50` shell
  region and `σ_Z(50)` inside it — so the per-reaction constructor sets `σZSL = 0` and takes the
  three parameters Table A does not carry, `SL50` and the two wing slopes, from the systematics.
  Parameters the evaluation assumed rather than determined are listed in each entry's `assumed`
  field, because they are not measurements and should not be read as such.

  Table A also places its own regions, in footnotes a and c, and `WahlSystematics` now carries
  `shell` and `symmetric` so that a source which supplies parameters can supply their regions with
  them; the systematics leaves both empty and eq. (10) places the regions as before. The two
  disagree most for ²⁵²Cf, where eq. (10) holds `σ_Z(50)` over `A' = 127–140` and the footnote
  gives `119–129`: taking the footnote drops the `σ_Z` difference from that reaction's table by a
  factor of 6.4. It costs a little on ²³³U and ²³⁹Pu, 0.0399 → 0.0444 and 0.0255 → 0.0281, which is
  the price of using one reaction's published regions rather than boundaries fitted across
  reactions, and is worth paying for parameters that are themselves per-reaction. The footnotes
  quote fission-product mass while this package indexes pre-neutron fragment mass; the two differ
  by `ν` at the sawtooth minimum, 0.2 to 0.7 neutrons here, so no shift is applied and a one-unit
  shift would move the summed difference by 6 %.

  Only the low-energy branch is implemented. The intermediate and high-energy branches of the 1980
  and 2002 reports describe multiple-chance fission, which this model does not treat.

### Changed

- `SegmentedCurve` is undefined below its first point as it is above the last, so a fit is
  never extrapolated past the data that supported it. `segment_coverage` reports a low-end gap
  the same way it reports one at the high end. No shipped curve starts above `heavy_mass_min`,
  so no shipped result moves.

- `wahl_charge_distribution` puts the systematics on the footing of the tables it stands in for.
  Wahl eq. (8) of the 1988 evaluation gives `rms = √(σ_Z² + 1/12)`: `σ_Z` is the width of the
  model's Gaussian, `rms` the second moment of the yields over the integer charge lattice, and it
  is `rms` that an evaluated table tabulates and that this package reads as a charge dispersion.
  The systematics returned the first where the second was wanted. Applying the relation lowers the
  mean absolute difference from every evaluated table — 34 % for ²³³U, 24 % for ²³⁵U, 5 % for
  ²⁵²Cf, 4 % for ²³⁹Pu — and makes none of them worse. `wahl_dispersion` still returns the model
  parameter, which is what its equations define.
- ²³³U reads an evaluated charge distribution rather than the Wahl systematics. The table is of
  the same construction as the other three — a plain Gaussian fitted to the `Zₚ` model evaluated
  with per-reaction parameters — which its own signatures confirm: the residual against the
  systematics peaks at lag 5, matching `2A_F/Z_F = 5.09`, and the yield-weighted averages are
  `⟨ΔZ⟩ = −0.499` and `⟨σ_Z⟩ = 0.616`. It restores the even-odd effect the systematics cannot
  carry.
- The ²³³U marginals are named for the EXFOR accessions they came from, 21981005 and 21981008,
  both Geltenbort 1985, and `Y(A)` is in percent rather than a fraction. The reconstruction
  renormalises, so no number moved; what changed is that the measurement is now citable.
- A temperature-ratio curve is taken through the producing run's own record, not through a path
  named in a configuration. `[partition] segmented_curve_file` is replaced, for
  `method = "temperature_ratio"`, by `[partition] subdirectory` — the staged run, under the data
  directory — and `[partition] curve_label`, which picks one of its manifest's
  `[[segmented_curve]]` entries. `method = "excitation_ratio"` keeps a direct path under
  `[partition] excitation_ratio_file`, because no package of this toolchain produces that curve
  and so no run record describes one.
  - The reason is physical rather than tidy. An extraction writes `R_T(A_H)` and the multiplicity
    ratio `r_ν(A_H)` side by side, in the same column of two files; `R_T ≈ 1.2` and `r_ν ≈ 0.5`
    are both plausible, so a curve staged by hand can be the wrong quantity and the model will run
    to completion and report a result nobody can tell is wrong.
  - `read_temperature_ratio_manifest` refuses a manifest that does not declare
    `[run] ordinate = "R_T"` and `[run] abscissa = ["A_H"]`, that repeats a
    `[[segmented_curve]]` label, that gives a curve a `kind` outside
    `MANIFEST_CURVE_KINDS`, or that names one file as both `temperature_ratio_file` and
    `multiplicity_ratio_pivots_file`. Each refusal names the key at fault.
  - `temperature_ratio_path` is the only route from a manifest to a file and resolves
    `temperature_ratio_file`. There is deliberately no counterpart for the pivots file.
  - `[run] columns` is read and kept as documentation. It is never consulted to find a column:
    every table in this package is read by position, and a header is discarded unexamined.
  - A run of one system may still take its curve from another system's manifest, which is what
    the resonance-channel case does with the thermal system's trend. That is recorded in the run
    metadata, not refused.
- `SegmentedCurve` carries the file it was read from, and a run records it as
  `[model] partition_curve`, so which tabulation was partitioned with is answerable from the
  result rather than from whoever staged it.
- `TableSpec` column symbols name their quantity, like the headers they describe:
  `[:A_H, :R_T, :R_T_uncertainty]` and `[:A_H, :excitation_ratio, :excitation_ratio_uncertainty]`
  in place of one spec naming its ordinate `value`, and `Y`, `Y_uncertainty`, `TKE_uncertainty`,
  `E_def` in place of bare `yield`, `uncertainty` and `energy`. No file on disk changed: these are
  internal names, and the headers were already correct.
- `DataTable` is `DelimitedTable` and `read_table` is `read_delimited_table`. `Data` as a prefix
  carried nothing; `Delimited` says which kind of table it is, and so keeps it apart from the
  `Table` suffix that `MassExcessTable` and `ShellCorrectionTable` carry, which asserts a lookup
  keyed by nucleon number.
- "Parameterization" is retired as the noun for a curve, as `SegmentedCurve` already was: a
  curve is fitted, segmented, or tabulated, and which of those a reader is looking at is worth
  saying.

- One name per quantity, from configuration key to column header, applied across the package.
  The vocabulary is shared with `ExforFissionData.jl` and `FissionTemperatureRatio.jl` and is
  set out in `docs/src/naming.md`. This is a breaking change to the configuration schema, the
  data layout and the exported names, with no compatibility aliases — the schema changed
  incompatibly regardless, and no alias can bridge a key that no longer exists.
  - Types: `MassTable` → `MassExcessTable`, `ShellCorrections` → `ShellCorrectionTable`,
    `ChargePolarization` → `ChargeDistribution`, `ExperimentalYield` → `MassEnergyYield`,
    `FragmentYields` → `FragmentYield`, `JointedSegments` → `SegmentedCurve`,
    `ExcludedConfiguration` → `Exclusion`, `ConfigurationValues` → `FragmentQuantity`. A yield is
    named by its index space; `Configuration` keeps its meaning as the validated TOML, and the
    tensor over fission configurations no longer shares the word.
  - Readers return the type they are named for: `read_mass_excess_table`,
    `read_shell_correction_table`, `read_charge_distribution`, `read_mass_energy_yield`,
    `read_segmented_curve`. Accessors are named for the quantity:
    `charge_polarization(distribution, A)` and `charge_dispersion(distribution, A)` in place of
    `polarization_at` and `width_at`.
  - `rms` is retired everywhere for the dispersion of the isobaric charge distribution:
    `sigma_Z` in files and headers, `σ_Z` in Julia. `LEGACY_MASS_EXCESS_SPEC` is
    `HEADERLESS_MASS_EXCESS_SPEC`, which is what it actually describes.
- A fissioning system is declared by an entrance channel, `[system] channel = "sf" | "nth" |
  "nres" | "nfast"`, which replaces `reaction`. The reaction code follows from it through
  `CHANNEL_REACTION` and the compound nucleus through `REACTIONS`, so neither is written down.
  `system_label` returns the one token that names a system — `Cf252_sf`, `U235_nres` — and
  `system_notation` the typeset form for a figure; one name may not mean both.
- Configuration keys: `digits` → `significant_digits` (a rename only — it has always meant
  significant figures here, and no number moved), `prescription` → `model`,
  `segments_file` → `segmented_curve_file`, `deformation_file` → `deformation_energy_file`,
  `kinetic_energy_width` → `kinetic_energy_dispersion`, and the three `[yield]` file keys named
  for the quantity each holds. Configurations are named for the system they declare, and the
  test suite asserts it.
- `data/` is one directory per fissioning system and one for the evaluations that belong to no
  system, replacing the split between a yield tree and an experimental tree. File names are
  `<quantity>_vs_<abscissa>`, headers are `<quantity>_uncertainty`, and file extensions no
  longer carry a nuclide. Every data row survived the move byte-identically; only header lines
  were rewritten, which readers do not consult.
- Results are written to `results/<system>/<run>`, with `nu_vs_A.dat`, `E_star_vs_TKE.dat`,
  `eps_vs_A.dat`, `E_r_vs_A.dat`, `spectrum_vs_E.dat` and `Y_vs_Z.dat` in place of the former
  names.

### Removed

- `[output] subdirectory`. It named the directory a run wrote to, and with the entrance channel
  in the file that directory follows from `target_A`, `target_Z` and `channel`, so the key could
  only repeat the system section or contradict it. It existed because two runs of one target
  could not otherwise be told apart; the channel is what tells them apart now.

### Added

- `wahl_charge_distribution(...; even_odd = true)` folds in the even-odd factors `F_Z` and `F_N`
  that the `Zₚ` model carries and this package otherwise omits, with `wahl_even_odd_factors` and
  `wahl_charge_yield` exposing the pieces. `p(Z,A)` becomes a Gaussian modulated by `F(A)` and
  renormalized — LA-13928 requires the renormalization once `F(A)` is present — and a plain
  Gaussian is fitted to it per mass, which is how an evaluated table is constructed.
  - It reproduces the ripple an evaluated table carries: amplitude 0.36 in `ΔZ` and 0.30 in
    `σ_Z` against 0.32 and 0.30 present, correlating at 0.62 for ²³⁵U, 0.36 for ²³⁹Pu and 0.07 for
    ²⁵²Cf, which follows the even-odd effect's own decline with the mass of the fissioning system.
  - It does **not** close the gap to an evaluated table, because that gap is the systematics' own
    parameter error rather than the missing factors — Fig. 18 of LA-13928 gives reduced `χ²` of
    7.9 for the systematics against 2.9 for a per-reaction fit. Better of two estimates where no
    table exists; not a route to reproducing one.
  - Off by default: a tabulated distribution already carries the effect inside its fitted width,
    so folding the factors in on top of one would count it twice.
- `factorized_yield` takes a measured `σ_TKE(A)` instead of one width everywhere, through
  `[yield] kinetic_energy_dispersion_file` and the new `KineticEnergyDispersion`, which falls back
  to `kinetic_energy_dispersion` at masses the table does not reach — the same tabulated-with-a-
  fallback shape `ChargeDistribution` has. The dispersion of `P(TKE|A)` runs from about 6 to 13
  MeV across a fragment range, so a single width was a summary of a real spread.
- A run writes the isobaric charge distribution it was weighted by, beside the two it could have
  taken instead: `charge_distribution_vs_A.dat`, `charge_distribution_vs_A_wahl_systematics.dat`
  where `is_wahl_applicable` says the system is inside the fitted range, and
  `charge_distribution_vs_A_mean.dat`. Each is the three-column `A dZ sigma_Z` table an
  evaluation is supplied in, over exactly the heavy masses swept, so the three are directly
  comparable and the plotting stage needs nothing but a reader. `plotting/figures.jl` draws them
  as `figures/charge_distribution`.
  - Which distribution a run used has a measurable effect on `ν(A)` and was recorded only as a
    path in a configuration. The diagnostic also shows where a sweep runs past the end of its
    evaluated table onto the flat fallback, which ²⁵²Cf does for `A_H` 170 to 174.
- `is_wahl_applicable`, so a caller offering the systematics as one option among several can
  leave it out where the system is outside `WAHL_VALIDITY` rather than catching the constructor.
- `TemperatureRatioManifest`, `ManifestCurve`, `ManifestSystem`,
  `read_temperature_ratio_manifest`, `staged_manifest`, `curve_labels`, `manifest_curve` and
  `temperature_ratio_path`: the run record of a temperature-ratio extraction, and the resolution
  of one of its curves to the table that holds `R_T`.
- The sequential emission model: the energy conservation equation
  `Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²` solved sequence by sequence at every fragment
  configuration `(A, Z, TKE)`, under a constant or a variable evaporation cross section.
- Three TXE partitions — scission modelling, the excitation energy ratio `E*_H/TXE (A_H)`, and
  the temperature ratio `R_T (A_H)` — behind one `partition` interface. `read_segmented_curve` reads
  what `FissionTemperatureRatio.jl` writes, so the partition can come from measured `ν(A)`
  rather than from a curve typed into a configuration.
- Two level density models, Gilbert-Cameron with both the deformed and the undeformed
  branch, and the Egidy-Bucurescu back-shifted Fermi gas.
- The Wahl `Zₚ` systematics, LA-13928 (2002), low-energy branch, as a fallback charge
  distribution for systems with no evaluated `ΔZ(A)`, `rms(A)`. It replaces a single mean value
  applied at every mass; a per-reaction evaluated table remains authoritative wherever one
  exists, and a generated table says so in its name and header.
- Fragment yields `Y(A,Z,TKE)` and their projections, and `factorized_yield`, which reconstructs
  a joint distribution from the two marginals `Y(A)` and `⟨TKE⟩(A)` where no joint measurement
  exists. Keeping the mass dependence of `⟨TKE⟩` is what makes the reconstruction worth doing: a
  single mean applied at every mass turns the rising heavy branch of `ν(A)` into a falling one.
- The averaging layer: every observable is a contraction of a quantity tabulated over fragment
  configurations against a weight, with uncertainties propagated through `Measurements.jl` so
  that correlated terms stay correlated.
- The prompt fission neutron spectrum, centre-of-mass and laboratory, in closed form through the
  upper incomplete gamma function.
- A configuration-driven entry point: `config/<system>.toml` validated key by key against the
  bounds its own comment declares, `scripts/run.jl`, and results written under
  `results/<system>/` with a metadata file recording the configuration, the package version, the
  commit and the platform. A directory that already holds a run is moved aside rather than
  overwritten.
- Five shipped configurations: ²⁵²Cf(sf), ²³⁵U(nth,f), ²³⁹Pu(nth,f), ²³³U(nth,f), and
  ²³⁵U(nres,f) on a 580 eV resonance beam.
- `run_identifier` and `RUN_IDENTIFIER_ABBREVIATIONS`: a run is named by the configuration keys
  that change its result, abbreviated in one exported table rather than at the point of use.
  `run_identifier` throws on a key with no entry, and the test suite asserts coverage and that no
  two keys share a token. `incident_energy` is covered — it changes the compound excitation, and
  two runs of one system on one channel may differ in nothing else.
- Validation figures in `plotting/`, carrying their own environment and reading the written
  tables rather than solving the model again.
- `check.jl`, applying the same formatting gate as CI and then the tests.

### Validated

- The naming migration moves no number. All five production cases were re-run against the
  migrated tree and every result table compared cell by cell with its predecessor: identical
  throughout, to the last written digit.
- Both archived thesis runs reproduced in their entirety: 9819 and 9788 fragment
  configurations, and every one of the 44 679 and 45 070 emission sequences they contain, each
  quantity to 5 × 10⁻⁹ MeV — exact agreement at the eight decimals the archive holds. Both
  sequence counts match at every configuration. The fragment yield matrix is reproduced cell for
  cell, 48 985 of them.
- The averaged observables close the chain: `ν(A)`, `T(A)`, `E*(A)` and `E*(TKE)` over every
  mass and kinetic energy of the archived run, to between 4 × 10⁻¹⁵ and 7 × 10⁻¹⁴.
- Total multiplicities within −7 % to +5 % of the evaluated values across four fissioning
  systems, with no parameter adjusted per system.
- Below the second-chance fission threshold the model is flat in incident energy: between
  thermal and 580 eV, `ν(A)` moves by at most 5 × 10⁻⁴ neutrons and the total by 10⁻⁷.
- The prompt neutron spectrum has no archived counterpart, so it is validated against the
  identities its derivation implies — unit area, `⟨E⟩ = E_f + ⟨ε⟩`, and the second moment — over
  the whole domain.

### Known limitations

- Total multiplicity is 7 % low for ²⁵²Cf and ²³⁹Pu and 4–5 % high for ²³⁵U and ²³³U. The runs
  are internally consistent and the inputs are the best available; the spread is not diagnosed.
- The laboratory spectrum is about 8 % soft in `⟨E⟩` for ²⁵²Cf, 1.96 MeV against a measured
  2.13, with the deficit above 3 MeV.
- Post-emission distributions are not yet ported.
