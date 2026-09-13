# Changelog

Notable changes to DeterministicSequentialEmission.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

The code as submitted with the 2023 MSc thesis is on the `original` branch and is not versioned
here. This changelog begins with the rewrite.

## [Unreleased]

### Changed

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
