# Changelog

Notable changes to DeterministicSequentialEmission.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

The code as submitted with the 2023 MSc thesis is on the `original` branch and is not versioned
here.

## [Unreleased]

### Added

- The sequential emission model: the energy conservation equation
  `Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²` solved sequence by sequence at every fragment
  configuration `(A, Z, TKE)`, under a constant or a variable evaporation cross section.
- Three TXE partitions behind one `partition` interface: scission modelling, `E*_H/TXE (A_H)` and
  `R_T (A_H)`. The temperature-ratio curve is read through the producing run's manifest, which is
  refused unless it declares ordinate `R_T` over abscissa `A_H`.
- Two level density models: Gilbert-Cameron, deformed and undeformed, and the Egidy-Bucurescu
  back-shifted Fermi gas.
- The Wahl `Zₚ` charge distribution, selected by `charge_distribution_file = "wahl"`, in three
  layers after an evaluated table: per-reaction Table A parameters, the eq. (17) systematics, then
  the conventional means `ΔZ = −0.5`, `σ_Z = 0.6`. Even-odd factors are optional and off by
  default; the layer used is recorded in the run metadata.
- Fragment yields `Y(A,Z,TKE)`, their projections, and `factorized_yield`, which reconstructs the
  joint distribution from `Y(A)` and `⟨TKE⟩(A)` with a measured `σ_TKE(A)`.
- The averaging layer, with uncertainties propagated through `Measurements.jl`.
- The prompt fission neutron spectrum, centre-of-mass and laboratory, in closed form.
- `[spectrum] pair_weighting = "multiplicity" | "excitation"` selects the light and heavy
  weights within a pair, and `[spectrum] average = "per_fission" | "per_neutron"` the weight of
  a pair spectrum in the total, `Y` or `Y ν_pair`. The defaults are `multiplicity` and
  `per_neutron`, the form of evaluated spectra and of the Los Alamos and Hauser-Feshbach
  prescriptions; `per_fission` is the thesis formulation. Both keys enter the run identifier and
  the metadata.
- A configuration-driven entry point, `scripts/run.jl`, validating every key. Runs are written to
  `data/sims/<system>/<run>` with `metadata.toml` and copies of the configuration and the resolved
  `Manifest.toml`; a directory that already holds a run is moved aside with a numbered suffix.
- Per-run charge-distribution diagnostic tables: the distribution used, the Wahl systematics where
  applicable, and the conventional means.
- Five shipped configurations: ²⁵²Cf(sf), ²³⁵U(nth,f), ²³⁹Pu(nth,f), ²³³U(nth,f), and
  ²³⁵U(nres,f) at the 580 eV resonance.
- `run_identifier`, built from the configuration keys that change a result.
- Figures in `plotting/`, redrawn from the written tables.
- `check.jl`: the formatting gate, then the test suite.
- A References page in the documentation listing every source with its DOI, through
  `DocumenterCitations.jl`; docstrings cite by DOI as well.

### Changed

- One name per quantity from configuration key to column header, set out in
  `docs/src/naming.md`; configuration schema, data layout and exported names change, with no
  compatibility aliases.
- The prompt neutron spectrum weights a mass-symmetric fragmentation once in its pair sum.
- `SegmentedCurve` is undefined outside its first and last points.
- The `Zₚ` systematics return the lattice dispersion `rms = √(σ_Z² + 1/12)`, the quantity an
  evaluated table tabulates.
- `[system] channel = "sf" | "nth" | "nres" | "nfast"` replaces `reaction`.
- For `method = "temperature_ratio"`, `[partition] subdirectory` and `curve_label` replace a
  direct curve path; `excitation_ratio_file` keeps a direct path.

### Removed

- `[output] subdirectory`; the run directory follows from the system section.

### Validated

- The naming scheme moves no number: all five configurations reproduce their result tables cell
  by cell, to the last written digit.
- Both archived thesis runs reproduced: 9819 and 9788 configurations, all 44 679 and 45 070
  sequences to 5 × 10⁻⁹ MeV, and the yield matrix cell for cell, 48 985 of them.
- The averaged observables `ν(A)`, `T(A)`, `E*(A)` and `E*(TKE)` agree with the archived run to
  between 4 × 10⁻¹⁵ and 7 × 10⁻¹⁴.
- Total multiplicities within −7 % to +5 % of the evaluated values across four fissioning systems,
  with no parameter adjusted per system.
- Between thermal and 580 eV, `ν(A)` moves by at most 5 × 10⁻⁴ neutrons and the total by 10⁻⁷.
- The prompt neutron spectrum satisfies unit area, `⟨E⟩ = E_f + ⟨ε⟩` and the second-moment
  identity over the whole domain.

### Known limitations

- Total multiplicity is 7 % low for ²⁵²Cf and ²³⁹Pu and 4–5 % high for ²³⁵U and ²³³U; undiagnosed.
- The laboratory spectrum is soft for ²⁵²Cf: `⟨E⟩` is 2.01 MeV against a measured 2.13 with the
  per-neutron average, 1.96 MeV with the per-fission one. The remainder is not a weighting
  effect.
- Post-emission distributions are not yet ported.
- Only the low-energy Wahl branch is implemented, so multiple-chance fission is not covered.
