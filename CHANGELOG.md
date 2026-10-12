# Changelog

Notable changes to DeterministicSequentialEmission.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

The code as submitted with the 2023 MSc thesis is on the `original` branch and is not versioned
here.

## [0.1.1] - 2026-10-12

### Changed

- The README and the documentation cite by DOI the physics the package implements: the
  deterministic sequential emission model, the point-by-point partition and evaporation
  framework it builds on, the extraction of the temperature ratio from measured multiplicities,
  and the spectrum treatments, so that the published formulation is told from this
  implementation of it. The documentation pages cite through the References page, and the
  bibliography gains the 2022 and 2024 papers.

## [0.1.0] - 2026-10-12

### Added

- The sequential emission model: the energy conservation equation
  `Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²` solved sequence by sequence at every fragment
  configuration `(A, Z, TKE)`, under a constant or a variable evaporation cross section.
- Three TXE partitions behind one `partition` interface: scission modelling, `E*_H/TXE (A_H)` and
  `R_T (A_H)`. The temperature-ratio curve is read through the producing run's manifest, which is
  refused unless it declares ordinate `R_T` over abscissa `A_H`.
- The fragmentation domain from
  [FissionFragmentsDomain.jl](https://github.com/PaulGoG/FissionFragmentsDomain.jl), resolved from
  its repository and pinned to a release. It provides nuclides and systems, the shipped AME2020
  table (`mass_excess_file = "ame2020"`), the charge distribution, level density parameters, the
  energy balance, fragment yields and averages, and the temperature-ratio run record.
- A staged `R_T(A_H)` run whose manifest records its `[domain]` is checked against the
  configuration: a different level density model or Gilbert–Cameron branch, number of charges per
  mass, charge model or mass table is refused, naming the key. A curve inverted through an effective
  level density ratio loads with a warning, since the partition here applies `R_T` charge by charge.
  A manifest without the record loads with a warning.
- The shipped configurations take their fragment yields from EXFOR, through ExforFissionData.jl
  (v0.3.0 at the time of writing), in place of staged tables.
  - ²⁵²Cf: the measured matrix of Göök et al. (2014), on its 1-MeV bin centres, so the TKE grid is
    now 130.5–229.5 MeV.
  - ²³⁵U: Al-Adili et al. (2016), with Baba's width.
  - ²³⁹Pu: Wagemans et al. (1984), with the ²⁴⁰Pu(sf) width of Schillebeeckx et al. (1992). With
    Wagemans's marginals it reproduces the measured ²³⁹Pu(nth,f) integral width to 0.02 MeV.
  - ²³³U: Geltenbort (1985), with Baba's width.
  - The run record keeps each input's path.
- A reconstructed yield takes the width of `P(TKE|A)` from `kinetic_energy_dispersion_file`, or
  from a width the configuration states; none is assumed otherwise. The run record says whether
  the width is from the experiment that measured `⟨TKE⟩(A)` or borrowed from another,
  `yield_kinetic_energy_width`, judged by EXFOR entry. It also lists the masses without a width,
  which sit at their mean TKE; a width-free run writes no tables resolved in TKE.
  - ²³³U(nth,f) borrows Baba's (1997) pre-neutron width, EXFOR 23014002, since Geltenbort (1985)
    reports none.
  - Without any width, `⟨ν⟩` moves by +1.2 % (²⁵²Cf) and +3.9 % (²³⁵U).
- `[yield] symmetrize`, true by default, imposes the pre-neutron identities `Y(A) = Y(A₀ − A)`,
  `⟨TKE⟩(A) = ⟨TKE⟩(A₀ − A)` and `σ_TKE(A) = σ_TKE(A₀ − A)` by averaging measured complements. It
  enters the run identifier as `Ysym`.
- The run record compares the input yield's mean TKE with the energy standard of the system,
  `recommended_mean_total_kinetic_energy`, and `scripts/run.jl` prints the offset.
- `charge_distribution_file = "wahl"` takes Wahl's `Zₚ` model directly: the 1988 model of the
  reaction where it exists, which reproduces that evaluation's own tables, else the 2002
  systematics, else the conventional means. An evaluated table named in the configuration
  outranks all three, and the layer used is recorded in the run metadata. The per-run
  charge-distribution diagnostics reduce each model to effective Gaussians. The shipped
  configurations take `"wahl"`, so every input they need traces to a public source: for the
  five systems ⟨ν⟩ moves by +0.03 % to +0.76 % and ν(A) by at most 0.035 over the yield peaks,
  against the supplied tables used before.
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
- `[spectrum] sequence_weighting = "emission_probability" | "uniform"` selects how the emission
  sequences of a fragment are weighted in its spectra and its mean neutron energy. The default,
  `emission_probability`, weights the `k`-th neutron of a chain of `n` by `(n+1−k)/n`, its
  probability of emission under `ν = (n+1)/2`, and averages `⟨ε⟩(A)` per emitted neutron;
  `uniform` is the average of the thesis and of Tudora *et al.* (2018), with `⟨ε⟩(A)` averaged
  per fission. Multiplicities, residual temperatures and excitation energies do not depend on
  it. The key is a token of the run identifier.

### Changed

- The staged temperature-ratio runs are those of FissionTemperatureRatio.jl v0.3.0 on the
  ExforFissionData v0.3.0 retrieval. Its trend covariance takes the error of a pooled dataset as
  correlated along the mass axis, Basova (1979) and Zamyatnin (1979) are pooled as one
  experiment, and Nishio (1995) leaves the ²³⁹Pu pool: the trend `⟨R_T⟩` moves from 1.0878 to
  1.0821 for ²⁵²Cf and from 1.0268 to 1.0388 for ²³⁹Pu, `⟨ν⟩` by −0.0005 and −0.0014, `⟨E⟩` by
  −0.001 and +0.005 MeV; ²³⁵U and ²³³U are unchanged.
- `shell_correction_file = "gilbert_cameron"` takes the Gilbert–Cameron table the domain package
  now ships, Table III of the 1965 paper verified against the print and against Geant4. The five
  shipped configurations name it; a path under the data directory still works.
- FissionFragmentsDomain 0.2.0. Wahl's 1988 model now has `ΔZ = 0` at `A_F/2` for ²⁵²Cf, which
  moves `ν(126)` by 0.025 and `⟨ν⟩` by 7 × 10⁻⁵; no other system changes.

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
- Total multiplicities within 5.5 % of the evaluated values across four fissioning systems, with
  no parameter adjusted per system.
- Between thermal and 580 eV, `ν(A)` moves by at most 5 × 10⁻⁴ neutrons and the total by 10⁻⁷.
- The prompt neutron spectrum satisfies unit area, `⟨E⟩ = E_f + ⟨ε⟩` and the second-moment
  identity over the whole domain.

### Known limitations

- Total multiplicity is 5 % low for ²⁵²Cf, 2 % low for ²³⁹Pu and 5 % high for ²³⁵U and ²³³U. The
  multiplicity rule `ν = (n+1)/2` responds to excitation energy too weakly: at fixed mass
  `−dTKE/dν` is 11.5 MeV per neutron for ²⁵²Cf against 7.3 measured, and the distribution of
  the neutron number is too broad in all four systems.
- The ²⁵²Cf laboratory spectrum follows the measured shape between 0.3 and 5 MeV, lies below the
  data under 0.3 MeV and above them beyond 6 MeV; `⟨E⟩` is 2.16 MeV against a measured 2.13.
  The spectra of the thermal systems relative to that of ²⁵²Cf are reproduced less well than
  under `sequence_weighting = "uniform"`, and `⟨ε⟩` falls too steeply with TKE at fixed mass.
- Post-emission distributions are not yet ported.
- Only the low-energy Wahl branch is implemented, so multiple-chance fission is not covered.
