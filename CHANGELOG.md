# Changelog

Notable changes to DeterministicSequentialEmission.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

The code as submitted with the 2023 MSc thesis is on the `original` branch and is not versioned
here.

## [Unreleased]

### Added

- `yield_within_domain` in `totals.dat` and in the run metadata: the share of the input yield
  whose fragment masses lie in the mass domain, recorded because what lies outside it is
  renormalised away without a trace: of the shipped inputs, 0.60 % of the ²³⁹Pu yield, beyond
  `A_H = 159`, 0.034 % of ²³⁵U, 0.009 % of ²⁵²Cf and none of ²³³U.

### Changed

- `plotting/figures.jl` draws ν(A), ⟨ε⟩(A) and the product mass yield of each system over its
  own mass window, the masses where the run's `Y(A)` holds at least 0.01 % of the fissions,
  in place of the shared 76–176, each panel labelling its own mass axis; the ordinate limits of
  the ν(A) and ⟨ε⟩(A) grids ignore error bars wider than 0.5 neutrons and 0.3 MeV.
- The staged inputs are the ExforFissionData.jl v0.6.0 retrieval and the FissionTemperatureRatio.jl
  v0.4.1 runs. Every table of v0.5.0 is byte-identical through 0.5.1, 0.5.2 and 0.6.0; the records
  carry the publication windows and ranges of 0.5.1, the stated prompt-γ conventions of 0.5.2 and
  the row-by-row limits of 0.6.0, which admits eight post-neutron yield sets the earlier versions
  refused on the kind of their first row: the independent yields of Dubey 2020 (²³⁵U, ²³⁹Pu) and
  the chain yields of Erdal 1969 (²³⁵U) enter the product figures, the other five hold no point on
  the elements drawn. The FissionTemperatureRatio.jl manifests record
  `zero_polarization_at_symmetry`. The skip reasons of the prompt-γ energy figure name the two
  Billnert 2013 measurements the entry of Oberstedt 2015 reprints.

## [0.2.0] - 2026-10-15

### Breaking

- Results change by default. Each `P(TKE|A,Z)` of the fragment yield is bounded at the physical
  limit of its split, `TKE < Q + E*_CN`, and the measured `⟨TKE⟩(A)` is restored over the
  retained cells (`[yield] kinetic_energy_limit = "physical"`, `kinetic_energy_moments = "mean"`).
  The solved `⟨TKE⟩` now equals the measured input in every system, where it fell 0.12–0.22 MeV
  short for the thermal ones; `⟨TXE⟩` falls by as much and `⟨ν⟩` by 0.9 % (²³⁵U, 2.5353 →
  2.5129), 0.4 % (²³⁹Pu, 2.8132 → 2.8019) and 0.6 % (²³³U, 2.6313 → 2.6156), the mean laboratory
  energy by 0.1–0.3 %. ²⁵²Cf(sf), whose measured matrix keeps its own row means, moves by
  2 × 10⁻⁴ in `⟨ν⟩`. `kinetic_energy_limit = "grid"` reproduces the 0.1.x results.
- Under the physical bound a yield reconstructed from marginals is refused when
  `kinetic_energy_max + kinetic_energy_step` falls short of the largest `Q + E*_CN` of the
  domain; the error names the least `kinetic_energy_max` that suffices. The thermal
  configurations of 0.1.1, at 200 MeV for ²³⁵U and ²³³U and 210 MeV for ²³⁹Pu, are refused: raise
  `kinetic_energy_max`, as the shipped files do to 205 and 216 MeV, or select `"grid"`.
- Run identifiers gain the tokens `pc`, for `[products] count`, and `TKEb`, for the bound of the
  yield on TKE, so new runs land beside the 0.1.x runs rather than replacing them;
  `plotting/figures.jl` draws the most recent run of each system.
- `plotting/figures.jl` writes to `results/`, which replaces `plots/`.
- FissionFragmentsDomain.jl 0.2.4 is required.

No configuration key and no output table is renamed or removed, and a configuration without
`[products]` takes its defaults.

### Added

- The post-neutron product distributions: the mass yield `Y_vs_A_p.dat`, the independent yields
  `Y_vs_Z_A_p.dat` in the columns `Z A_p Y`, `Y_vs_N_p.dat`, the product kinetic energies
  `E_K_p_vs_A_p.dat` and `E_K_p_vs_Z.dat`, and `TKE_p_vs_A_p.dat` against the heavy product
  mass, normalised to the total of the input yield; `mean_product_mass_loss` and `mean_TKE_p`
  in `totals.dat`. A table is named, and its columns ordered, as ExforFissionData.jl names the
  same quantity, so that model and measurement share one layout.
- `[products]`, with `enabled` and `count`. `"emission_probability"`, the default, ends a chain of
  `n` possible sequences on each of its residuals with probability `1/n`, so that the products
  lose the multiplicity; `"sequence_count"` subtracts all `n`, the published prescription, and
  reproduces the post-neutron mass yields of the two ²⁵²Cf(sf) runs of the thesis.
  `ProductCount`, `EmissionProbabilityCount`, `SequenceCount` and `product_masses` implement it.
- The energy a fragment's neutrons leave for prompt γ rays, `gamma_energy`: the excitation less
  what the count's neutrons carry, `E* − Σ_k P(N ≥ k)(Sₙ^(k) + ⟨ε⟩_k)`, the emissions beyond
  `k_max`, which carry no temperature, charged their separation energy alone. A run writes it
  against fragment mass, `E_gamma_vs_A.dat`, and summed over the two fragments of a split against
  the mass of either and against TKE, `E_gamma_tot_vs_A.dat` and `E_gamma_tot_vs_TKE.dat`, the
  names ExforFissionData.jl gives the measured quantities; `E_gamma_per_fragment` and
  `E_gamma_per_fission` go to `totals.dat`, and the second to the run metadata as
  `E_gamma_per_fission_MeV`.
- `[yield] kinetic_energy_limit` and `kinetic_energy_moments`: `"physical"` bounds each
  `P(TKE|A,Z)` at `TKE < Q + E*_CN` and restores the measured moments over the retained cells,
  `"mean"` the mean, `"mean_and_width"` `σ_TKE(A)` too and `"none"` neither, through
  `bounded_fragment_yield` of FissionFragmentsDomain.jl 0.2.4; `"grid"` is the bound of 0.1.x. A
  bounded run writes `TKE_bound_vs_A.dat`, and every run records the solved `⟨TKE⟩` and `σ_TKE`
  in its metadata.
- The `[domain]` check of a staged temperature-ratio run refuses a curve extracted with
  `zero_polarization_at_symmetry = true`, which FissionFragmentsDomain.jl 0.2.4 records.
- `results/`, with a README that shows the model against measurement: every figure beside the
  datasets it draws, each with its EXFOR accession and reference, and what the runs write that
  is not drawn yet.
- The figures of the post-neutron products: the mass yield `Y(A_p)` of the four systems against
  the evaluated independent yields of ENDF/B-VIII.0 and JEFF-3.3, summed over charge, with the
  run's own pre-neutron `Y(A)` as a guide; the elemental yield `Y(Z)` with the global proton
  even–odd effect of model and evaluations; the isotopic distributions of the four most abundant
  light and heavy elements, one figure per system; the neutron-number distribution `Y(N_p)`; and
  the product kinetic energies `E_K_p(A_p)`, beside the pre-neutron `E_K(A)`, `E_K_p(Z)`, and
  `TKE_p(A_p)`, beside the pre-neutron `⟨TKE⟩(A)`.
- The measured post-neutron products of EXFOR in those figures: the product mass yields
  `MASS,SEC,FY` by name and the chain yields `MASS,CHN,FY` pooled in `Y(A_p)`, the independent
  yields `ELEM/MASS,IND,FY` pooled in the isotopic distributions, and the product kinetic
  energies `MASS,SEC,KE` in `E_K_p(A_p)`. Each table is converted from the unit its
  `retrieval.toml` states, from fragments per fission to per cent for the yields, and a unit the
  script does not know stops it. Sets qualified `SPA` or `DERIV` are drawn with open markers;
  sets the record marks superseded, and the ²³⁵U mass yield of Koblik (2014), are left out, and
  the sidecars name them.
- Two figures of the prompt γ energy against the γ observables of ExforFissionData.jl v0.5.0:
  `results/prompt_gamma_energy_per_fission`, the total per fission of the four systems against
  every measured value, one row per dataset, the sets coded partial in γ energy or marked
  preliminary open, the dependent members of a group and a range-limited set left out; and
  `results/prompt_gamma_energy_cf252`, ²⁵²Cf(sf) against TKE, with both slopes fitted, and against
  fragment mass, per fragment and per fission.
- The mean centre-of-mass neutron energy against fragment mass, `⟨ε⟩(A)`, for the four systems
  against the retrieved centre-of-mass means.
- The spectra of the three thermal systems as ratios of the spectrum per fission to that of
  ²⁵²Cf(sf), which is how they were measured, with the model carrying its own multiplicities.
  Sets tabulated in the inverse orientation are inverted, the energy ranges of one EXFOR entry
  are joined, and superseded sets are skipped, as the retrieval record states.
- `scripts/stage_product_yields.jl`, which stages the evaluated independent yields of
  ENDF/B-VIII.0 and JEFF-3.3 from the public ENDF-6 files with their provenance, under
  `data/<system>/Y_vs_A_p/` and, resolved in charge, `data/<system>/Y_vs_Z_A_p/`
  (`Z A_p Y Y_uncertainty`).
- `plotting/figures.jl` takes the directory it writes to as an optional third argument.

### Changed

- The staged inputs are the ExforFissionData.jl v0.5.0 retrieval and the runs of
  FissionTemperatureRatio.jl v0.4.0 on its v0.4.0 retrieval. Every table the shipped
  configurations read holds the numbers it held in 0.1.1; the run records and the names of the
  temperature-ratio files change, the ²⁵²Cf spectrum of Blinov (1980) loses its 0.2209 MeV point,
  an archive defect, and the post-neutron and prompt-γ observables are added.
- The per-system comparison figures are laid out on a two-column grid rather than stacked, each
  panel with its own legend above it. In the product figures the legends of a grid row share one
  row of the layout, so that its axes align whatever the number of entries, and the row grows by
  each legend row beyond two; the measured sets there take colours other than the vermilion and
  blue of the evaluations.

### Fixed

- The evaluated multiplicities of the comparison figures are the prompt ones throughout,
  ENDF/B-VIII.0 MF1 MT456 at 0.0253 eV and, for ²⁵²Cf(sf), its decay-sublibrary constant 3.7590:
  ²³³U was compared with a total multiplicity (2.4968, now 2.4785) and ²³⁹Pu with 2.8750 (now
  2.8680), and ²⁵²Cf took 3.7573 from another evaluation.

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

[Unreleased]: https://github.com/PaulGoG/DeterministicSequentialEmission.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/PaulGoG/DeterministicSequentialEmission.jl/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/PaulGoG/DeterministicSequentialEmission.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/PaulGoG/DeterministicSequentialEmission.jl/releases/tag/v0.1.0
