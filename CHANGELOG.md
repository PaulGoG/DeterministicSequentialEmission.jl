# Changelog

Notable changes to DeterministicSequentialEmission.jl. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

The code as submitted with the 2023 MSc thesis is on the `original` branch and is not versioned
here. This changelog begins with the rewrite.

## [Unreleased]

### Added

- The sequential emission model: the energy conservation equation
  `Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²` solved sequence by sequence at every fragment
  configuration `(A, Z, TKE)`, under a constant or a variable evaporation cross section.
- Three TXE partitions — scission modelling, the excitation energy ratio `E*_H/TXE (A_H)`, and
  the temperature ratio `R_T (A_H)` — behind one `partition` interface. `read_segments` reads
  what `FissionTemperatureRatio.jl` writes, so the partition can come from measured `ν(A)`
  rather than from a curve typed into a configuration.
- Two level density prescriptions, Gilbert-Cameron with both the deformed and the undeformed
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
- Five shipped configurations: ²⁵²Cf(sf), ²³⁵U(n_th,f), ²³⁹Pu(n_th,f), ²³³U(n_th,f), and
  ²³⁵U(n,f) on a 580 eV resonance beam.
- Validation figures in `plotting/`, carrying their own environment and reading the written
  tables rather than solving the model again.
- `check.jl`, applying the same formatting gate as CI and then the tests.

### Validated

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
