# Input data

Almost nothing here is version-controlled. The files are third-party scientific data — an
atomic mass evaluation, charge distribution and shell correction systematics, a resonance
parameter library, and experimental fragment mass yields. None of it originates with this
package, and the MIT licence covering the source code does not extend to it: each collection
carries the terms and the attribution of its own source, recorded below.

The atomic mass evaluation is not here: it ships with
[FissionFragmentsDomain.jl](https://github.com/PaulGoG/FissionFragmentsDomain.jl), and a
configuration takes it with `mass_excess_file = "ame2020"`. Its provenance is recorded in that
package's `data/README.md`. The tests that need any of the tables below are skipped when a table
is absent, with a warning naming the file.

## Layout

One directory per fissioning system, named by the system token `<ElementSymbol><A>_<channel>`, and
one directory for the evaluations that belong to no system. A system is one token and one tree:
its yields, its charge distribution and its measured multiplicities sit together.

```
data/
├── reference/                      evaluations that belong to no one system
└── <system>/                       Cf252_sf, U235_nth, U235_nres, U233_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat
    ├── R_T_vs_A_H_systematic_trend.csv
    ├── Y_vs_A_TKE.dat
    └── nu_vs_A/                    one file per measurement
```

| Path | Content | Header |
|---|---|---|
| `reference/shell_corrections_gilbert_cameron.dat` | optional: a Gilbert–Cameron table other than the shipped one | `n S_N S_Z`, MeV |
| `reference/s_wave_resonances_ripl3.dat` | s-wave average resonance parameters, the source of `S₀(A)` | as distributed: commented header, `S₀` in units of 10⁻⁴ |
| `reference/ground_state_deformation_moller.dat` | ground-state quadrupole deformations | `Z A beta2` |
| `<system>/charge_distribution_vs_A.dat` | charge polarization and Gaussian dispersion, evaluated per reaction | `A dZ sigma_Z` |
| `U233_nth/charge_distribution_vs_A_wahl_systematics.dat` | the same, from the `Zₚ` systematics rather than a fit — **generated, and superseded by the table beside it** | `A dZ sigma_Z` |
| `<system>/R_T_vs_A_H_systematic_trend.csv` | temperature ratio `R_T(A_H)` for the TXE partition | `A_H,R_T,R_T_uncertainty` |
| `<system>/Y_vs_A_TKE/relative/*.dat` | measured pre-neutron `Y(A,TKE)`, in counts | `A TKE Y Y_uncertainty` |
| `<system>/Y_vs_A/*.dat` | pre-neutron mass yield `Y(A)`, a marginal | `A Y [Y_uncertainty]` |
| `<system>/TKE_vs_A/*.dat` | mean total kinetic energy `⟨TKE⟩(A)`, the other marginal | `A TKE [TKE_uncertainty]` |
| `<system>/sigma_TKE_vs_A/*.dat` | dispersion of `P(TKE|A)`, pre-neutron | `A sigma_TKE [sigma_TKE_uncertainty]`, MeV |
| `<system>/Y_vs_A_TKE.dat`, `U235_nth/Y_vs_A_H_TKE.dat`, `Pu239_nth/Y_vs_A_TKE_chiba.dat` | staged joint tables of earlier work; no longer read by a shipped configuration | `A TKE Y Y_uncertainty` |
| `U235_nth/E_def_vs_A_Z.dat` | `ΔE_def` for the scission-modelling partition | `A Z E_def`, MeV |
| `<system>/nu_vs_A/*.dat` | measured `ν(A)` per fragment | `A nu nu_uncertainty` |
| `Cf252_sf/spectrum_maxwellian_ratio_vs_E/*.dat` | measured prompt neutron spectrum as a ratio to a Maxwellian | `E spectrum_maxwellian_ratio spectrum_maxwellian_ratio_uncertainty` |

A file name is `<quantity>_vs_<abscissa>`, which is a statement rather than a list of symbols.
Where a quantity was measured more than once, the directory carries the quantity and the file
carries the provenance, `<accession>_<Author>_<year>.dat`. Where one table of a quantity is
authoritative and another is an alternative, the first takes the plain name and the second carries
a qualifier — `Y_vs_A_TKE.dat` and `Y_vs_A_TKE_chiba.dat`.

All files are whitespace-separated with a single header line, except the `.csv` tables, which are
comma-separated. **Readers take columns by position, not by header text**, and the uncertainty
column may be omitted rather than zero-filled. Renaming a header is therefore a no-op for code,
which is what made this layout safe to adopt.

Nucleon-number columns are read through `integer_column`, which accepts an integral value
written in floating point — some of these tables were emitted by a generator and write `118.0`
where they mean `118`.

The shipped configurations take their fragment yields from EXFOR through ExforFissionData.jl
v0.2.0. For each system, `Y(A)` and `⟨TKE⟩(A)` come from one experiment, or its measured matrix
does. A width that experiment did not report is borrowed from another on the same system, or on
the same compound nucleus, and the run record flags it:

| System | `Y(A)` and `⟨TKE⟩(A)` | `σ_TKE(A)` |
|---|---|---|
| ²⁵²Cf(sf) | Göök 2014, the matrix 23268002 (counts on 1-MeV bin centres) | in the matrix |
| ²³⁵U(nth,f), (nres,f) | Al-Adili 2016, 23164005/007 (the thermal reference run) | Baba 1997, 23014003 |
| ²³⁹Pu(nth,f) | Wagemans 1984, 21995028/038 | Schillebeeckx 1992, ²⁴⁰Pu(sf), 22273023 |
| ²³³U(nth,f) | Geltenbort 1985, 21981005/008 | Baba 1997, 23014002 |

The ²³⁹Pu width is that of spontaneous fission of the same compound nucleus. With Wagemans's
`Y(A)` and `⟨TKE⟩(A)` it gives a mass-integrated width of 11.79 MeV, against the 11.81 MeV measured
for ²³⁹Pu(nth,f) itself (EXFOR 22273003). The Asghar (1981) width, EXFOR 21771014, is against
post-neutron mass (`SEC`) and is not used for pre-neutron data. The staged tables of earlier work
remain beside them and can still be named by path.

`U235_nres` is the same target as `U235_nth` on the resonance channel. Its evaluated charge
distribution, its `R_T(A_H)` trend and its yield inputs are the thermal ones and are read from the
`U235_nth` tree rather than copied; only its measured `ν(A)` is its own.

## Where to obtain it

**Retrieved data is preferred wherever it can be had**, in this order:

1. A retrieval from the EXFOR archive through `ExforFissionData.jl`, which writes this layout
   directly and records what it kept and what it excluded.
2. Where the archive holds no joint `Y(A, TKE)` for the chosen experiment, the joint distribution
   is reconstructed from its marginals `Y_vs_A` and `TKE_vs_A`, with a width from
   `sigma_TKE_vs_A`, by `factorized_yield`. ²³⁵U, ²³⁹Pu and ²³³U take this route.
3. A staged table, only where neither of the above is available for a system.

Two of these collections come from companion packages, and that is the supported route: both
write this layout directly, and both record the provenance of what they retrieved.

**Experimental fragment mass yields and prompt neutron multiplicities** —
[`ExforFissionData.jl`](https://github.com/PaulGoG/ExforFissionData.jl) retrieves them from the
IAEA EXFOR archive into `data/<system>/<observable>/` and writes them beside a `retrieval.toml`
naming every dataset it kept or excluded, with the reason. A configuration and that package
reproduce a retrieval exactly, since EXFOR entries are immutable once published.

```bash
git clone --branch v0.2.0 https://github.com/PaulGoG/ExforFissionData.jl
cd ExforFissionData.jl
julia scripts/retrieve.jl config/Cf252_sf_Y_vs_A_TKE.toml
julia scripts/retrieve.jl config/U235_nth_Y_vs_A.toml
julia scripts/retrieve.jl config/U235_nth_TKE_vs_A.toml
julia scripts/retrieve.jl config/U235_nth_sigma_TKE_vs_A.toml
julia scripts/retrieve.jl config/Pu240_sf_sigma_TKE_vs_A.toml
```

```bash
julia scripts/retrieve.jl config/Cf252_sf_nu_vs_A.toml
julia scripts/retrieve.jl config/Cf252_sf_spectrum_maxwellian_ratio_vs_E.toml
```

Read the `retrieval.toml` it writes before using what it retrieved. It records, per dataset, the
reaction code, the units and the qualifiers the archive attached, and those qualifiers matter:
`DERIV` marks a dataset derived from other data rather than measured, `SPA` one averaged over a
fission spectrum rather than the thermal one it is filed under. The sets staged here exclude both.
A handful of `252-Cf` multiplicity subentries are miscoded in the archive — `MASS,PR,NU` without
the `/FRG` that would mark them per fragment — and are included deliberately, having been
established as per fragment from the complement invariance `ν(A) + ν(A₀−A) ≈ ν̄` in the data
itself. Those files hold `ν(A)` like every other file in `nu_vs_A/`, and their header now says so
rather than naming the query they came from.

**The temperature ratio `R_T(A_H)`** for the temperature-ratio TXE partition —
[`FissionTemperatureRatio.jl`](https://github.com/PaulGoG/FissionTemperatureRatio.jl) extracts it
from experimental `ν(A)` and tabulates it at every integer `A_H`.

**Stage the whole results directory of one run, not one table out of it.** A run writes

```
results/<system>/
├── manifest_<run>.toml                          the run record; this package reads it first
├── R_T_vs_A_H_segmented_<label>_<run>.csv       header: A_H,R_T,R_T_uncertainty
├── r_nu_vs_A_H_pivots_<label>_<run>.csv         header: A_H,r_nu,r_nu_uncertainty
├── total_average_R_T_<run>.csv
├── dataset_diagnostics_<run>.csv
└── metadata_<run>.toml
```

Copy that directory into `data/<system>/` and point `[partition] subdirectory` at it;
`[partition] curve_label` then picks one of the manifest's `[[segmented_curve]]` entries by its
`label`. Keep one staged run per system — two manifests in one directory is refused, because
nothing would say which run to take.

Copying one table by hand is what this replaces, and the reason is physical rather than tidy: the
two CSVs above hold different quantities in the same column position, `R_T ≈ 1.2` against
`r_ν ≈ 0.5`, and both are plausible. Staging the pivots file by mistake runs to completion and
reports a number that looks like a result. The manifest is what makes the quantity checkable —
this package refuses a manifest that does not declare `ordinate = "R_T"` and
`abscissa = ["A_H"]`, and resolves every curve through `temperature_ratio_file` and through no
other field. The `columns` the manifest states are documentation: tables are read by position.

Run the extraction with `ratio_averaging = "charge_resolved"` for this consumer; the reason is in
`read_segmented_curve`'s docstring, and the manifest's `[domain]` record is checked on load.

The remaining tables are one-off downloads:

**Gilbert–Cameron shell corrections** — Table III of A. Gilbert and A. G. W. Cameron, *Canadian
Journal of Physics* **43**, 1446 (1965), pp. 1453–1455, doi:10.1139/p65-139. It ships with
FissionFragmentsDomain.jl, verified against the printed table and against Geant4's electronic copy,
and `shell_correction_file = "gilbert_cameron"` selects it. A table of the same layout under this
directory can still be named by path.

**`reference/s_wave_resonances_ripl3.dat`** — s-wave average resonance parameters.
RIPL-3: R. Capote *et al.*, *Nuclear Data Sheets* **110**, 3107 (2009). The file is byte-identical
to <https://www-nds.iaea.org/RIPL-3/resonances/resonances0.dat>, header included. The IAEA also
publishes the library as JSON at <https://github.com/IAEA-NDS/ripl3_json>, and RIPL-4 is in
preparation at <https://github.com/IAEA-NDS/RIPL-4>; neither changes the numbers this file holds.

**`<system>/charge_distribution_vs_A.dat`** — charge polarization `ΔZ(A)` and dispersion
`σ_Z(A)`, `A dZ sigma_Z`. The four tables of `Cf252_sf`, `U235_nth`, `Pu239_nth` and `U233_nth` are
supplied evaluations of documented lineage: Tudora et al., *Phys. Procedia* **64**, 62 (2015),
doi:10.1016/j.phpro.2015.04.008, applying the `Zₚ` model of Wahl, *At. Data Nucl. Data Tables*
**39**, 1 (1988), doi:10.1016/0092-640X(88)90016-2. They are not regenerated here.
`charge_distribution_file = "wahl"` takes Wahl's model directly instead. FissionFragmentsDomain
reproduces the 1988 evaluation's own tables row by row, and records in its `data/README.md` how
far these supplied tables lie from it. `U233_nth/charge_distribution_vs_A_wahl_systematics.dat`
is an earlier tabulation of the systematics, kept for comparison.

**`reference/ground_state_deformation_moller.dat`** — ground-state quadrupole deformations from
the finite-range liquid-drop model. P. Möller, J. R. Nix, W. D. Myers, W. J. Swiatecki, *Atomic
Data and Nuclear Data Tables* **59**, 185 (1995), via RIPL.

**`U235_nth/E_def_vs_A_Z.dat`** — not a source table. It is generated from the ground-state
deformations and the deformation energies at scission, and is the input the scission-modelling
partition needs.

**The fragment mass yields.** `U235_nth/Y_vs_A_H_TKE.dat` is the ²³⁵U(nth,f) measurement of
C. Straede *et al.* and `Pu239_nth/Y_vs_A_TKE_chiba.dat` the ²³⁹Pu(nth,f) Chiba set, as each
file's own header records. `Cf252_sf/Y_vs_A_TKE.dat`, `U235_nth/Y_vs_A_TKE.dat` and
`Pu239_nth/Y_vs_A_TKE.dat` are the joint distributions on the whole `A` domain with an uncertainty
column, and are what those three configurations read. `Cf252_sf/Y_vs_A_TKE.dat` is the ²⁵²Cf(SF)
distribution behind every chapter 4 result of the thesis.

These three are staged tables under step 3 above. Moving them onto the retrieval route needs
`Y_vs_A` and `TKE_vs_A` retrieved for each system, after which the configuration takes
`mass_yield_file` and `mean_kinetic_energy_file` in place of `mass_energy_yield_file`, exactly as
`U233_nth.toml` does. That changes the yield weights and so every averaged observable, so it is a
run to compare rather than a substitution to make silently.

## Citing

Experimental data belongs to the groups that measured it. Any published result derived from the
yields or multiplicities here should cite the measurement, not this package — the EXFOR accession
number that opens each retrieved file name is what makes those citations recoverable, and the
`retrieval.toml` beside it records the reaction code and units.
