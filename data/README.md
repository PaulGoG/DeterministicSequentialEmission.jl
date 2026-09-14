# Input data

Almost nothing here is version-controlled. The files are third-party scientific data — an
atomic mass evaluation, charge distribution and shell correction systematics, a resonance
parameter library, and experimental fragment mass yields. None of it originates with this
package, and the MIT licence covering the source code does not extend to it: each collection
carries the terms and the attribution of its own source, recorded below.

**The one exception is `reference/mass_excess_ame2020.dat`**, which is shipped. It is an openly
published evaluation, it is small, and keeping it means the tests pinned against the archived
reference run — `Q`, the level density parameter, both TXE partitions, the emission sequences —
execute on a bare clone. The tests that need any of the rest are skipped when it is absent, with a
warning naming the file.

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
| `reference/mass_excess_ame2020.dat` | mass excesses, 3558 nuclides — **shipped** | `Z A element mass_excess mass_excess_uncertainty`, keV |
| `reference/shell_corrections_gilbert_cameron.dat` | `S(N)`, `S(Z)` against nucleon number, 11 … 150 | `n S_N S_Z`, MeV |
| `reference/s_wave_resonances_ripl3.dat` | s-wave average resonance parameters, the source of `S₀(A)` | as distributed: commented header, `S₀` in units of 10⁻⁴ |
| `reference/ground_state_deformation_moller.dat` | ground-state quadrupole deformations | `Z A beta2` |
| `<system>/charge_distribution_vs_A.dat` | charge polarization and Gaussian dispersion, evaluated per reaction | `A dZ sigma_Z` |
| `U233_nth/charge_distribution_vs_A_wahl_systematics.dat` | the same, from the `Zₚ` systematics rather than a fit — **generated, not a source table** | `A dZ sigma_Z` |
| `<system>/R_T_vs_A_H_systematic_trend.csv` | temperature ratio `R_T(A_H)` for the TXE partition | `A_H,R_T,R_T_uncertainty` |
| `<system>/Y_vs_A_TKE.dat` | pre-neutron fragment mass yield `Y(A,TKE)` | `A TKE Y Y_uncertainty` |
| `U235_nth/Y_vs_A_H_TKE.dat` | the same indexed by the heavy fragment alone | `A_H TKE Y Y_uncertainty` |
| `U233_nth/Y_vs_A.dat` | pre-neutron mass yield `Y(A)`, a marginal | `A Y` |
| `U233_nth/TKE_vs_A.dat` | mean total kinetic energy `⟨TKE⟩(A)`, the other marginal | `A TKE` |
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

`U235_nres` is the same target as `U235_nth` on the resonance channel. Its evaluated charge
distribution, its `R_T(A_H)` trend and its yield are the thermal ones and are read from the
`U235_nth` tree rather than copied; only its measured `ν(A)` is its own.

## Where to obtain it

**Retrieved data is preferred wherever it can be had**, in this order:

1. A retrieval from the EXFOR archive through `ExforFissionData.jl`, which writes this layout
   directly and records what it kept and what it excluded.
2. Where the retrieval's observable vocabulary has no joint `Y(A, TKE)` — it holds `Y_vs_A` and
   `TKE_vs_A` but no product of the two — the joint distribution is reconstructed from those two
   retrievable marginals by `factorized_yield`. Keeping the mass dependence of `⟨TKE⟩` is what
   makes that worth doing; `U233_nth` takes this route already.
3. A staged table, only where neither of the above is available for a system.

Two of these collections come from companion packages, and that is the supported route: both
write this layout directly, and both record the provenance of what they retrieved.

**Experimental fragment mass yields and prompt neutron multiplicities** —
[`ExforFissionData.jl`](https://github.com/PaulGoG/ExforFissionData.jl) retrieves them from the
IAEA EXFOR archive into `data/<system>/<observable>/` and writes them beside a `retrieval.toml`
naming every dataset it kept or excluded, with the reason. A configuration and that package
reproduce a retrieval exactly, since EXFOR entries are immutable once published.

```bash
git clone https://github.com/PaulGoG/ExforFissionData.jl
cd ExforFissionData.jl
julia --project scripts/retrieve.jl config/Cf252_sf_Y_vs_A.toml
```

```bash
julia --project scripts/retrieve.jl config/Cf252_sf_nu_vs_A.toml
julia --project scripts/retrieve.jl config/Cf252_sf_spectrum_maxwellian_ratio_vs_E.toml
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

Run the extraction with `ratio_averaging = "mean_of_ratios"` for this consumer; the reason is in
`read_segmented_curve`'s docstring.

The remaining tables are one-off downloads:

**`reference/mass_excess_ame2020.dat`** — the 2020 atomic mass evaluation.
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021).
Distributed by the Atomic Mass Data Center, <https://www-nds.iaea.org/amdc/>. The evaluation is
distributed with no header line at all; `HEADERLESS_MASS_EXCESS_SPEC` reads it in that form, and
the copy shipped here carries one.

**`reference/shell_corrections_gilbert_cameron.dat`** — shell corrections `S(N)` and `S(Z)`,
Table III of A. Gilbert and A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965),
pp. 1453–1455, as redistributed in the IAEA Reference Input Parameter Library. The paper
tabulates `S(Z)` only to `Z = 98`; the file is rectangular and pads the column with zeros beyond
that, which the reader discards rather than storing as a correction of zero.

**`reference/s_wave_resonances_ripl3.dat`** — s-wave average resonance parameters.
RIPL-3: R. Capote *et al.*, *Nuclear Data Sheets* **110**, 3107 (2009). The file is byte-identical
to <https://www-nds.iaea.org/RIPL-3/resonances/resonances0.dat>, header included. The IAEA also
publishes the library as JSON at <https://github.com/IAEA-NDS/ripl3_json>, and RIPL-4 is in
preparation at <https://github.com/IAEA-NDS/RIPL-4>; neither changes the numbers this file holds.

**`<system>/charge_distribution_vs_A.dat`** — charge polarization `ΔZ(A)` and the dispersion
`σ_Z(A)` of the isobaric charge distribution. A. C. Wahl, *Atomic Data and Nuclear Data Tables*
**39**, 1 (1988); the model is set out in full in Wahl, *Systematics of Fission-Product Yields*,
LA-13928 (2002), which is openly published by the Los Alamos National Laboratory. That report's
own bibliography gives the 1988 paper as volume **39**, as does the reference list of
IAEA-TECDOC-1168, in which Wahl's later account of the same models appears.

Two different objects share that format and must not be confused. The `Cf252_sf`, `U235_nth` and
`Pu239_nth` files are **per-reaction least-squares fits** of the `Zₚ` model to the data of each
reaction — Fig. 18 of LA-13928 quotes a reduced `χ²` of 2.9 for that against 7.9 for the
systematics — and they are authoritative wherever they exist.
`U233_nth/charge_distribution_vs_A_wahl_systematics.dat` is **the systematics**, eq. (17) and
Table 2 of the report evaluated by this package for a system that has no fit of its own; its name
and its header row both say so, and `wahl_charge_distribution` produces it without a file. It
replaces a single mean `ΔZ = −0.5`, `σ_Z = 0.6` applied at every mass; it does not replace an
evaluation.

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
