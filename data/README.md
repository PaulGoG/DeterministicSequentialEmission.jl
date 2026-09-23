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
| `U233_nth/charge_distribution_vs_A_wahl_systematics.dat` | the same, from the `Zₚ` systematics rather than a fit — **generated, and superseded by the table beside it** | `A dZ sigma_Z` |
| `<system>/R_T_vs_A_H_systematic_trend.csv` | temperature ratio `R_T(A_H)` for the TXE partition | `A_H,R_T,R_T_uncertainty` |
| `<system>/Y_vs_A_TKE.dat` | pre-neutron fragment mass yield `Y(A,TKE)` | `A TKE Y Y_uncertainty` |
| `U235_nth/Y_vs_A_H_TKE.dat` | the same indexed by the heavy fragment alone | `A_H TKE Y Y_uncertainty` |
| `U233_nth/Y_vs_A/*.dat` | pre-neutron mass yield `Y(A)`, a marginal | `A Y` |
| `U233_nth/TKE_vs_A/*.dat` | mean total kinetic energy `⟨TKE⟩(A)`, the other marginal | `A TKE` |
| `U233_nth/sigma_TKE_vs_A_*.dat` | dispersion of `P(TKE|A)`, for the reconstruction | `A sigma_TKE`, MeV |
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

²³³U now carries an evaluated charge distribution of the same construction as the other three, so
it no longer runs on the systematics. Its two marginals sit under their EXFOR accession numbers,
21981005 and 21981008, both Geltenbort 1985, and `Y(A)` is in percent as every other yield here
is. Its `σ_TKE(A)` comes in two processed sets, Baba and Asghar; the configuration reads the first
and the second is staged beside it for comparison.

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
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021),
[doi:10.1088/1674-1137/abddb0](https://doi.org/10.1088/1674-1137/abddb0);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021),
[doi:10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf).
The source is the unrounded mass table `mass.mas20` of the Atomic Mass Data Center, published at
<https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt> (file dated 3 March 2021, MD5
`6d28b75833cf53c7cc230223f63da6f6`). The shipped table carries that file's mass excess and
uncertainty for all 3558 nuclides, including the 1008 values the evaluation marks with `#` as
estimated — the mark itself is not carried, so the table does not tell an estimated mass from a
measured one — after one pass through single precision: every entry is the source rounded to
the nearest `Float32` and printed with eight decimals, verified row by row against the source on
23 September 2026. The pass moves a mass excess by at most 0.006 keV, and by at most 0.004 keV
over 70 ≤ A ≤ 180, where it exceeds the evaluation's own uncertainty for ⁸⁶Kr alone, by 0.04 eV;
elsewhere only for ¹H, ²H, ³H and ³He. No quantity computed here is sensitive to it, and the
table is kept as it is so that the reference values the tests pin keep their meaning.
`HEADERLESS_MASS_EXCESS_SPEC` reads the same five columns without a header line, the form the
table circulated in before the shipped copy was given one.

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

Two different objects share that format and must not be confused, and the difference is not one
of quality alone.

The four `<system>/charge_distribution_vs_A.dat` files — `Cf252_sf`, `U235_nth`, `Pu239_nth` and
`U233_nth` — are **effective Gaussian parameters**: the `Zₚ` model is evaluated per reaction with
its own least-squares parameters, *including* the even-odd factors `F_Z` and `F_N`, and a plain
Gaussian is then fitted to the resulting `p(Z,A)` at each mass. `ΔZ(A)` is the centre of that
fitted Gaussian less `Z_UCD(A)`, and `σ_Z(A)` is its width. Wahl's own least-squares `σ_Z` is a
different quantity — a model parameter that carries no even-odd structure.

That is why these tables oscillate with a period of about five mass units while the systematics
cannot: the even-odd effect is *inside* them. The period follows from the mechanism — `Zₚ`
advances by `Z_F/A_F` per mass unit, so even-charge dominance recurs once it has advanced by two,
every `2A_F/Z_F ≈ 5.1` masses.

**So a plain Gaussian `p(Z,A)` is the right partner for these tables**, which is what this package
uses: the even-odd modulation is already absorbed into the fitted width and centre, and nothing
further has to be applied or renormalized away. Averaged over each system's own measured `Y(A)`
they return

| | `⟨ΔZ⟩` | `⟨σ_Z⟩` | masses |
|---|---|---|---|
| ²⁵²Cf(sf) | −0.492 | 0.649 | 126–169 |
| ²³⁵U(n_th,f) | −0.503 | 0.618 | 118–160 |
| ²³⁹Pu(n_th,f) | −0.515 | 0.614 | 120–160 |
| ²³³U(n_th,f) | −0.499 | 0.615 | 117–160 |

which is where the conventional fallbacks of `|0.5|` and `0.6` come from, and is the sense in
which those fallbacks are not arbitrary.

They are authoritative wherever they exist, and where they do not, the package builds the same
model rather than reaching for the means. It does so in two layers, and the difference between
them is which reaction the parameters were fitted to. Table A of the 1988 paper carries
least-squares parameters for exactly these four reactions; eq. (17) and Table 2 of LA-13928
estimate the same parameters from `Z_F`, `A_F` and the excitation energy for any reaction in
`90 ≤ Z_F ≤ 98`, `230 ≤ A_F ≤ 252`, `PE ≤ 8 MeV`. Fig. 18 of that report quotes a reduced `χ²` of
2.9 for the first against 7.9 for the second, on ²³⁵U(n_th,f). Neither needs a file: both are in
`src/wahl.jl` and `src/wahl_reactions.jl`, and `charge_distribution_file = "wahl"` selects
whichever the system reaches.

`U233_nth/charge_distribution_vs_A_wahl_systematics.dat` is **the systematics**, the second of
those two layers, tabulated here as the comparison that shows what the first layer and the fitted
table add. Its name and its header row both say so. It is generated, not measured, and it does not
replace an evaluation.

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
