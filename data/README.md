# Input data

Almost nothing here is version-controlled. The files are third-party scientific data: the
measurements retrieved from EXFOR (fragment yields and kinetic energies, prompt neutron
multiplicities, spectra and mean energies, post-neutron product yields and kinetic energies,
prompt γ observables), the evaluated independent fission product yields of ENDF/B-VIII.0 and
JEFF-3.3, a resonance parameter library, ground-state deformations, the tables among the inputs
of my 2023 MSc thesis, and the temperature-ratio runs FissionTemperatureRatio.jl extracts from the
measured multiplicities. None of the data originates with this package, and the MIT licence
covering the source code does not extend to it: each collection carries the terms and the
attribution of its own source, recorded below.

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
└── <system>/                       Cf252_sf, U235_nth, U235_nres, U233_nth, Pu239_nth, Pu240_sf
    ├── <observable>/               one retrieved observable: a file per measurement, retrieval.toml
    ├── temperature_ratio/          one staged FissionTemperatureRatio.jl run: its manifest and tables
    └── charge_distribution_vs_A.dat, Y_vs_A_TKE.dat    tables from the thesis inputs, optional
```

`Pu240_sf/` is not a shipped system: it holds the width `σ_TKE(A)` of ²⁴⁰Pu(sf) that the
²³⁹Pu(n_th,f) configuration borrows, and the joint `Y(A,TKE)` of Demattè (1997), EXFOR 22413013,
which no configuration reads.

| Path | Content | Header |
|---|---|---|
| `reference/shell_corrections_gilbert_cameron.dat` | optional: a Gilbert–Cameron table other than the shipped one | `n S_N S_Z`, MeV |
| `reference/s_wave_resonances_ripl3.dat` | s-wave average resonance parameters, the source of `S₀(A)` | as distributed: commented header, `S₀` in units of 10⁻⁴ |
| `reference/ground_state_deformation_moller.dat` | ground-state quadrupole deformations | `Z A beta2` |
| `<system>/charge_distribution_vs_A.dat` | charge polarization and Gaussian dispersion, evaluated per reaction | `A dZ sigma_Z` |
| `<system>/temperature_ratio/` | the staged FissionTemperatureRatio.jl run: `manifest_<run>.toml`, the `R_T(A_H)` curves it names, its diagnostics and metadata | curves: `A_H,R_T,R_T_uncertainty` |
| `<system>/Y_vs_A_TKE/relative/*.dat` | measured pre-neutron `Y(A,TKE)`, in counts | `A TKE Y Y_uncertainty` |
| `<system>/Y_vs_A/*.dat` | pre-neutron mass yield `Y(A)`, a marginal | `A Y [Y_uncertainty]` |
| `<system>/TKE_vs_A/*.dat` | mean total kinetic energy `⟨TKE⟩(A)`, the other marginal | `A TKE [TKE_uncertainty]` |
| `<system>/sigma_TKE_vs_A/*.dat` | dispersion of `P(TKE|A)`, pre-neutron | `A sigma_TKE [sigma_TKE_uncertainty]`, MeV |
| `<system>/Y_vs_A_TKE.dat`, `U235_nth/Y_vs_A_H_TKE.dat`, `Pu239_nth/Y_vs_A_TKE_chiba.dat` | joint yields from the thesis inputs; read by no shipped configuration | `A TKE Y Y_uncertainty`; `A_H TKE Y Y_uncertainty` and `A TKE Y` for the second and third |
| `U235_nth/E_def_vs_A_Z.dat` | `ΔE_def` for the scission-modelling partition | `A Z E_def`, MeV |
| `<system>/nu_vs_A/*.dat` | measured `ν(A)` per fragment | `A nu nu_uncertainty` |
| `<system>/nu_vs_A_TKE/*.dat` | measured `ν` against fragment mass and total kinetic energy | `A TKE nu nu_uncertainty` |
| `<system>/eps_vs_A/*.dat` | mean centre-of-mass energy of the prompt neutrons against pre-neutron fragment mass | `A eps eps_uncertainty`, MeV |
| `<system>/eps_vs_A_TKE/*.dat` | the same against mass and total kinetic energy; ²⁵²Cf and ²³³U only | `A TKE eps eps_uncertainty` |
| `Cf252_sf/T_eff_vs_A/*.dat`, `Cf252_sf/lambda_vs_A/*.dat` | effective temperature and exponent of the form fitted to the centre-of-mass spectrum of a fragment's neutrons, against fragment mass; the record states the form | `A T_eff T_eff_uncertainty`, MeV; `A lambda lambda_uncertainty` |
| `<system>/P_vs_nu/*.dat` | probability of emitting `ν` prompt neutrons in a fission, as tabulated | `nu P P_uncertainty` |
| `<system>/spectrum_vs_E/*.dat` | measured prompt neutron spectrum in the laboratory frame | `E spectrum spectrum_uncertainty` |
| `<system>/spectrum_maxwellian_ratio_vs_E/*.dat` | the spectrum as a ratio to a Maxwellian, whose temperature the record states | `E spectrum_maxwellian_ratio spectrum_maxwellian_ratio_uncertainty` |
| `<system>/spectrum_cf252_ratio_vs_E/*.dat` | the ratio of the system's spectrum and that of ²⁵²Cf(sf) at the same outgoing energy, in the orientation the record states | `E spectrum_cf252_ratio spectrum_cf252_ratio_uncertainty` |
| `<system>/Y_vs_A_p/<library>.dat` | evaluated post-neutron product mass yield, the independent yields summed over charge; `ENDF-B-VIII.0` and `JEFF-3.3`, staged by `scripts/stage_product_yields.jl` with a `provenance.toml` beside them | `A_p Y Y_uncertainty`, per cent per fission, 200 over both products |
| `<system>/Y_vs_Z_A_p/<library>.dat` | the same resolved in charge: the evaluated independent yields before β decay, isomeric states summed | `Z A_p Y Y_uncertainty` |
| `<system>/Y_vs_A_p/*.dat` | measured post-neutron product mass yield, EXFOR `MASS,SEC,FY`, and chain yields `MASS,CHN,FY`, retrieved beside the evaluations | `A_p Y Y_uncertainty`, fragments per fission, 2 over both products |
| `<system>/Y_vs_Z_A_p/*.dat` | measured independent yields `ELEM/MASS,IND,FY`, isomeric states as the record states | `Z A_p Y Y_uncertainty`, fragments per fission |
| `<system>/E_K_p_vs_A_p/*.dat` | measured mean kinetic energy of the products against product mass, `MASS,SEC,KE`; none for ²⁵²Cf | `A_p E_K_p E_K_p_uncertainty`, MeV |
| `<system>/M_gamma_bar/`, `E_gamma_tot/`, `eps_gamma/` | measured prompt γ multiplicity, total γ energy and mean photon energy per fission, one row to a file | `M_gamma_bar`, `E_gamma_tot`, `eps_gamma`, each with its uncertainty where one is stated; MeV |
| `Cf252_sf/` and `U235_nth/`: `E_gamma_vs_A/`, `E_gamma_tot_vs_A/`, `E_gamma_tot_vs_TKE/`, `E_gamma_tot_vs_nu/`, `M_gamma_vs_A/`, `M_gamma_bar_vs_A/`, `M_gamma_bar_vs_TKE/`, `eps_gamma_vs_A/` | the same resolved against fragment mass, TKE or neutron number, per fragment or per fission; *Prompt γ rays* below says which system holds which | abscissa first, e.g. `TKE E_gamma_tot E_gamma_tot_uncertainty`; MeV |

A file name is `<quantity>_vs_<abscissa>`, which is a statement rather than a list of symbols.
Where a quantity was measured more than once, the directory carries the quantity and the file
carries the provenance, `<accession>_<Author>_<year>.dat`. Where one table of a quantity is
authoritative and another is an alternative, the first takes the plain name and the second carries
a qualifier — `Y_vs_A_TKE.dat` and `Y_vs_A_TKE_chiba.dat`.

All files are whitespace-separated with a single header line, except the `.csv` tables, which are
comma-separated. **Readers take columns by position, not by header text**, and the uncertainty
column may be omitted rather than zero-filled. Renaming a header is therefore a no-op for code.

Yields are taken in the normalisation the retrieval writes (per fission, percent or counts, as
its `retrieval.toml` states) and renormalised to `[yield] total`, so only their shape enters.

Nucleon-number columns are read through `integer_column`, which accepts an integral value
written in floating point — some of these tables were emitted by a generator and write `118.0`
where they mean `118`.

The shipped configurations take their fragment yields from EXFOR through ExforFissionData.jl
v0.5.0. For each system, `Y(A)` and `⟨TKE⟩(A)` come from one experiment, or its measured matrix
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
post-neutron mass (`SEC`) and is not used for pre-neutron data. The joint yields from the thesis
inputs remain beside them and can still be named by path.

`U235_nres` is the same target as `U235_nth` on the resonance channel. Its `R_T(A_H)` trend and
its yield inputs are the thermal ones and are read from the `U235_nth` tree rather than copied;
only the measurements made on the resonance beam are its own, `ν(A)`, `ν(A, TKE)` and `⟨ε⟩(A)`.

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

Each retrieval writes into the checkout of that package unless an output root is given as the
second argument; the root of this repository, `julia scripts/retrieve.jl <configuration>
<this repository>`, puts the tables under `data/` here. The yield inputs of the shipped
configurations:

```bash
git clone --branch v0.5.0 https://github.com/PaulGoG/ExforFissionData.jl
cd ExforFissionData.jl
julia scripts/retrieve.jl config/Cf252_sf_Y_vs_A_TKE.toml
julia scripts/retrieve.jl config/U235_nth_Y_vs_A.toml
julia scripts/retrieve.jl config/U235_nth_TKE_vs_A.toml
julia scripts/retrieve.jl config/U235_nth_sigma_TKE_vs_A.toml
julia scripts/retrieve.jl config/Pu239_nth_Y_vs_A.toml
julia scripts/retrieve.jl config/Pu239_nth_TKE_vs_A.toml
julia scripts/retrieve.jl config/Pu240_sf_sigma_TKE_vs_A.toml
julia scripts/retrieve.jl config/U233_nth_Y_vs_A.toml
julia scripts/retrieve.jl config/U233_nth_TKE_vs_A.toml
julia scripts/retrieve.jl config/U233_nth_sigma_TKE_vs_A.toml
```

The neutron measurements the figures compare against:

```bash
julia scripts/retrieve.jl config/Cf252_sf_nu_vs_A.toml
julia scripts/retrieve.jl config/U235_nth_nu_vs_A.toml
julia scripts/retrieve.jl config/Pu239_nth_nu_vs_A.toml
julia scripts/retrieve.jl config/U233_nth_nu_vs_A.toml
julia scripts/retrieve.jl config/Cf252_sf_eps_vs_A.toml
julia scripts/retrieve.jl config/U235_nth_eps_vs_A.toml
julia scripts/retrieve.jl config/Pu239_nth_eps_vs_A.toml
julia scripts/retrieve.jl config/U233_nth_eps_vs_A.toml
julia scripts/retrieve.jl config/Cf252_sf_spectrum_maxwellian_ratio_vs_E.toml
julia scripts/retrieve.jl config/U235_nth_spectrum_cf252_ratio_vs_E.toml
julia scripts/retrieve.jl config/Pu239_nth_spectrum_cf252_ratio_vs_E.toml
julia scripts/retrieve.jl config/U233_nth_spectrum_cf252_ratio_vs_E.toml
```

The tree is the ExforFissionData.jl v0.5.0 retrieval of every configuration that package ships,
the resonance channel of ²³⁵U and ²⁴⁰Pu(sf) included. Beyond the pre-neutron observables, it
holds the three post-neutron directories and the prompt-γ directories described below, and its
records state the unit of every subentry, the sum of a yield table and the handling of isomeric
states. From the ²⁵²Cf spectrum of Blinov (1980), 40535004, a ratio to a Maxwellian, the line at
220.9 keV is left out, an archive defect its record states. Four of the configurations, the joint
`Y(A,TKE)` of ²³⁵U, ²³⁹Pu and ²³³U and the product kinetic energies of ²⁵²Cf, find no dataset in
the archive.

**Evaluated product yields.** `<system>/Y_vs_A_p/` and `<system>/Y_vs_Z_A_p/` hold the
independent fission product yields of ENDF/B-VIII.0 (Brown *et al.*, *Nucl. Data Sheets* **148**,
1 (2018), doi:10.1016/j.nds.2018.02.001) and JEFF-3.3 (Plompen *et al.*, *Eur. Phys. J. A* **56**,
181 (2020), doi:10.1140/epja/s10050-020-00141-9), which the post-neutron product yields of the
model are compared against. Each directory has one file per library, `ENDF-B-VIII.0.dat` and
`JEFF-3.3.dat`, in per cent per fission, 200 over both products: `A_p Y Y_uncertainty` against
product mass, and `Z A_p Y Y_uncertainty` resolved in charge. They are taken from MF8 MT454 at
0.0253 eV, or at the single point of ²⁵²Cf(sf). Isomeric states are summed and uncertainties
added in quadrature, the files stating no correlations; the ternary particles JEFF-3.3 lists
among its products are left out. They are not measurements, and the figures draw them as lines.

```bash
julia scripts/stage_product_yields.jl
```

downloads the NFY and SFY sublibraries from the NNDC and from the NEA Data Bank
(doi:10.82555/nhfqy-hvz09) and writes both directories for the four systems; a `provenance.toml`
beside the tables records the URL, the SHA-256 of what was downloaded, the material and the energy
point. `--library` and `--system` restrict the staging; if a server refuses a scripted download,
fetch the file in a browser and pass it with `--file <library>=<path>`. A table that differs from
the one in place never replaces it; the old one moves aside as `<library>#1.dat`.

**Measured post-neutron products.** ExforFissionData.jl retrieves three post-neutron
observables, which the figures draw as markers against the model and the evaluations:
`Y_vs_A_p/`, the product mass yields EXFOR codes `MASS,SEC,FY` and the chain yields it codes
`MASS,CHN,FY`, which β decay, conserving `A`, makes post-neutron mass yields up to the
delayed-neutron transfer between neighbouring chains; `Y_vs_Z_A_p/`, the independent yields
`ELEM/MASS,IND,FY`; and `E_K_p_vs_A_p/`, the mean kinetic energy of the products against their
mass, `MASS,SEC,KE`, of which the archive holds none for ²⁵²Cf. Each holds
`<accession>_<Author>_<year>.dat` files with the `retrieval.toml` and the `subentries/` of the
retrieval. The yields are written in fragments per fission, 2 over both products, against the
per cent of the evaluations; `unit_written` in the record states the unit of each table, and
`plotting/figures.jl` converts from it and stops on a unit it does not know. The kinetic energies
are in MeV.

```bash
git clone --branch v0.5.0 https://github.com/PaulGoG/ExforFissionData.jl
cd ExforFissionData.jl
julia scripts/retrieve.jl config/<system>_Y_vs_A_p.toml
julia scripts/retrieve.jl config/<system>_Y_vs_Z_A_p.toml
julia scripts/retrieve.jl config/<system>_E_K_p_vs_A_p.toml
```

for each of `Cf252_sf`, `U235_nth`, `Pu239_nth` and `U233_nth`. In `Y_vs_A_p/` and
`Y_vs_Z_A_p/` the retrieval sits beside the evaluations: the measured tables, `retrieval.toml`
and `subentries/` next to `<library>.dat` and `provenance.toml`, with no name in common. A
retrieval never writes into a directory that exists, and the staging script writes only its own
files, so the two are merged by hand: the contents of the directory the retrieval writes are
copied into the one holding the evaluations.

**Prompt γ rays.** ExforFissionData.jl v0.5.0 retrieves the prompt γ observables: per fission,
tabulated against nothing, the multiplicity `M_gamma_bar/`, the total energy `E_gamma_tot/` and the
mean photon energy `eps_gamma/` of the four systems, one row to a file; for ²⁵²Cf also the energy
of one fragment against its mass, `E_gamma_vs_A/`, the total energy against fragment mass,
`E_gamma_tot_vs_A/`, against TKE, `E_gamma_tot_vs_TKE/`, and against the number of neutrons,
`E_gamma_tot_vs_nu/`, and the multiplicities `M_gamma_vs_A/` (one fragment), `M_gamma_bar_vs_A/` and
`M_gamma_bar_vs_TKE/` (per fission); for ²³⁵U `M_gamma_vs_A/`, `M_gamma_bar_vs_A/` and
`eps_gamma_vs_A/`. A measured γ energy is defined by its γ-energy range, its time window and its
corrections, which differ from set to set by more than the stated uncertainties; each record states
them as the entry or its publication does, flags a set coded `PAR` (`gamma_range_partial`), and its
`[datasets]` warnings name the sets whose definition departs. A range-limited set is not a total.

```bash
julia scripts/retrieve.jl config/<system>_M_gamma_bar.toml
julia scripts/retrieve.jl config/<system>_E_gamma_tot.toml
julia scripts/retrieve.jl config/<system>_eps_gamma.toml
```

for each of the four systems, and the configurations named for the directories above for ²⁵²Cf
and ²³⁵U.

Read the `retrieval.toml` it writes before using what it retrieved. It records, per dataset, the
reaction code, the units and the qualifiers the archive attached, and those qualifiers matter:
`DERIV` marks a dataset derived from other data rather than measured, `SPA` one averaged over a
fission spectrum rather than the thermal one it is filed under. No yield input carries either, and
the neutron figures draw no such set; the product figures draw them with open markers, and
`results/README.md` lists them. Of the post-neutron sets, those the record marks superseded are
not drawn, and neither is the ²³⁵U product mass yield of Koblik (2014), 31738002, whose heavy
products alone sum to 182 % per fission. An uncertainty the archive does not state is written
`NaN`, never 0; a 0 appears only where EXFOR states it, and `uncertainty_zero_rows` counts those
rows.

The neutron-energy and spectrum-ratio tables need their record more than the others do:
- `eps_vs_A`, `eps_vs_A_TKE`: each dataset states `ordinate_frame`, the centre of mass for every
  set held, and `ordinate_frame_basis`, whether that reading rests on the heading of the
  subentry, on its text or on the publication. `mean_formed_from` says how the mean was formed:
  the first moment of the measured spectrum, of a form fitted to it, of the spectrum completed
  by such a form, or unstated. `mean_threshold_mev` and `mean_threshold_frame` give the
  low-energy threshold of the measurement, 0.2 to 0.7 MeV where one is stated; a calculated mean
  compares with a measured one under the same cut. The sets of Batenkov (2004) carry a
  `preliminary:` qualifier, as their publication calls its results. `MXW`, on the sets of
  Nishio (1998) and Tsuchiya (2000), is the average over the incident thermal spectrum, as on
  their ν(A) sets.
- `T_eff_vs_A`, `lambda_vs_A`: `fit_form` gives the form each parameter belongs to and
  `fit_evidence` the equation of the publication.
- `P_vs_nu`: a distribution is written as tabulated and never renormalised; `distribution_sum`,
  `mean_multiplicity` and its deviation from `ν̄`, `mean_deviation`, are recorded per dataset.
  The Boldeman (1967) sets of ²³³U, ²³⁵U and ²³⁹Pu are superseded by the 1985 reanalysis of the
  same measurement, and the ²⁵²Cf sets 30046011 and 30772015 are four runs of one experiment at
  different discriminator biases; each names the others in `correlated_with`. The evaluations
  of Holden (1988), V0045011, V0045012 and V0045013 for ²³³U, ²³⁵U and ²³⁹Pu, are written like a
  measurement and flagged `evaluation:`, with `evaluation_warning` naming them under
  `[datasets]`; the archive holds none for ²⁵²Cf.
- `spectrum_cf252_ratio_vs_E`: both orientations occur in one directory and none is inverted;
  `ratio_orientation` says which a dataset is. The ²³⁵U sets 41516017 and 41597002 are one
  measurement published twice and name each other in `correlated_with`; the first is superseded
  by the second.
- `spectrum_maxwellian_ratio_vs_E`: `maxwellian_temperature_mev` is the temperature of the
  reference Maxwellian, read by the retrieval from the subentry. The Starostov (1985) sets of
  ²³³U, ²³⁵U and ²³⁹Pu are four parts of one spectrum each; for ²⁵²Cf, 40418006 and 40875003
  are superseded by 40418008 and 41158003.

Three statements of the archive about the standing of a dataset are carried into its record,
and none removes it:
- `preliminary:` among the `qualifiers` marks a set whose subentry, or the common subentry of
  its entry, carries the status `PRELM`. Of the sets the figures draw, the ²³⁵U `ν(A)` of
  Vorobyev (2010), 41516012, is one.
- `superseded:` names the later subentry that replaces a set (status `SPSDD`). Both are
  written, and the later one is not always the fuller table: Lajtai (1990), 41158003, quotes no
  uncertainties where the superseded 40875003 does.
- `correlated_with` lists the sets that share an experiment with a set, and
  `correlation_relation` says how: `republication`, `alternative_analysis`, `repeated_run`,
  `complementary_range` or `dependent`, where one member is formed from another. Any
  combination of them is one measurement. The `ν(A)` of Basova (1979) and Zamyatnin (1979), for
  ²⁵²Cf and for ²³⁹Pu, are two reductions of one experiment.

Four further qualifiers say what a table holds where its directory alone would mislead:
- `multiplicity_from_masses:` on the ²³⁹Pu `ν(A)` of Nishio (1995), 23012008, which the figure
  draws: the multiplicity is the difference of the pre- and post-neutron fragment masses of one
  measurement, not a count of detected neutrons, and the paper puts the total it implies about
  10 % above the evaluated one.
- `mean_threshold_unsettled:` on the ²³⁹Pu `⟨ε⟩(A)` of Tsuchiya (2000), 22650008: the
  publication averages the points above 0.5 MeV and does not say whether the range below
  enters the mean, so `mean_threshold_mev` is not to be applied as a cut.
- `spectrum_form:` on four ²⁵²Cf sets under `spectrum_vs_E` that tabulate `N(E)/√E`.
- `evaluation:` on the ²⁵²Cf spectrum of Mannhart (1987), V0101002, and on the P(ν) of Holden
  (1988), V0045011, V0045012 and V0045013, which are evaluations and not measurements.

Eight ν(A) subentries are coded without the `/FRG` that marks a multiplicity per fragment: for ²⁵²Cf
Britt (1964), Piksaykin (1977), Budtz-Jørgensen (1988), Zeynalov (2011) and Göök (2014); for ²³⁹Pu
Tsuchiya (2000) and Batenkov (2004); for ²³⁵U Batenkov (2004). So is the ν(A,TKE) of
Budtz-Jørgensen. The retrieval reads them by a complement test on the data. For a per-fragment
`ν(A)` the pair sum `ν(A) + ν(A₀−A)` varies along the sawtooth, and only its yield-weighted mean
equals `ν̄`; a multiplicity per fission is the same for both complements. The record of each dataset
gives the basis of the reading, the deviation of the weighted pair sum from `ν̄`, and the paper that
corroborates it where one was consulted. For two of them the pair sum does not return `ν̄`: the
²³⁹Pu sets of Tsuchiya, +3.9 ± 0.6 %, which its publication normalises to 2.88 neutrons per fission,
and Batenkov, −10.7 ± 1.3 %. They are written uncorrected, flagged `scale_consistent = false`, and
not drawn in the figures. Both ν(A) sets of Batenkov also carry the `preliminary:` qualifier.

**The temperature ratio `R_T(A_H)`** for the temperature-ratio TXE partition —
[`FissionTemperatureRatio.jl`](https://github.com/PaulGoG/FissionTemperatureRatio.jl) extracts it
from experimental `ν(A)` and tabulates it at every integer `A_H`.

**Stage the whole results directory of one run, not one table out of it.** `scripts/run.jl`
of that package writes `data/sims/<system>/<run>/`, some seventy-five files for ²⁵²Cf; the ones this
package reads are

```
data/sims/<system>/<run>/
├── manifest_<run>.toml                      the run record; this package reads it first
├── R_T_vs_A_H_segmented_<curve>.csv         one per [[segmented_curve]]; A_H,R_T,R_T_uncertainty
├── r_nu_vs_A_H_pivots_<curve>.csv           its pivots; A_H,r_nu,r_nu_uncertainty; not read here
├── total_average_R_T_<run>.csv              ⟨R_T⟩ per curve and yield distribution
├── dataset_diagnostics.csv
├── leave_one_out_<run>.csv                  the refits with each dataset left out; not read here
├── metadata.toml                            commit, platform, every accession used
└── configuration.toml, Manifest.toml        the run's own configuration and environment
```

beside the point-by-point and consensus tables of every dataset. `<curve>` is
`<accession>_<Author>_<year>` for a dataset and `systematic_trend` for the trend; the manifest
also records the accession of each dataset curve. Copy the whole directory to
`data/<system>/temperature_ratio/` and leave `[partition] subdirectory` pointing there;
`[partition] curve_label` then picks one of the manifest's `[[segmented_curve]]` entries by its
`label`, `"systematic_trend"` for the shipped configurations. Keep one staged run per system —
two manifests in one directory is refused, because nothing would say which run to take.

The staged runs are those of the configurations FissionTemperatureRatio.jl v0.4.0 ships, on the
ExforFissionData.jl v0.4.0 retrieval, with each system's primary `Y(A)` and `⟨TKE⟩(A)` as above;
the `ν(A)`, `Y(A)` and `⟨TKE⟩(A)` tables of that retrieval are identical to those of v0.5.0.
`U235_nres` reads the thermal run. Each staged directory carries the configuration and the
metadata of its run. To produce the ²⁵²Cf(sf) run again:

```bash
git clone --branch v0.4.0 https://github.com/PaulGoG/ExforFissionData.jl
git clone --branch v0.4.0 https://github.com/PaulGoG/FissionTemperatureRatio.jl
cd ExforFissionData.jl
julia scripts/retrieve.jl config/Cf252_sf_nu_vs_A.toml ../FissionTemperatureRatio.jl
julia scripts/retrieve.jl config/Cf252_sf_Y_vs_A.toml ../FissionTemperatureRatio.jl
julia scripts/retrieve.jl config/Cf252_sf_TKE_vs_A.toml ../FissionTemperatureRatio.jl
cd ../FissionTemperatureRatio.jl
julia scripts/run.jl config/Cf252_sf.toml
```

and copy the directory it writes, `data/sims/Cf252_sf/<run>/`, whole to
`data/Cf252_sf/temperature_ratio/` here. `U235_nth`, `Pu239_nth` and `U233_nth` follow the same
three retrievals and their own configuration.

A table copied by hand cannot be checked: the segmented `R_T` curve and the pivots file of the
tree above hold different quantities in the same column position, `R_T ≈ 1.2` against
`r_ν ≈ 0.5`, and both are plausible. Staging the pivots file by mistake runs to completion and
reports a number that looks like a result. The manifest is what makes the quantity checkable —
this package refuses a manifest that does not declare `ordinate = "R_T"` and
`abscissa = ["A_H"]`, and resolves every curve through `temperature_ratio_file` and through no
other field. The `columns` the manifest states are documentation: tables are read by position.

Run the extraction with `ratio_averaging = "charge_resolved"` and
`zero_polarization_at_symmetry = false`, as the configurations of that package do; the reason for
the first is in `read_segmented_curve`'s docstring, and the manifest's `[domain]` record, the
second included, is checked on load.

The remaining tables:

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
doi:10.1016/j.phpro.2015.04.009, applying the `Zₚ` model of Wahl, *At. Data Nucl. Data Tables*
**39**, 1 (1988), doi:10.1016/0092-640X(88)90016-2. They are not regenerated here. The first three
are thesis inputs, recoverable as below; the ²³³U table is not, and is not distributed.
`charge_distribution_file = "wahl"` takes Wahl's model directly instead. FissionFragmentsDomain.jl
reproduces the 1988 evaluation's own tables row by row, and records in its `data/README.md` how
far these supplied tables lie from it.

**`reference/ground_state_deformation_moller.dat`** — ground-state quadrupole deformations from
the finite-range liquid-drop model. P. Möller, J. R. Nix, W. D. Myers, W. J. Swiatecki, *Atomic
Data and Nuclear Data Tables* **59**, 185 (1995), via RIPL.

**`U235_nth/E_def_vs_A_Z.dat`** — not a source table: the extra deformation energy at scission
`ΔE_def`, generated from the ground-state deformations above by
`Data conversion/generateDeltaEdeformation.jl` on the `original` branch. It is the input the
scission-modelling partition needs.

**The fragment mass yields.** `U235_nth/Y_vs_A_H_TKE.dat` is the ²³⁵U(n_th,f) measurement of
C. Straede *et al.*, *Nucl. Phys. A* **462**, 85 (1987), doi:10.1016/0375-9474(87)90381-2, as its
header records, and `Pu239_nth/Y_vs_A_TKE_chiba.dat` a ²³⁹Pu(n_th,f) set its header labels Chiba,
25 July 2023. `Cf252_sf/Y_vs_A_TKE.dat`, `U235_nth/Y_vs_A_TKE.dat` and `Pu239_nth/Y_vs_A_TKE.dat`
are the joint distributions on the whole `A` domain with an uncertainty column;
`Cf252_sf/Y_vs_A_TKE.dat` is the ²⁵²Cf(sf) distribution behind every result of chapter 4 of the
thesis. No shipped configuration reads any of them, every system taking its yield from EXFOR as
tabulated above, but a configuration can name one by path.

**The thesis inputs.** These tables, the charge distributions of ²⁵²Cf, ²³⁵U and ²³⁹Pu, the joint
yields, `E_def_vs_A_Z.dat`, and the deformation, shell-correction and resonance tables under
`reference/`, are inputs of my 2023 MSc thesis, whose code and text are on the `original` branch.
`input_data/README.md` there gives the origin of each, and the import commit `bdba1d9` holds them
under the names the thesis code reads:

| Here | `input_data/` at `bdba1d9` |
|---|---|
| `Cf252_sf/Y_vs_A_TKE.dat` | `Cf-252/CF52YATKE.VES` |
| `U235_nth/Y_vs_A_TKE.dat` | `U-235/U5YATKE.ANA` |
| `U235_nth/Y_vs_A_H_TKE.dat` | `U-235/U5YATKE.SRE` |
| `Pu239_nth/Y_vs_A_TKE.dat` | `Pu-239/P9YATKE.ANA` |
| `Pu239_nth/Y_vs_A_TKE_chiba.dat` | `Pu-239/P9YATKE1.CHI` |
| `U235_nth/E_def_vs_A_Z.dat` | `U-235/EXTRADEF.U5` |
| `Cf252_sf/`, `U235_nth/`, `Pu239_nth/charge_distribution_vs_A.dat` | `Cf-252/DeltaZ_rms_A.CF52`, `U-235/DeltaZ_rms_A.U5`, `Pu-239/DeltaZ_rms_A.PU39` |
| `reference/ground_state_deformation_moller.dat` | `Aux files/B2MOLLER.ANA` |
| `reference/shell_corrections_gilbert_cameron.dat` | `Aux files/SZSN.GC` |
| `reference/s_wave_resonances_ripl3.dat` | `Aux files/resonances0.dat` |

```bash
git show "bdba1d9:input_data/Cf-252/CF52YATKE.VES" > data/Cf252_sf/Y_vs_A_TKE.dat
```

restores one. Every data row here is byte-identical to that commit's; only the header lines were
rewritten to the names of this package, and readers take columns by position.

## Citing

Experimental data belongs to the groups that measured it. Any published result derived from the
yields or multiplicities here should cite the measurement, not this package — the EXFOR accession
number that opens each retrieved file name is what makes those citations recoverable, and the
`retrieval.toml` beside it records the reaction code and units.
