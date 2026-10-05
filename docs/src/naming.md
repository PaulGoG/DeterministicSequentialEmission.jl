```@meta
CurrentModule = DeterministicSequentialEmission
```

# Naming

A quantity has exactly one name, and that name is used everywhere it appears: Julia identifier,
struct field, configuration key, directory name, file name, column header, figure label, prose.
The spelling changes register between contexts — `charges_per_mass` in a configuration, `Z` in a
header — but the *quantity* does not.

The corollary is the rule that decides the hard cases: **a name says what the thing is, not what
was convenient to type.** `Y_vs_A_TKE` is a statement; `Y_A_TKE` is a list of symbols with no
relation between them. `Data` as a suffix carries nothing. `err` names a different concept from
the one it was used for — an error is not an uncertainty. And a name that is relative to a moment
ages: a table with no header line is `HEADERLESS_MASS_EXCESS_SPEC`, not a legacy one.

This page is the vocabulary as this package applies it. The retrieval that supplies its
experimental input and the extraction that supplies its temperature ratio share the quantity table
and the file layout, so a name learned here is the name there.

## Two registers

| Register | Spelling | Where |
| :--- | :--- | :--- |
| word | `heavy_mass_max`, `temperature_ratio` | configuration keys and values, prose |
| symbol | `A_H`, `R_T`, `nu` | directory names, file names, column headers, figure labels, Julia identifiers of the equations |

A configuration is edited by hand and has to explain itself, so it spells a quantity out. A path
and a header are read at a glance and are the field's own nomenclature, so they carry the symbol.
This is why the configuration key is `charge_distribution_file` while the file it names is
`charge_distribution_vs_A.dat` with the header `A dZ sigma_Z`, and why
[`ChargeDistribution`](@ref) has the fields `ΔZ` and `σ_Z`, mirroring the equations.

## The quantities

| Symbol | ASCII token | Julia | Quantity |
| :--- | :--- | :--- | :--- |
| *A* | `A` | `A` | pre-neutron fragment mass |
| *A*_H, *A*_L | `A_H`, `A_L` | `A_H`, `A_L` | heavy, light fragment mass |
| *A*_p | `A_p` | `A_p` | post-neutron product mass |
| *A*₀, *Z*₀ | `A_0`, `Z_0` | `A₀`, `Z₀` | compound nucleus mass, charge |
| *Z* | `Z` | `Z` | fragment charge |
| *N* | `N` | `N` | neutron number |
| *Z*_p | `Z_p` | `Zₚ` | most probable charge |
| Δ*Z* | `dZ` | `ΔZ` | charge polarization |
| σ_*Z* | `sigma_Z` | `σ_Z` | charge dispersion |
| *p*(*Z*,*A*) | `p_Z` | `p` | isobaric charge probability |
| TKE | `TKE` | `TKE` | total kinetic energy, pre-neutron |
| TXE | `TXE` | `TXE` | total excitation energy |
| *E** | `E_star` | `excitation` | fragment excitation energy |
| *Q* | `Q` | `Q` | Q value |
| Δ | `mass_excess` | `Δ` | mass excess |
| *S*_n | `S_n` | `Sₙ` | neutron separation energy |
| *S*(*Z*), *S*(*N*) | `S_Z`, `S_N` | `S_Z`, `S_N` | shell correction |
| Δ*E*_def | `E_def` | `ΔE_def` | extra deformation energy at scission |
| *Ē*_r | `E_r` | `residual_energy` | residual energy before an emission |
| *T* | `T` | `T` | residual temperature |
| ⟨ε⟩ | `eps` | `ε` | mean centre-of-mass neutron energy |
| *E* | `E` | `E` | laboratory neutron energy |
| *a* | `a` | `a` | level density parameter |
| ν | `nu` | `ν` | prompt multiplicity per fragment |
| *Y* | `Y` | `Y` | fragment yield |
| *N*(*E*) | `spectrum` | `spectrum` | prompt neutron spectrum |
| *R*_T | `R_T` | `R_T` | temperature ratio, T_L/T_H |

`rms` is retired. It named a root-mean-square of nothing stated; the quantity is the Gaussian
dispersion of the isobaric charge distribution, and the field writes σ_Z.

A yield type is named by its **index space**, which is the distinction that matters when more than
one is in scope: [`MassEnergyYield`](@ref) is the measured *Y*(*A*, TKE) and
[`FragmentYield`](@ref) the model weight *Y*(*A*, *Z*, TKE). Neither `y` nor `Y` may stand for the
other: one quantity, one spelling, and here it is `Y`.

### Uncertainty

One spelling, one position: **`<quantity>_uncertainty`**, immediately after the quantity it
belongs to. `A nu nu_uncertainty`, never `A nu errnu`, and never a bare `uncertainty` column that
leaves the reader to work out which quantity it belongs to. In Julia the uncertainty travels with
the value inside a `Measurement` instead, so it needs no name of its own.

## Julia identifiers

| Kind | Case | Example |
| :--- | :--- | :--- |
| module, type, abstract type | `CamelCase` | [`SegmentedCurve`](@ref), [`LevelDensityModel`](@ref) |
| function, macro | `lowercase_snake_case` | [`level_density_parameter`](@ref), [`system_label`](@ref) |
| variable, field, keyword | `lowercase_snake_case` | `charges_per_mass`, `significant_digits` |
| constant | `SCREAMING_SNAKE_CASE` | [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref) |

Unicode is used where an identifier *is* a symbol of the governing equations — `ΔZ`, `σ_Z`, `Zₚ`,
`Sₙ`, `ε`, `Δ` — and only there, so the implementation reads like the formulation. The exported
API stays ASCII-typeable: every public function can be called, and every keyword argument of one
given, without a compose key. That is why [`read_charge_distribution`](@ref) takes
`default_charge_dispersion` while the field it fills is `default_σ_Z`.

Verb prefixes are a small closed set — `read_`, `write_`, `build_`, `load_`, `run_`, `is_` — and a
function that fits none of them is named for what it returns. `read_<thing>` must return the type
named `<Thing>`, which is the rule that decides [`read_mass_excess_table`](@ref) against
`read_mass_excess`: the type is [`MassExcessTable`](@ref), so the reader carries `table` too, and
the same rule turns `read_segments` into [`read_segmented_curve`](@ref). `load_<thing>` reads
*and* validates, which is why the configuration entry point is [`load_configuration`](@ref) and
not `read_configuration`.

An accessor is named for the quantity it returns, not for the access:
[`charge_polarization`](@ref)`(distribution, A)`, not `polarization_at(distribution, A)`.

Type names carry no `Data` suffix and no abstract noun standing in for a description: the
piecewise-linear curve a partition is read from is a [`SegmentedCurve`](@ref), not a
`JointedSegments` — "segmented" says how, and a plural of a part is not the name of the whole. A
dropped input is an [`Exclusion`](@ref). `Table` survives only where the thing genuinely is a
lookup table keyed by nucleon number: [`MassExcessTable`](@ref),
[`ShellCorrectionTable`](@ref).

One word may not carry two meanings in one module. [`Configuration`](@ref) is the validated TOML
and is shared with the codes upstream; a quantity tabulated over fission configurations
(*A*, *Z*, TKE) is a [`FragmentQuantity`](@ref). "Fission configuration" survives in prose for the
triple, where it is unambiguous.

The result type of a package is `<Verb>Result`, **one per package**. Here it is
[`RunResult`](@ref): this package runs a model. The codes upstream name theirs for what they do,
so `using` any pair of them collides on nothing.

Abbreviations are permitted only where the field itself uses them — `TKE`, `TXE`, `Q`, `BSFG`,
`GC`, `EXFOR`, `AME`, `RIPL`, `sf` — and never as `err`, `param`, `calc`, `val`, `idx`, `num`,
`cfg`, `pts`. The one place short forms are unavoidable is a run identifier, and they are
collected there in a single table; see below.

## Configuration files

A configuration is named `<system>.toml`, the token that names the data directories it reads and
the output directory it writes. A system carrying more than one configuration would name them
`<system>_<case>.toml`; none does. The test suite asserts the name against the system each file
declares, so a file cannot drift from what it runs.

Sections are `[system]`, `[fragmentation]`, `[energetics]`, `[level_density]`, `[evaporation]`,
`[partition]`, `[yield]`, `[spectrum]` and `[output]`, and **a key is never prefixed with its
section's name**: `[level_density] model`, not `level_density_model`.

A key's comment states three things and nothing else: what the key is, in one short line; the
enumerated choices where it is constrained; and the bounds and the unit. No physics narrative, no
run rationale — those belong in the documentation. The loader enforces exactly what the comment
declares, so a run cannot start from a configuration it cannot honour, and it names the offending
key when it refuses.

**A key that can only be redundant or wrong is not in the file.** The reaction code follows from
the entrance channel through [`CHANNEL_REACTION`](@ref), so it is not written down; the compound
nucleus follows from the target and the channel through [`REACTIONS`](@ref), so it is not written
down either. Nor is the output directory: it is [`system_label`](@ref) of the system the file
already declares. It existed as `[output] subdirectory` only because two runs of one target could
not otherwise be told apart, and the channel is what tells them apart now.

`subdirectory` names a folder under a root; `directory` only ever names a root path. So every
input key in a configuration is a path under the data root, and the root itself is an argument to
[`load_configuration`](@ref).

`significant_digits` is significant figures, and has always been so here. It is not decimal
places, and the distinction is not cosmetic: a multiplicity near 3 and a yield near `10⁻⁴⁰` appear
in the same run, and a fixed number of decimals cannot carry both to the same meaning.

## The system identifier

One token for a fissioning system, used as a directory name, a configuration file name and an
output directory:

```
<ElementSymbol><A>_<entrance channel>
```

with the channel spelled as the field spells it and not as a reaction code writes it — `Cf252_sf`,
`U235_nth`, `U235_nres`, `U233_nth`, `Pu239_nth`. `0f` is a reaction code, not a name, and it
belongs in the run record where the reaction code already is. The channel is also what
distinguishes a thermal from a resonance run of one target, which `n,f` alone cannot.

[`system_label`](@ref) returns that token. The typeset form for a figure or a caption —
`²⁵²Cf(sf)`, `²³³U(nth,f)` — is [`system_notation`](@ref), a separate function: one name may not
mean both.

## Directories, files and headers

```
data/
├── reference/                      evaluations that belong to no one system
│   ├── shell_corrections_gilbert_cameron.dat
│   ├── s_wave_resonances_ripl3.dat
│   └── ground_state_deformation_moller.dat
└── <system>/                       Cf252_sf, U235_nth, U235_nres, U233_nth, Pu239_nth
    ├── charge_distribution_vs_A.dat
    ├── E_def_vs_A_Z.dat
    ├── temperature_ratio/          one staged FissionTemperatureRatio.jl run
    ├── Y_vs_A_TKE.dat
    └── nu_vs_A/
        └── <accession>_<Author>_<year>.dat
```

One directory per system, one file per quantity, and a subdirectory where a quantity was measured
more than once — there the **directory** carries the quantity and the **file** carries the
provenance, which is what makes the archive accession number recoverable from a figure legend.
There is no second axis: a system is one token and one tree, so the yields of a system and its
multiplicities are not kept apart under different roots.

A tabular file is named `<quantity>_vs_<abscissa>[_<abscissa>]`. `vs` is what makes the name a
statement rather than a list of symbols, and it is precisely what `Y_A_TKE` is missing. Where a
system holds two tables of one quantity in one index space, the authoritative one takes the plain
name and the other carries a qualifier, exactly as a generated table does.

Extensions state a format and nothing else: `.dat` whitespace-separated, `.csv` comma-separated,
`.toml` metadata. **A file extension never carries a nuclide**, because the nuclide is already in
the path.

A generated table is marked as generated, in the name and in the header row —
`charge_distribution_vs_A_wahl_systematics.dat` — so it can never be mistaken for an evaluation.
Where several generations of one quantity sit beside each other the qualifier names the model that
produced it and where its parameters came from, not its rank: `_wahl_per_reaction`,
`_wahl_systematics`, `_mean`.

One header line, ASCII, abscissae first, then the ordinate, then its uncertainty:

```
A dZ sigma_Z
A TKE Y Y_uncertainty
A nu nu_uncertainty
A_H,R_T,R_T_uncertainty
```

**No column is named `value`, and no uncertainty column is named `uncertainty`.** Naming a column
for its quantity is what frees the reader from having to know the file name to know what the file
is holding, and it matters most for the files that cross a package boundary. This binds a
`TableSpec`'s own column symbols as much as a header on disk: a spec reading `[:A_H, :value]`
states nothing, which is why there are two curve specs — [`TEMPERATURE_RATIO_CURVE_SPEC`](@ref)
and [`EXCITATION_RATIO_CURVE_SPEC`](@ref) — rather than one that could not name its ordinate.

**A quantity that crosses a package boundary is taken through the producer's run record, not
through a path.** `R_T(A_H)` and `r_ν(A_H)` are written side by side, in the same column of two
files, and are both plausible floats; a configuration naming one of them by hand can select the
wrong physics and finish without complaint. [`TemperatureRatioManifest`](@ref) is what makes the
quantity checkable, and [`temperature_ratio_path`](@ref) is the only route from it to a file.

Readers take columns **by position**, not by header text. That is what makes a header rename a
no-op for code, upstream and downstream alike, and it must stay that way: every reader in this
package says so in its docstring, and nothing may be changed to a lookup by name. The one
evaluation held as distributed, the s-wave resonance library, keeps its own commented header
untouched, and [`HEADERLESS_MASS_EXCESS_SPEC`](@ref) reads an atomic mass table with no header line
at all — the same rule seen from the other side.

## Output files and run identifiers

A run writes to `data/sims/<system>/<run>`, and an output follows the same
`<quantity>_vs_<abscissa>` rule as an input, so it can be fed back in without translation:

| File | Content |
| :--- | :--- |
| `nu_vs_A.dat` | prompt multiplicity against fragment mass |
| `T_vs_A.dat` | residual temperature |
| `E_star_vs_A.dat` | fragment excitation energy |
| `E_r_vs_A.dat` | residual energy before emission |
| `eps_vs_A.dat` | mean centre-of-mass neutron energy, per emitted neutron under the default sequence weighting |
| `nu_vs_TKE.dat`, `E_star_vs_TKE.dat` | the same against total kinetic energy |
| `spectrum_vs_E.dat` | the laboratory prompt neutron spectrum |
| `Y_vs_A.dat`, `Y_vs_Z.dat`, `Y_vs_TKE.dat` | the fragment yield projections |
| `totals.dat` | one row per quantity averaged over the whole run |
| `metadata.toml` | how it was produced: configuration, commit, environment, machine |
| `configuration.toml` | the configuration the run loaded, copied |
| `Manifest.toml` | the resolved environment the run executed in, copied |

**A run-identifier token is the configuration key it came from.** Spelling every key out in full
would put some two hundred characters into every directory name, so the keys are abbreviated — and
abbreviated in one place, [`RUN_IDENTIFIER_ABBREVIATIONS`](@ref), rather than at the point of use.
[`run_identifier`](@ref) refuses a key with no entry rather than inventing a token, and the test
suite asserts that every token it uses has one and that no two keys share a token.

The identifier says which keys it covers rather than claiming all of them. What it leaves out — the
input file paths, the spectrum grid, `significant_digits` — is written in full into the metadata
beside the result, and [`run_parameters`](@ref) records the reason for each.
