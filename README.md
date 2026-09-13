# DeterministicSequentialEmission.jl — `original`

This branch is the MSc thesis code exactly as submitted to the Faculty of Physics,
University of Bucharest in September 2023, together with the thesis it belongs to. It is kept
as a reference point. `main` is a rewrite.

Nothing on this branch has been repaired. There is no project environment, the files are
top-level scripts communicating through globals rather than a module, and the defects listed
below are left in place deliberately.

## What it computes

Prompt neutron emission from fission fragments, treated sequence by sequence. Where global
treatments integrate an evaporation spectrum over a residual temperature distribution and
return averages, the DSE model solves the energy conservation equation of every emission
sequence separately:

```
Ē_r^(k-1) − S_n^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … n(A,Z,TKE)
```

with Weisskopf-Ewing evaporation spectra in the centre-of-mass system and fragment level
densities in the Fermi-gas regime. Because the number of sequences `n` is an integer at every
configuration, the post-neutron fragment mass `A_p = A − n` is an integer too, without rounding.

The fragmentation domain is built deterministically: symmetric split to maximum probable
asymmetry in steps of one mass unit, an odd number of charge numbers per `A` centred on
`Z_p(A) = Z_UCD(A) + ΔZ(A)` with Wahl polarization, each weighted by the Gaussian isobaric
charge distribution `p(Z,A)`, crossed with the full TKE grid. Results are matrices
`q_k(A,Z,TKE)`, optionally averaged over sequences and over an experimental `Y(A,TKE)` to give
`q(A)`, `q(Z)`, `q(TKE)` and total averages comparable with measurement.

Implemented:

| | |
|---|---|
| TXE partition | modelling at scission with extra-deformation energies; `E*_H/TXE (A_H)` by jointed segments; `R_T = T_L/T_H (A_H)` by jointed segments |
| Level density parameter | Gilbert-Cameron with shell corrections; Egidy-Bucurescu BSFG |
| Evaporation cross section | constant, giving `T_k` analytically; Iwamoto's `σ_0(1 + α/√ε)`, giving a transcendental equation solved by bisection |
| s-wave strength function | `S_0(A)` parameterized by jointed segments over RIPL-3 resonance data |
| Derived quantities | `ν(A,Z,TKE)`, `ν_pair`, `A_p`, `Y(A_p)`, `Y(Z,A_p)`, `Y(N_p)`, `KE_p`, `TKE_p`, `R_T`, `T_0`, `Q`, `TXE`, `E*` |
| Distributions | `P(T)`, `P_k(T)`, `P(E_r)`, `P(⟨ε⟩)`, `P(ν)`, `P(E*)` |
| Prompt fission neutron spectrum | laboratory spectrum per fragment and sequence-averaged CM spectrum, with the ratio to an equivalent Maxwellian at `T_M = (2/3)⟨E⟩` |

## Files

| File | Role |
|---|---|
| `main_DSE.jl` | Entry point. Includes the others in order and issues the pipeline calls |
| `init.jl` | Startup, `Distribution` and `Distribution_unidym` definitions |
| `input.jl` | Fissioning system, fragment and TKE ranges, method selectors |
| `options.jl` | Input file names and layouts, output flags, segment points, histogram steps |
| `error_handling.jl` | Validation of the two files above |
| `aux_func.jl` | Data loading, `Q`, `S_n`, `TXE`, fragmentation domain, sequence averaging |
| `density_parameters.jl` | `a(A,Z)` by either systematic |
| `txe_partitioning.jl` | The three partition methods and the jointed-segment helper |
| `dse_eq_solvers.jl` | The sequential equations for both cross sections, `S_0(A)`, the transcendental solver |
| `write_output.jl` | `readme.out` and the `main_DSE.out` record file |
| `secondary_outputs.jl` | `Y(A,Z,TKE)`, all yield averaging, post-emission quantities, the `.dat` exports |
| `neutron_spectrum.jl` | Laboratory and CM spectra, Maxwellian normalization |
| `plots.jl` | Figures |
| `input_data/` | AME2020 mass excesses, Wahl `ΔZ(A)` and `rms(A)`, Gilbert-Cameron shell corrections, Möller `β₂`, RIPL-3 resonances, and `Y(A,TKE)` for ²⁵²Cf, ²³⁵U and ²³⁹Pu — documented in `input_data/README.md`, files not tracked |
| `Data conversion/` | One-off scripts that built some of the input files |

## Input data

The data files under `input_data/` are third-party scientific data — mass evaluations, charge
distribution and shell correction systematics, a resonance parameter library, and experimental
fragment mass yields. None of it originates with this work and the MIT licence on the code does
not extend to it, so the files are not tracked here. `input_data/README.md` lists every one of
them with its source, and they remain in the import commit:

```
git checkout bdba1d9 -- input_data/
```

That restores the exact bytes the thesis calculations ran against.

## Running it

Restore `input_data/` as above, edit `input.jl` and `options.jl`, then `julia main_DSE.jl` from
this directory.

Expect friction. There is no `Project.toml`, so `DataFrames`, `CSV`, `Roots`, `QuadGK`,
`Trapz`, `Plots` and `LaTeXStrings` must already be present. The code `cd`s into `input_data/`
and back, leaving the working directory changed. `A₀` is incremented in place for `(n,f)`, so
re-including the chain in one session increments it again. Output directories are created with
`mkdir` on a single existence check, so a partially present tree fails.

Runtime is dominated by boolean-mask filtering, not by the physics — see below.

## Known defects

Recorded because they motivate the rewrite, not because they were fixed.

**Performance.** The dominant pattern throughout is
`q.Value[(q.A .== A) .& (q.Z .== Z) .& (q.TKE .== TKE)][1]` inside triple loops, and the
`isassigned(x[mask], 1)` idiom, which allocates a whole filtered array to test emptiness. It
recurs roughly a hundred times in `secondary_outputs.jl` alone and makes the averaging routines
O(N²)–O(N³). `Probability_of_occurrence` runs a full nested scan *per histogram bin*; at
`ΔT_k = 0.01` MeV it is the hottest routine in the program. `Separation_energy` and
`density_parameter` are re-evaluated per sequence per configuration although they depend only on
`(A−k+1, Z)`.

**Type instability.** `Distribution.Value` and `Distribution_unidym.Argument` are untyped
fields, so inference fails in the core data structure. Everything executes at top level, so
every global read inside the hot loops is boxed and non-inferable. `using` statements and
`const` declarations sit inside `if` blocks.

**Correctness.** `Fragmentation_domain` sorts the `A` and `Z` columns independently and then
looks values up by `(A,Z)` taking the first match, which can mispair rows and silently drops
duplicates. The symmetric-fragment guard tests `A_L != A_H_min` where `A_L != A_H` is meant;
equivalent only when `A_H_min = A₀/2`. `Q_value_released`, `Separation_energy`,
`Total_excitation_energy` and `density_parameter` return a `(value, σ)` tuple on success and a
bare `NaN` on failure, which works only because scalar indexing of a `Float64` is legal.
`Sort_q_Argument` writes back only as many entries as there are unique arguments, leaving a
stale tail when duplicates exist.

**The output flags are not independent.** `y_A` is built only under `secondary_output_Yield`
but consumed under `secondary_output_TXE_Q`, `_Ap` and `_E_excitation`; `Eᵣ_A_Z_TKE` exists only
if `secondary_output_T` ran but is read by `secondary_output_Eᵣ`. Some combinations error, and
some silently pick up a value left over from another loop.

**Output handling.** Every writer truncates. `Write_seq_output_Y` overwrites the `main_DSE.out`
that the `p(Z,A)` path wrote earlier in the same run, and the `Average_quantities` file is
opened outside every flag guard, so an all-flags-false run replaces a good file with an empty
one. Nothing carries a run identifier, so two runs differing only in method overwrite each other
unless `file_output_identifier` is changed by hand.

**Uncertainties stop early.** `σ` is computed for `Q`, `S_n` and `TXE`, but `TXE` is always
evaluated with `σ_TKE = 0` and the `σ` field is left empty through the solver chain. The
uncertainties that reach the output come from the mass excess table alone.

**Figures.** `Plots.jl` with the GR backend, PNG only at a hardcoded 750×500 px, a title on
nearly every panel, `minorgrid` on, no font sizing, automatic series colours and markers,
annotations anchored to magic indices, and two backends active in one run.

## Deliberate, not defects

Four things that look like bugs and are not. They are preserved in the rewrite.

- **`A_p = A − n` uses `n`, the integer sequence count**, while `ν = (n+1)/2`. These are
  different quantities: `n` is the number of emissions at one configuration, `ν` is the prompt
  neutron multiplicity, an observable. Integer `n` is the reason for treating emission
  sequentially at all.
- **`SeqAvg_A_Z_TKE` divides by `k_max`**, the sequences that admit a temperature solution, not
  by `n`. Padding with zeros over the remainder would dilute a temperature with values that are
  not temperatures.
- **The per-sequence figures stop two sequences short of the maximum.** The last sequences carry
  negligible probability; plotting them adds invisible curves and dead legend entries.
- **`T_k` below 0.01 MeV is rejected.** Below roughly that value the level density leaves the
  Fermi-gas regime for the constant-temperature one, so a Fermi-gas temperature there has no
  meaning.

## The thesis

`MSc_thesis_2023.pdf` is the submitted thesis rebuilt from its LaTeX source with orthographic
and grammatical corrections, Romanian quotation marks, a repointed source-code citation, and one
corrected intermediate expression in the appendix proof of root uniqueness — verified not to
affect that proof's conclusion. The argument, structure, figures and numerical results are as
submitted.

Two typographic passes follow. The first: margin kerning and font expansion, vertical spacing
left ragged rather than stretched, the one side-by-side figure pair that overran the text block
brought back to the 8 cm the document uses everywhere else, and the two moment integrals of the
appendix broken after their substitution annotation. The submitted source carried 8 overfull
lines, the worst nearly a centimetre into the margin, and 6 pages padded to justify their
height; the rebuild carries none.

The second indents the first line of each paragraph. The submitted source sets no indentation
and no vertical separation between paragraphs, so they run together — three consecutive
paragraphs read as one block in the chapter on code structure. Indenting restores the separator
without changing vertical spacing.

The document is still 48 pages and no wording changed.

## Licence

Code under the MIT Licence, see `LICENSE`. The thesis PDF is licensed CC BY 4.0.
