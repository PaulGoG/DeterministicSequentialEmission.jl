# Reference data

Output of the implementation submitted with the 2023 MSc thesis, kept here so the rewrite can
be pinned against it. The configuration these files come from:

| | |
|---|---|
| Fissioning system | ²⁵²Cf(SF) |
| Heavy fragment masses | 126 … 174 |
| TKE grid | 130 … 230 MeV, step 2 MeV |
| Charge numbers per mass | 5 |
| Mass excesses | the 2020 atomic mass evaluation |
| Charge distribution | the per-reaction `ΔZ(A)`, `σ_Z(A)` fit for ²⁵²Cf |
| Evaporation cross section | variable |
| Level density parameter | back-shifted Fermi gas |

Two runs were archived, one per TXE partition, and they differ in nothing else. Their
parameterizations, transcribed from their own provenance headers:

| Partition | Segments as run |
|---|---|
| `E*_H/TXE (A_H)` | `(126, 0.5) (130, 0.14) (150, 0.607) (174, 0.865)` |
| `R_T (A_H)` | `(126, 1.0) (130, 1.734) (134, 1.28) (138, 1.28) (150, 0.9) (167, 0.68) (174, 0.59)` |

They are laid out one directory per fissioning system, and named by the same
`<quantity>_vs_<abscissa>` rule as every other table in the toolchain. The partition a file was
produced under is part of its name, because two of them differ in nothing else.

| File, under `Cf252_sf/` | Content |
|---|---|
| `Q_vs_A_H_Z_H.dat` | `Q(A_H, Z_H)` in MeV, 245 rows — one per fragmentation |
| `Q_vs_A_H.dat` | `Q(A_H)` in MeV, averaged over the isobaric charge distribution |
| `a_vs_A_H_Z_H.dat` | `a(A, Z)` in MeV⁻¹ over the fragment domain |
| `E_star_vs_A_H_Z_H_TKE_excitation_ratio.dat` | `E*_L`, `E*_H` in MeV under the `E*_H/TXE (A_H)` partition |
| `E_star_vs_A_H_Z_H_TKE_temperature_ratio.dat` | `E*_L`, `E*_H` in MeV under the `R_T (A_H)` partition |
| `sequences_excitation_ratio.dat`, `sequences_temperature_ratio.dat` | per-sequence `Sₙ`, `a`, `T`, `Ē_r`, `⟨ε⟩` |
| `nu_vs_A.dat`, `T_vs_A.dat`, `E_star_vs_A.dat`, `E_star_vs_TKE.dat` | the averaged observables |
| `Y_vs_A.dat`, `Y_vs_Z.dat`, `Y_vs_TKE.dat`, `TKE_vs_A_H.dat` | the fragment yield projections |
| `Y_vs_A_p_excitation_ratio.dat`, `Y_vs_A_p_temperature_ratio.dat` | the post-neutron mass yield `Y_p(A_p)` under each partition, `A_p = A − n`, 98 rows, summing to 200 |

The two excitation files sample the same configurations — `A_H` ∈ {126, 134, 144, 154, 164,
174}, all five charge numbers, the whole TKE grid, 1200 rows each — so the partitions can be
compared row for row. The values are rounded to eight decimals on write, which is why the tests
assert agreement below `1e-8` rather than to machine precision: `5e-9` is exact agreement.

Only the header lines were rewritten when the naming convention was applied; every data row is
byte-identical to what the archived run produced, and `.gitattributes` keeps these files exempt
from line-ending normalisation so they stay that way. The two post-neutron tables were written
with CRLF terminators and are kept with LF; their data rows are otherwise the archived ones.

The sample is a regression guard, not the extent of what has been checked. Both partitions have
been swept over **every** record of their archived run — 9819 and 9788 — and reproduce each one
to the same `5e-9` MeV.

These are not inputs. They exist so that a change in the physics shows up as a failing test
rather than as a number nobody checked.
