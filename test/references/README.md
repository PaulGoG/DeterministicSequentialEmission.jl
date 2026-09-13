# Reference data

Output of the implementation submitted with the 2023 MSc thesis, kept here so the rewrite can
be pinned against it. The configuration these files come from:

| | |
|---|---|
| Fissioning system | ²⁵²Cf(SF) |
| Heavy fragment masses | 126 … 174 |
| TKE grid | 130 … 230 MeV, step 2 MeV |
| Charge numbers per mass | 5 |
| Mass excesses | `AME2020.ANA` |
| Charge polarization | `DeltaZ_rms_A.CF52` |
| Evaporation cross section | variable |
| Level density parameter | back-shifted Fermi gas |

Two runs were archived, one per TXE partition, and they differ in nothing else. Their
parameterizations, transcribed from their own provenance headers:

| Partition | Segments as run |
|---|---|
| `E*_H/TXE (A_H)` | `(126, 0.5) (130, 0.14) (150, 0.607) (174, 0.865)` |
| `R_T (A_H)` | `(126, 1.0) (130, 1.734) (134, 1.28) (138, 1.28) (150, 0.9) (167, 0.68) (174, 0.59)` |

| File | Content |
|---|---|
| `CF52_Q_AH_Z.dat` | `Q(A_H, Z_H)` in MeV, 245 rows — one per fragmentation |
| `CF52_Q_AH.dat` | `Q(A_H)` in MeV, averaged over the isobaric charge distribution |
| `CF52_level_density.dat` | `a(A, Z)` in MeV⁻¹ over the fragment domain |
| `CF52_excitation_PARAM.dat` | `E*_L`, `E*_H` in MeV under the `E*_H/TXE (A_H)` partition |
| `CF52_excitation_RT.dat` | `E*_L`, `E*_H` in MeV under the `R_T (A_H)` partition |
| `CF52_sequences.dat` | per-sequence `Sₙ`, `a`, `T`, `Ē_r`, `⟨ε⟩` |

The two excitation files sample the same configurations — `A_H` ∈ {126, 134, 144, 154, 164,
174}, all five charge numbers, the whole TKE grid, 1200 rows each — so the partitions can be
compared row for row. The values are rounded to eight decimals on write, which is why the tests
assert agreement below `1e-8` rather than to machine precision: `5e-9` is exact agreement.

The sample is a regression guard, not the extent of what has been checked. Both partitions have
been swept over **every** record of their archived run — 9819 and 9788 — and reproduce each one
to the same `5e-9` MeV.

These are not inputs. They exist so that a change in the physics shows up as a failing test
rather than as a number nobody checked.
