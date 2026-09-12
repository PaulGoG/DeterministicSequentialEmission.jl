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
| TXE partition | `E*_H/TXE (A_H)` by jointed segments |

| File | Content |
|---|---|
| `CF52_Q_AH_Z.dat` | `Q(A_H, Z_H)` in MeV, 245 rows — one per fragmentation |
| `CF52_Q_AH.dat` | `Q(A_H)` in MeV, averaged over the isobaric charge distribution |

These are not inputs. They exist so that a change in the physics shows up as a failing test
rather than as a number nobody checked.
