# Input data

Almost nothing here is version-controlled. The files are third-party scientific data — an
atomic mass evaluation, charge distribution and shell correction systematics, a resonance
parameter library, and experimental fragment mass yields. None of it originates with this
package, and the MIT licence covering the source code does not extend to it: each collection
carries the terms and the attribution of its own source, recorded below.

**The one exception is `mass_excess/ame2020.dat`**, which is shipped. It is an openly published
evaluation, it is small, and keeping it means the tests pinned against the archived reference
run — `Q`, the level density parameter, both TXE partitions, the emission sequences — execute on
a bare clone. The tests that need any of the rest are skipped when it is absent, with a warning
naming the file.

## Layout

| Path | Content | Format |
|---|---|---|
| `mass_excess/ame2020.dat` | mass excesses, 3558 nuclides — **shipped** | header, then `Z A element Δ σΔ`, keV |
| `level_density/shell_corrections_gilbert_cameron.dat` | `S(N)`, `S(Z)` against nucleon number, 11 … 150 | header, then `n S(N) S(Z)`, MeV |
| `level_density/resonances_s_wave.dat` | s-wave average resonance parameters, the source of `S₀(A)` | commented header; `S₀` in units of 10⁻⁴ |
| `charge_distribution/{Cf252,U235,Pu239}.dat` | charge polarization and Gaussian dispersion | header, then `A ΔZ rms` |
| `deformation/ground_state_beta2_moller.dat` | ground-state quadrupole deformations | header, then `Z A β₂` |
| `deformation/extra_deformation_U235.dat` | `ΔE_def` for the scission-modelling partition | header, then `A Z ΔE_def`, MeV |
| `yield/<nucleus>/*.dat` | pre-neutron fragment mass yields `Y(A,TKE)` | header, then `A TKE Y [σY]` |

All files are whitespace-separated with a single header line. Readers take columns by position,
not by header text, and the uncertainty column may be omitted rather than zero-filled.

Nucleon-number columns are read through `integer_column`, which accepts an integral value
written in floating point — some of these tables were emitted by a generator and write `118.0`
where they mean `118`.

## Where to obtain it

Two of these collections come from companion packages, and that is the supported route: both
write this layout directly, and both record the provenance of what they retrieved.

**Experimental fragment mass yields and prompt neutron multiplicities** —
[`ExforFissionData.jl`](https://github.com/PaulGoG/ExforFissionData.jl) retrieves them from the
IAEA EXFOR archive and writes them beside a `retrieval.toml` naming every dataset it kept or
excluded, with the reason. A configuration and that package reproduce a retrieval exactly, since
EXFOR entries are immutable once published.

```bash
git clone https://github.com/PaulGoG/ExforFissionData.jl
cd ExforFissionData.jl
julia --project scripts/retrieve.jl config/Cf252_0f_yield_A.toml
```

**The temperature ratio `R_T(A_H)`** for the temperature-ratio TXE partition —
[`FissionTemperatureRatio.jl`](https://github.com/PaulGoG/FissionTemperatureRatio.jl) extracts it
from experimental `ν(A)` and tabulates it at every integer `A_H`. Feed the resulting
`R_T_parameterized_<set>_<run>.csv` to `read_segments`. Run it with
`ratio_averaging = "mean_of_ratios"` for this consumer; the reason is in that function's
docstring.

The remaining tables are one-off downloads:

**`mass_excess/ame2020.dat`** — the 2020 atomic mass evaluation.
W. J. Huang, M. Wang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030002 (2021);
M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi, *Chinese Physics C* **45**, 030003 (2021).
Distributed by the Atomic Mass Data Center, <https://www-nds.iaea.org/amdc/>.

**`level_density/shell_corrections_gilbert_cameron.dat`** — shell corrections `S(N)` and `S(Z)`,
Table III of A. Gilbert and A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965),
pp. 1453–1455, as redistributed in the IAEA Reference Input Parameter Library. The paper
tabulates `S(Z)` only to `Z = 98`; the file is rectangular and pads the column with zeros beyond
that, which the reader discards rather than storing as a correction of zero.

**`level_density/resonances_s_wave.dat`** — s-wave average resonance parameters.
RIPL-3: R. Capote *et al.*, *Nuclear Data Sheets* **110**, 3107 (2009). The file is byte-identical
to <https://www-nds.iaea.org/RIPL-3/resonances/resonances0.dat>.

**`charge_distribution/*.dat`** — charge polarization `ΔZ(A)` and the dispersion `rms(A)` of the
isobaric charge distribution. A. C. Wahl, *Atomic Data and Nuclear Data Tables* **38**, 1 (1988);
also tabulated in Wahl, *Systematics of Fission-Product Yields*, LA-13928 (2002).

**`deformation/ground_state_beta2_moller.dat`** — ground-state quadrupole deformations from the
finite-range liquid-drop model. P. Möller, J. R. Nix, W. D. Myers, W. J. Swiatecki, *Atomic Data
and Nuclear Data Tables* **59**, 185 (1995), via RIPL.

**`deformation/extra_deformation_U235.dat`** — not a source table. It is generated from the
ground-state deformations and the deformation energies at scission, and is the input the
scission-modelling partition needs.

## Citing

Experimental data belongs to the groups that measured it. Any published result derived from the
yields or multiplicities here should cite the measurement, not this package — the EXFOR accession
number that opens each retrieved file name is what makes those citations recoverable, and the
`retrieval.toml` beside it records the reaction code and units.
