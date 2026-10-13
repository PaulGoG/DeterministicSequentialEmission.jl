# Results: the model against measurement

The figures here are drawn by `plotting/figures.jl` from the tables one run per system wrote
under `data/sims/`, with nothing solved again and nothing fitted to the data shown. Each figure
is exported as PDF and PNG beside a `.toml` sidecar naming, per system, the run drawn, the commit
its metadata records and the measured files, so every panel traces to a configuration and a
commit. The measured side comes from EXFOR through ExforFissionData.jl v0.3.0, laid out as
`data/README.md` describes; the temperature ratios the runs were solved with are the
FissionTemperatureRatio.jl v0.3.0 trends staged beside them. The two diagnostic panels the
script also writes here, the residual quantities and the charge distribution, have no measured
counterpart and are not tracked.

```
julia scripts/run.jl config/<system>.toml    # one run per system, under data/sims/
julia plotting/figures.jl                    # redraws everything here from those runs
```

| System | ⟨ν⟩ | evaluated | | ⟨E⟩ [MeV] | ν̄/ν̄(²⁵²Cf) |
|---|---|---|---|---|---|
| ²⁵²Cf(sf) | 3.557 | 3.757 | −5.3 % | 2.164 (measured 2.13) | 1 |
| ²³⁵U(n_th,f) | 2.535 | 2.414 | +5.0 % | 1.954 | 0.713 (evaluated 0.643) |
| ²³⁹Pu(n_th,f) | 2.813 | 2.875 | −2.1 % | 2.046 | 0.791 (evaluated 0.765) |
| ²³³U(n_th,f) | 2.631 | 2.497 | +5.4 % | 2.007 | 0.740 (evaluated 0.665) |

The evaluated multiplicities are ENDF/B-VIII.0 for the induced systems and the Axton standard
for ²⁵²Cf. Datasets are labelled by author and the year of their EXFOR entry; the references
below give the publication, which for conference proceedings can be a later year.

## Prompt neutron multiplicity against fragment mass

![Prompt neutron multiplicity against fragment mass, four fissioning systems](multiplicity_vs_experiment.png)

The sawtooth is the signature observable and it is not fitted: the total excitation energy is
shared between the fragments by a temperature ratio `R_T(A_H)` extracted from measured `ν(A)`
by the companion package, the emission chain is solved configuration by configuration, and the
result is contracted against the measured fragment yield `Y(A, TKE)`. Each panel states the
total multiplicity against the evaluated value. The sets drawn are a fixed selection per system,
for legibility; `data/README.md` lists the whole retrieval.

| System | Dataset | EXFOR | Reference |
|---|---|---|---|
| ²⁵²Cf(sf) | Zeynalov 2011 | 23118006 | J. Korean Phys. Soc. **59**, 1396 (2011), doi:10.3938/jkps.59.1396 |
| | Budtz-Jørgensen 1988 | 23175008 | Nucl. Phys. A **490**, 307 (1988), doi:10.1016/0375-9474(88)90508-8 |
| | Göök 2014 | 23268005 | Phys. Rev. C **90**, 064611 (2014), doi:10.1103/PhysRevC.90.064611 |
| | Al-Adili 2020 | 23764010 | Phys. Rev. C **102**, 064610 (2020), doi:10.1103/PhysRevC.102.064610 |
| | Vorobiev 2001 | 41425014 | ISINN-9, JINR E3-2001-192, p. 288 (2001) |
| ²³⁵U(n_th,f) | Maslin 1967 | 21095002 | AWRE O-43/67 (1967) |
| | Nishio 1998 | 22464004 | Nucl. Phys. A **632**, 540 (1998), doi:10.1016/S0375-9474(98)00008-6 |
| | Al-Adili 2020 | 23764004 | as above |
| | Boldeman 1971 | 30909004 | Aust. J. Phys. **24**, 821 (1971), doi:10.1071/PH710821 |
| | Vorobyev 2010 | 41516012 | ISINN-17, p. 60 (2010) |
| ²³⁹Pu(n_th,f) | Fraser 1966 | 14369002 | Annu. Rev. Nucl. Sci. **16**, 379 (1966), doi:10.1146/annurev.ns.16.120166.002115 |
| | Nishio 1995 | 23012008 | J. Nucl. Sci. Technol. **32**, 404 (1995), doi:10.1080/18811248.1995.9731725 |
| | Apalin 1965 | 41397002 | Nucl. Phys. **71**, 553 (1965), doi:10.1016/0029-5582(65)90765-0 |
| ²³³U(n_th,f) | Fraser 1966 | 14369004 | as above |
| | Nishio 1998 | 22660005 | J. Nucl. Sci. Technol. **35**, 631 (1998), doi:10.1080/18811248.1998.9733919 |
| | Takamiya 1999 | 23595005 | J. Radioanal. Nucl. Chem. **239**, 117 (1999) |
| | Apalin 1965 | 41397004 | as above |

## Prompt fission neutron spectra

![Prompt fission neutron spectra: 252-Cf as a ratio to a Maxwellian, the thermal systems as ratios to 252-Cf](prompt_neutron_spectrum.png)

The ²⁵²Cf(sf) spectrum is drawn as the ratio of the model's unit-area laboratory spectrum to
a Maxwellian at 1.42 MeV, which is how every retrieved set reports it, so model and measurement
are both absolute. The three thermal systems were measured as the ratio of their spectrum per
fission to that of ²⁵²Cf(sf) with one detector, and the model is drawn the same way,
`ν̄ N(E) / ν̄_Cf N_Cf(E)` with its own multiplicities: the offset at low energy is the
multiplicity error of the pair, stated in each panel, and the fall with energy is the thermal
spectrum being the softer. Sets tabulated as ²⁵²Cf over the system are inverted, the two energy
ranges of one entry are joined under one label, and a set its record flags as superseded is
skipped.

| System | Dataset | EXFOR | Reference |
|---|---|---|---|
| ²⁵²Cf(sf) | Poenitz 1982 | 14278002 | Proc. Int. Conf. Nuclear Data for Science and Technology, Antwerp 1982, p. 465 |
| | Boldeman 1986 | 30775002 | Nucl. Sci. Eng. **93**, 181 (1986), doi:10.13182/NSE86-2 |
| | Blinov 1973 | 40418008 | Neutron Physics, Kiev 1973, vol. 4, p. 138 |
| | Lajtai 1990 | 41158003 | Nucl. Instrum. Methods A **293**, 555 (1990), doi:10.1016/0168-9002(90)90323-X |
| ²³⁵U(n_th,f) | Nefedov 1983 | 40871011, 40871012 | Neutron Physics, Kiev 1983, vol. 2, p. 285 |
| | Starostov 1983 | 40872007 | Neutron Physics, Kiev 1983, vol. 2, p. 290 |
| | Vorobyev 2013 | 41597002 | Vopr. At. Nauki Tekh., Ser. Yad. Konst. (1–2), 37 (2013); supersedes 41516017 |
| ²³⁹Pu(n_th,f) | Nefedov 1983 | 40871009, 40871010 | as above |
| | Starostov 1983 | 40872006 | as above |
| ²³³U(n_th,f) | Nefedov 1983 | 40871013 | as above |
| | Starostov 1983 | 40872008 | as above |

## Mean centre-of-mass neutron energy against fragment mass

![Mean centre-of-mass energy of the prompt neutrons against fragment mass, four fissioning systems](mean_neutron_energy_vs_experiment.png)

Both sides are means in the centre of mass of the emitting fragment, as the retrieval records
for every set drawn. The model's `⟨ε⟩(A)` is the sequence average the energy balance returns,
over all neutron energies; the measured means are formed from spectra above a detector
threshold — 0.52 MeV in the laboratory for Bowman, 0.3 MeV for Budtz-Jørgensen, whose mean is
the first moment of the form fitted to the spectrum, 0.2 to 0.5 MeV for the others — which
raises them slightly where the fragment's neutrons are slow. The model follows the measured
values along both wings and fails in two places: the band `A ≈ 126`–`134` of heavy fragments
near the `N = 82` closure, which barely emit in the model and whose few neutrons are slow, and
the lightest fragments, `A < 90`, where it lies below the data.

| System | Dataset | EXFOR | Reference |
|---|---|---|---|
| ²⁵²Cf(sf) | Bowman 1963 | 14065003 | Phys. Rev. **129**, 2133 (1963), doi:10.1103/PhysRev.129.2133 |
| | Budtz-Jørgensen 1988 | 23175012 | Nucl. Phys. A **490**, 307 (1988), doi:10.1016/0375-9474(88)90508-8 |
| | Göök 2014 | 23268011 | Phys. Rev. C **90**, 064611 (2014), doi:10.1103/PhysRevC.90.064611 |
| | Piksaykin 1977 | 41689005 | Yad. Fiz. **25**, 723 (1977) |
| ²³⁵U(n_th,f) | Nishio 1998 | 22464003 | Nucl. Phys. A **632**, 540 (1998), doi:10.1016/S0375-9474(98)00008-6 |
| | Batenkov 2004 | 41502008 | AIP Conf. Proc. **769**, 1003 (2005), doi:10.1063/1.1945175 |
| ²³⁹Pu(n_th,f) | Tsuchiya 2000 | 22650008 | J. Nucl. Sci. Technol. **37**, 941 (2000), doi:10.1080/18811248.2000.9714976 |
| | Batenkov 2004 | 41502009 | as above |
| ²³³U(n_th,f) | Fraser 1966 | 14369005 | Annu. Rev. Nucl. Sci. **16**, 379 (1966), doi:10.1146/annurev.ns.16.120166.002115 |
| | Nishio 1998 | 22660003 | J. Nucl. Sci. Technol. **35**, 631 (1998), doi:10.1080/18811248.1998.9733919 |

## Not drawn yet

The runs write more than the figures show, and the retrieval holds measurements the figures do
not yet use. Each of these enters this file with its figure and its dataset table when it is
ready.

- Post-emission distributions: the fragment mass and charge yields after neutron emission,
  against the independent and post-neutron yields of the evaluations. The model gives an integer
  `A_p = A − n` at every configuration, which is what makes it suitable for them; the port from
  the thesis code is pending, as the README's status table says.
- The neutron number distribution `P(ν)`, against the retrieved `P_vs_nu` sets; the model's is
  known to be too broad in all four systems.
- `ν(TKE)` and `⟨ε⟩(TKE)`, which the runs write and the retrieval holds as `nu_vs_A_TKE` and
  `eps_vs_A_TKE` matrices for ²⁵²Cf and ²³³U.
- The thermal systems' spectra as ratios to a Maxwellian, which the retrieval holds for all
  three.
