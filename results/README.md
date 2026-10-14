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

## Post-neutron product mass yield

![Post-neutron product mass yield of the four systems against the evaluated independent yields](product_mass_yield.png)

The product mass yield `Y(A_p)` is the pre-neutron yield carried through the emission chains.
A fragment whose chain allows `n` sequences ends on each of its `n` residuals with probability
`1/n`, so its products lie on `A − 1 … A − n`, and the mass its products have lost equals its
multiplicity exactly — 1.78 per fragment for ²⁵²Cf, the ⟨ν⟩ = 3.56 per fission of the sawtooth
figure. No parameter enters. The evaluated independent yields of ENDF/B-VIII.0 and JEFF-3.3,
summed over charge, are lines with their uncertainty as a band: they are evaluations, not
measurements. The run's own pre-neutron `Y(A)`, dotted, shows what the emission did: it moves
both humps inward by the mass lost, sharpens them, and builds the even-`Z` structure the
evaluations carry — the peaks at `A_p = 134` and `138` of ²³⁵U, at `134`–`136` of ²³⁹Pu, at
`136` of ²³³U — out of a pre-neutron yield that has none. The centroids of the two humps agree
with the evaluations within 0.8 u in every system.

| System | ⟨A_p⟩ light: model / ENDF/B-VIII.0 / JEFF-3.3 | ⟨A_p⟩ heavy | largest `Y(A_p)`: model / ENDF/B-VIII.0 |
|---|---|---|---|
| ²⁵²Cf(sf) | 106.9 / 106.1 / 106.2 | 141.7 / 141.9 / 142.1 | 106 (6.4 %) / 107 (6.6 %) |
| ²³⁵U(n_th,f) | 94.9 / 94.9 / 94.9 | 138.6 / 138.7 / 138.7 | 134 (7.7 %) / 134 (7.9 %) |
| ²³⁹Pu(n_th,f) | 98.9 / 99.0 / 99.0 | 138.3 / 138.1 / 138.2 | 134 (6.9 %) / 134 (7.7 %) |
| ²³³U(n_th,f) | 93.5 / 93.4 / 93.3 | 137.9 / 138.2 / 138.2 | 136 (7.7 %) / 93 (7.0 %) |

Where it departs: the heavy hump of ²³⁹Pu peaks at 6.9 % against 7.6–7.7 % evaluated, the
light hump of ²³³U at 7.5 % against 7.0 %, and the outer flank of the heavy hump of ²⁵²Cf falls
below the evaluations beyond `A_p ≈ 150`. The published prescription of the model, which
subtracts the full number of possible sequences from every fragment, is one key away,
`[products] count = "sequence_count"`; for ²⁵²Cf it removes 2.58 neutrons per fragment, 45 %
more than the multiplicity, and moves the humps a further 0.8 u inward.

| Library | Reference |
|---|---|
| ENDF/B-VIII.0, NFY and SFY sublibraries, MF8 MT454 at 0.0253 eV or spontaneous | Brown *et al.*, Nucl. Data Sheets **148**, 1 (2018), doi:10.1016/j.nds.2018.02.001 |
| JEFF-3.3, NFY and SFY, the same | Plompen *et al.*, Eur. Phys. J. A **56**, 181 (2020), doi:10.1140/epja/s10050-020-00141-9; NEA Data Bank record doi:10.82555/nhfqy-hvz09 |

Isomeric states are summed and uncertainties added in quadrature; JEFF-3.3's ternary particles
are left out. The measured post-neutron mass yields and independent yields of EXFOR are not
retrieved yet; when they are, they enter this figure as markers.

The figures below draw the independent yields `Y(Z, A_p)` the runs write as `Y_vs_Z_A_p.dat`,
regrouped by charge, by isotope and by neutron number, against the same two evaluations,
`Y_vs_Z_A_p/<library>.dat` summed the same way, and the kinetic energies of the products. All
of them share the evaluations' colours and line styles, and the model is the black line.

## Elemental yield of the products

![Elemental yield of the post-neutron products of the four systems against the evaluated independent yields](product_charge_yield.png)

`Y(Z) = Σ_{A_p} Y(Z, A_p)`, on a logarithmic axis down to 0.001 %. Emission leaves the charge
alone, so the model's `Y(Z)` is the charge distribution of the solved configurations: the `Z_p`
model of Wahl, At. Data Nucl. Data Tables **39**, 1 (1988), doi:10.1016/0092-640X(88)90016-2,
with the reaction's own parameters, carried through the pre-neutron yield. Its even–odd
staggering is that model's proton factor `F_Z`, fitted to measured independent yields, and not a
result of the emission. Each panel states the global proton even–odd effect
`δ_p = (Y_e − Y_o)/(Y_e + Y_o)`, `Y_e` and `Y_o` the yields summed over the even and the odd
charges.

| System | `δ_p`: model / ENDF/B-VIII.0 / JEFF-3.3 | `Y(Z₀/2)` [%]: model / ENDF/B-VIII.0 / JEFF-3.3 |
|---|---|---|
| ²⁵²Cf(sf) | 5 / 5 / 13 % | 0.77 / 0.073 / 0.31 |
| ²³⁵U(n_th,f) | 23 / 24 / 28 % | 0.052 / 0.035 / 0.047 |
| ²³⁹Pu(n_th,f) | 14 / 14 / 17 % | 0.13 / 0.071 / 0.12 |
| ²³³U(n_th,f) | 24 / 24 / 26 % | 0.013 / 0.041 / 0.060 |

The model's `Y(Z)` is symmetric about `Z₀/2`, each split giving one light and one heavy charge,
and the evaluations nearly so. Over both humps the model lies with the two evaluations across
three decades. The valley is where it departs: for ²⁵²Cf the model gives the symmetric charge
`Z = 49` ten times the yield of ENDF/B-VIII.0, and for ²³³U a third of it.

## Isotopic distributions

![Isotopic distributions of the four most abundant light and heavy elements, 252-Cf(sf)](product_independent_yield_Cf252_sf.png)
![Isotopic distributions of the four most abundant light and heavy elements, 235-U(nth,f)](product_independent_yield_U235_nth.png)
![Isotopic distributions of the four most abundant light and heavy elements, 239-Pu(nth,f)](product_independent_yield_Pu239_nth.png)
![Isotopic distributions of the four most abundant light and heavy elements, 233-U(nth,f)](product_independent_yield_U233_nth.png)

One figure per system: in a 2 × 2 grid of systems the eight panels of one system would be some
140 px wide at the 1200 px of the grid, too narrow for their tick labels. Each panel is `Y(Z, A_p)` at fixed `Z` against `A_p`,
for the four light and the four heavy elements of largest evaluated yield, the mean of the two
libraries, so that the choice does not move with the model. The model is the line, the
evaluations are markers: circles ENDF/B-VIII.0, squares JEFF-3.3. The model places the centroid
and the width of every distribution drawn and follows the evaluations over three decades. Its
distributions are smooth where the evaluations stagger from one isotope to the next: Wahl's
neutron factor `F_N` acts on the pre-neutron fragments, and emission, which moves each fragment
by its own number of neutrons, washes it out. The far wings, below 0.01 %, are where the two
libraries part from each other as well.

## Neutron number distribution of the products

![Yield against product neutron number of the four systems against the evaluated independent yields](product_neutron_number_yield.png)

`Y(N_p) = Σ_{A_p − Z = N_p} Y(Z, A_p)`, the evaluations with their uncertainty as a band and the
closed shells `N = 50` and `82` as guides. The heavy products of the three thermal systems pile
up at `N_p = 82`: a heavy fragment just above the closure loses its weakly bound neutrons down to
it, and one at the closure barely emits. The model reproduces the peak for ²³⁵U and ²³⁹Pu and
overshoots it for ²³³U. ²⁵²Cf, whose heavy fragments lie further above the closure, has its heavy
hump at `N_p ≈ 86`, and there too the model puts too much yield on the closure itself.

| System | `Y(N_p = 82)` [%]: model / ENDF/B-VIII.0 / JEFF-3.3 | ⟨N_p⟩ light, heavy: model / ENDF/B-VIII.0 |
|---|---|---|
| ²⁵²Cf(sf) | 9.0 / 6.9 / 7.9 | 64.2, 86.4 / 63.6, 86.5 |
| ²³⁵U(n_th,f) | 14.7 / 14.4 / 13.6 | 56.9, 84.6 / 56.9, 84.7 |
| ²³⁹Pu(n_th,f) | 13.7 / 14.5 / 14.2 | 59.1, 84.2 / 59.1, 84.0 |
| ²³³U(n_th,f) | 15.4 / 12.5 / 13.5 | 55.7, 83.7 / 55.6, 84.0 |

## Kinetic energy of the products

![Mean kinetic energy of the post-neutron products against product mass, four fissioning systems](product_kinetic_energy.png)

`E_K_p(A_p)`, the mean over the products of mass `A_p` of `E_K_p = E_K A_p/A`,
`E_K = TKE (A₀ − A)/A₀` being the kinetic energy of the fragment, beside the pre-neutron
`E_K(A)` of the yield the run was solved with as a dashed guide, against the fragment mass. That
guide is `⟨TKE⟩(A) (A₀ − A)/A₀`, `⟨TKE⟩(A)` read from the input the run's configuration names
and averaged with its complement as the run did. Emission isotropic in the fragment frame
leaves the velocity unchanged on average, so a product carries the fraction `A_p/A` of its
fragment's energy, and the curve is the guide moved down by 1 to 2 % and to lower mass by the
multiplicity. Masses whose product yield is below 0.01 % are not drawn. The measured
post-neutron kinetic energies, EXFOR `MASS,SEC,KE`, are not retrieved yet; ExforFissionData.jl
writes them under `E_K_p_vs_A_p/`, none for ²⁵²Cf, and `plotting/figures.jl` draws a set as
markers once its accession is listed in `SHOWN`.

![Mean kinetic energy of the post-neutron products against product charge, four fissioning systems](product_kinetic_energy_by_charge.png)

`E_K_p(Z)`, the same mean over the products of charge `Z`, in a figure of its own: it would
share the cell of each system with `E_K_p(A_p)` only with a second abscissa and its labels in
every cell. It has no measured or evaluated counterpart here.

## Total kinetic energy of the products

![Mean post-neutron total kinetic energy against heavy product mass, four fissioning systems](product_total_kinetic_energy.png)

`TKE_p(A_p)`, the mean total kinetic energy of the two products against the heavy product mass,
with the pre-neutron `⟨TKE⟩(A)` input as a dashed guide against the heavy fragment mass. The
two products of a pair are taken as independent, so that `TKE_p` is the mean over both ends of
both chains. Masses whose heavy-product yield is below 0.01 % are not drawn. Each panel states
both means over the run: emission takes 2.1 to 2.7 MeV from the total kinetic energy.

| System | ⟨TKE⟩ [MeV] | ⟨TKE_p⟩ [MeV] |
|---|---|---|
| ²⁵²Cf(sf) | 184.2 | 181.5 |
| ²³⁵U(n_th,f) | 170.4 | 168.3 |
| ²³⁹Pu(n_th,f) | 177.6 | 175.4 |
| ²³³U(n_th,f) | 169.8 | 167.7 |

## Not drawn yet

The runs write more than the figures show, and the retrieval holds measurements the figures do
not yet use. Each of these enters this file with its figure and its dataset table when it is
ready.

- The measured post-neutron yields and kinetic energies, which ExforFissionData.jl v0.4.0
  retrieves as `Y_vs_A_p`, `Y_vs_Z_A_p` and `E_K_p_vs_A_p`.
- The neutron number distribution `P(ν)`, against the retrieved `P_vs_nu` sets; the model's is
  known to be too broad in all four systems.
- `ν(TKE)` and `⟨ε⟩(TKE)`, which the runs write and the retrieval holds as `nu_vs_A_TKE` and
  `eps_vs_A_TKE` matrices for ²⁵²Cf and ²³³U.
- The thermal systems' spectra as ratios to a Maxwellian, which the retrieval holds for all
  three.
