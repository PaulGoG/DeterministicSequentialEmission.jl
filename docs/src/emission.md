```@meta
CurrentModule = DeterministicSequentialEmission
```

# Emission and spectrum

The evaporation cross section, the sequence-by-sequence solution it drives, and the prompt
neutron spectrum that follows from it. The sequential solution and the sequence average are the
model of Tudora, Hambsch and Tobosaru [TudoraHambschTobosaru2018](@cite), with the systematics of
its residual quantities [Tudora2020](@cite); the evaporation spectrum in closed form is Madland
and Nix [MadlandNix1982](@cite), its laboratory transformation and the probability of emission
along a chain follow Terrell [Terrell1959](@cite), Märten and Seeliger [Maerten1984](@cite) and
Kawano et al. [Kawano2021](@cite), and the weights of the total spectrum are the conventions of
the evaluated spectra [Capote2016](@cite).

```@autodocs
Modules = [DeterministicSequentialEmission]
Pages = ["evaporation.jl", "sequences.jl", "spectrum.jl"]
```

## Post-emission distributions

A fragment's emission chain ends on a product of the same charge, before β decay, and of mass
`A_p ≤ A`. The number of possible sequences `n` is an integer at every configuration
`(A, Z, TKE)`, so every product mass is an integer and the products are distributed without
rounding. What remains to be chosen is how many of the `n` possible neutrons a product has lost;
[`ProductCount`](@ref) carries that choice, and `[products] count` selects it.

- `emission_probability`, the default, [`EmissionProbabilityCount`](@ref). The chain ends on each
  of its `n` residuals with probability `1/n`,

  ```
  P(A_p = A − k) = 1/n,    k = 1 … n,
  ```

  and a fragment that emits nothing stays at `A_p = A`. This is the reading the multiplicity and
  the spectrum already rest on: `ν = (n+1)/2` is the mean of a count equiprobable on `1 … n`, and
  the `k`-th neutron is emitted with probability `(n+1−k)/n`. Applied to the product, it makes the
  mass the products lose equal to the multiplicity, `⟨A − A_p⟩ = ν`, configuration by
  configuration. It is the uniform average over the `n` endings of a chain that eq. (12) of
  [Tudora2022](@cite) takes for the mean product kinetic energy, applied to the mass as well.
- `sequence_count`, [`SequenceCount`](@ref). `A_p = A − n` with probability one: the
  prescription of [Tudora2022](@cite), of its applications to ²³⁹Pu(n_th,f)
  [TudoraFujio2024](@cite) and ²³³U(n_th,f) [Tudora2025](@cite), and of the code submitted with
  my thesis.

The default departs from the published count, and does so deliberately. Subtracting all `n`
removes `⟨n⟩ ≈ 2ν − 1` neutrons per fragment on average: 2.50 against `ν = 1.74` in the two
²⁵²Cf(sf) runs of my thesis, so that the product yields of one run imply 44 % more neutrons per
fission than its own multiplicity. Without fine-grained excitation-energy input the sequence average
stands in for the distribution, and a chain that can emit `n` neutrons emits `(n+1)/2` of them on
average; a product distribution that always removes `n` is that of the hottest ending alone. The
published count stays one key away, and with it the test suite reproduces the post-neutron mass
yields of both thesis runs to `1e-14`.

The yields follow from the fragment yield `Y(A, Z, TKE)` and the product probabilities `P`,

```
Y_p(Z, A_p) = 𝒩 Σ_{A, TKE} Y(A, Z, TKE) P(A_p | A, Z, TKE),
Y_p(A_p) = Σ_Z Y_p(Z, A_p),    Y_p(N_p) = Σ_{A_p − Z = N_p} Y_p(Z, A_p),
```

the sums running over the configurations solved. `𝒩` makes `Σ Y_p(A_p)` equal to the total of
the input yield, `[yield] total`, 200 for both fragments of every fission: the solved
configurations stand for the whole distribution, as they do in the average multiplicity. Emission
leaves the charge alone, so `Y_p(Z)` is the `Y(Z)` of the solved configurations and is not
written. At `A₀/2` a configuration is the heavy fragment of one labelling of a split and the light
fragment of the mirror labelling, and the self-complementary nuclide `(A₀/2, Z₀/2)` is both
fragments of one; each configuration enters once, through the chain its multiplicity is taken
from, as the input yield counts it once.

The kinetic energies are those of eqs. (6)–(9) of [Tudora2022](@cite),

```
E_K = TKE (A₀ − A)/A₀,    E_K_p = E_K A_p/A,    TKE_p = E_K_p,L + E_K_p,H,
```

emission isotropic in the fragment frame leaving the fragment velocity unchanged on average.
`E_K_p(A_p)` and `E_K_p(Z)` are means over `Y P`. `TKE_p` is taken over the two products of a pair,
independent of each other, with the joint probability `P_L P_H` and the pair weight that halves a
mass-symmetric split, as in the spectrum; it is written against the heavy product mass. Under the
sequence count these are eqs. (8) and (9) of the paper; under the default the mean of `E_K_p` over
the endings of a chain is its eq. (12).

A run writes `Y_vs_A_p.dat`, the independent yields `Y_vs_Z_A_p.dat`, `Y_vs_N_p.dat`,
`E_K_p_vs_A_p.dat`, `E_K_p_vs_Z.dat` and `TKE_p_vs_A_p.dat`, and adds `mean_product_mass_loss`,
`⟨A − A_p⟩` per fragment, and `mean_TKE_p` to `totals.dat`.

```@autodocs
Modules = [DeterministicSequentialEmission]
Pages = ["products.jl"]
```
