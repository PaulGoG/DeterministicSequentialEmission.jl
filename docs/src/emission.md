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
