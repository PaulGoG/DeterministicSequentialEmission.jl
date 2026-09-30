```@meta
CurrentModule = DeterministicSequentialEmission
```

# The excitation-energy partition

How the total excitation energy of a fragmentation is divided between its fragments: by a
tabulated `E*_H/TXE (A_H)`, by a temperature ratio `R_T(A_H)` through the relation
`E*_H/TXE = 1/(1 + R_a R_T²)`, or by modelling at scission. The curves, the level density
parameters and the relation come from the domain package; see [The fragmentation domain](domain.md).

```@autodocs
Modules = [DeterministicSequentialEmission]
Pages = ["partition.jl"]
```
