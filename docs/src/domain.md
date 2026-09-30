# The fragmentation domain

Everything this package takes from the fission fragments before emission begins lives in
[FissionFragmentsDomain.jl](https://github.com/PaulGoG/FissionFragmentsDomain.jl), which the
temperature-ratio extraction shares:

- nuclides, fissioning systems and their token;
- the tabular input contract and the shipped AME2020 mass table;
- the isobaric charge distribution, evaluated tables or Wahl's `Zₚ` model, and the fragmentation
  domain built from it;
- level density parameters and the relation between the temperature ratio and the
  excitation-energy partition;
- the energy balance and fragment kinematics;
- fragment yields and averages over them;
- the curves of the heavy-fragment mass and the temperature-ratio run record.

Its documentation is at <https://PaulGoG.github.io/FissionFragmentsDomain.jl/stable/>.

A configuration reaches it through three keys:
- `mass_excess_file = "ame2020"` takes the shipped table;
- `charge_distribution_file = "wahl"` takes Wahl's model at the best layer the system reaches;
- `[partition]` takes its `R_T(A_H)` curve through the run record of a temperature-ratio
  extraction.
