"""
    DeterministicSequentialEmission

Prompt neutron emission from fission fragments, treated sequence by sequence.

The model solves the energy conservation equation of every emission sequence separately,

```
Ē_r^(k-1) − Sₙ^(k-1) − ⟨ε⟩_k = a_k T_k²,    k = 1 … n(A,Z,TKE)
```

with Weisskopf-Ewing evaporation spectra in the centre-of-mass system and fragment level
densities in the Fermi-gas regime. Because the number of sequences `n` is an integer at every
configuration, the post-neutron fragment mass `A_p = A − n` is an integer without rounding.
"""
module DeterministicSequentialEmission

using Measurements: Measurement, measurement, value
using SpecialFunctions: gamma_inc

export Nuclide, NEUTRON, neutron_number, element_symbol
export FissioningSystem, spontaneous_fission, neutron_induced_fission, symmetric_mass
export TableSpec, DataTable, read_table, column, column_or, integer_column, has_column
export MassTable,
    MASS_EXCESS_SPEC,
    LEGACY_MASS_EXCESS_SPEC,
    read_mass_table,
    mass_excess,
    separation_energy,
    neutron_separation_energy,
    q_value,
    complementary_fragment
export ChargePolarization,
    read_charge_polarization,
    mean_charge_polarization,
    polarization_at,
    width_at,
    unchanged_charge_distribution,
    most_probable_charge,
    charge_probability,
    charge_numbers
export Fragmentation,
    FragmentationDomain, fragmentation_domain, fragments, mass_range, is_symmetric
export Energetics,
    ExcludedConfiguration,
    energetics,
    sweep_energetics,
    exclusion_summary,
    compound_nucleus_excitation,
    total_excitation_energy,
    kinetic_energy
export LevelDensityModel,
    ShellCorrections,
    read_shell_corrections,
    GilbertCameron,
    is_deformed_gilbert_cameron,
    BackShiftedFermiGas,
    level_density_parameter,
    shell_correction,
    liquid_drop_energy
export JointedSegments,
    SEGMENT_SPEC,
    read_segments,
    PartitionMethod,
    ExcitationEnergyRatio,
    TemperatureRatio,
    ScissionModelling,
    read_extra_deformation,
    ExcitationEnergies,
    partition
export EvaporationCrossSection,
    ConstantCrossSection,
    VariableCrossSection,
    neutron_strength_function,
    alpha_parameter,
    mean_neutron_energy,
    fermi_gas_energy,
    residual_temperature,
    MINIMUM_TEMPERATURE
export ExperimentalYield,
    EXPERIMENTAL_YIELD_SPEC,
    read_experimental_yield,
    fragment_charge_probabilities,
    FragmentYields,
    DEFAULT_YIELD_TOTAL,
    build_fragment_yields,
    coverage,
    yield_over,
    mean_kinetic_energy_by_mass
export centre_of_mass_spectrum,
    fragment_energy_per_nucleon,
    laboratory_spectrum,
    chain_spectrum,
    chain_centre_of_mass_spectrum,
    SpectrumWeighting,
    MultiplicityWeighting,
    ExcitationWeighting,
    pair_weights
export EmissionSequence,
    EmissionChain,
    emission_chain,
    is_solved,
    neutron_multiplicity,
    post_neutron_mass,
    sequence_average

include("nuclides.jl")
include("tables.jl")
include("masses.jl")
include("charge.jl")
include("domain.jl")
include("leveldensity.jl")
include("energetics.jl")
include("partition.jl")
include("evaporation.jl")
include("sequences.jl")
include("yields.jl")
include("spectrum.jl")

end
