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

using Dates: Dates
using Measurements: Measurement, measurement, value
using SpecialFunctions: gamma_inc
using TOML: TOML

export Nuclide, NEUTRON, neutron_number, element_symbol
export FissioningSystem,
    CHANNELS,
    CHANNEL_REACTION,
    REACTIONS,
    spontaneous_fission,
    neutron_induced_fission,
    is_spontaneous,
    reaction,
    symmetric_mass,
    system_label,
    system_notation
export TableSpec,
    DelimitedTable, read_delimited_table, column, column_or, integer_column, has_column
export MassExcessTable,
    MASS_EXCESS_SPEC,
    HEADERLESS_MASS_EXCESS_SPEC,
    read_mass_excess_table,
    mass_excess,
    separation_energy,
    neutron_separation_energy,
    q_value,
    complementary_fragment
export ChargeDistribution,
    CHARGE_DISTRIBUTION_SPEC,
    DEFAULT_CHARGE_POLARIZATION,
    DEFAULT_CHARGE_DISPERSION,
    read_charge_distribution,
    mean_charge_distribution,
    charge_polarization,
    charge_dispersion,
    unchanged_charge_distribution,
    most_probable_charge,
    charge_probability,
    charge_numbers
export WahlSystematics,
    WAHL_LOW_ENERGY,
    WAHL_VALIDITY,
    wahl_polarization,
    wahl_dispersion,
    wahl_charge_distribution
export Fragmentation,
    FragmentationDomain,
    fragmentation_domain,
    fragments,
    distinct_splits,
    mass_range,
    is_symmetric
export Energetics,
    Exclusion,
    energetics,
    sweep_energetics,
    exclusion_summary,
    compound_nucleus_excitation,
    total_excitation_energy,
    kinetic_energy
export LevelDensityModel,
    LEVEL_DENSITY_MODELS,
    ShellCorrectionTable,
    read_shell_correction_table,
    GilbertCameron,
    is_deformed_gilbert_cameron,
    BackShiftedFermiGas,
    level_density_parameter,
    shell_correction,
    liquid_drop_energy
export ManifestCurve,
    ManifestSystem,
    TemperatureRatioManifest,
    MANIFEST_CURVE_KINDS,
    MANIFEST_ORDINATE,
    MANIFEST_ABSCISSA,
    read_temperature_ratio_manifest,
    staged_manifest,
    curve_labels,
    manifest_curve,
    temperature_ratio_path
export SegmentedCurve,
    TEMPERATURE_RATIO_CURVE_SPEC,
    EXCITATION_RATIO_CURVE_SPEC,
    read_segmented_curve,
    PartitionMethod,
    ExcitationEnergyRatio,
    TemperatureRatio,
    ScissionModelling,
    segment_coverage,
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
export FragmentQuantity, average, total_average, weights_from
export MassEnergyYield,
    MASS_ENERGY_YIELD_SPEC,
    MASS_YIELD_SPEC,
    MEAN_KINETIC_ENERGY_SPEC,
    read_mass_energy_yield,
    factorized_yield,
    fragment_charge_probabilities,
    FragmentYield,
    DEFAULT_YIELD_TOTAL,
    build_fragment_yield,
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
export Configuration, load_configuration
export RunResult, run_pipeline, write_results, run_identifier, RUN_IDENTIFIER_ABBREVIATIONS

"""
    PACKAGE_ROOT

The directory the package was loaded from. Default input and output locations are taken relative
to it, so a pipeline runs unchanged from any working directory and on any machine.
"""
const PACKAGE_ROOT = dirname(@__DIR__)

include("nuclides.jl")
include("tables.jl")
include("masses.jl")
include("charge.jl")
include("domain.jl")
include("leveldensity.jl")
include("energetics.jl")
include("wahl.jl")
include("partition.jl")
include("manifest.jl")
include("evaporation.jl")
include("sequences.jl")
include("yields.jl")
include("averaging.jl")
include("spectrum.jl")
include("configuration.jl")
include("pipeline.jl")

end
