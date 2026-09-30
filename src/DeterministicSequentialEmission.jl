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

using DrWatson: savename
using FissionFragmentsDomain:
    AME2020_MASS_EXCESS_FILE,
    BackShiftedFermiGas,
    CHANNELS,
    ChargeDistribution,
    ChargeModel,
    DEFAULT_KINETIC_ENERGY_DISPERSION,
    DEFAULT_YIELD_TOTAL,
    EXCITATION_RATIO_CURVE_SPEC,
    FissioningSystem,
    FragmentQuantity,
    FragmentYield,
    Fragmentation,
    GILBERT_CAMERON_SHELL_CORRECTION_FILE,
    GilbertCameron,
    LEVEL_DENSITY_MODELS,
    LevelDensityModel,
    MassExcessTable,
    Nuclide,
    SegmentedCurve,
    TableSpec,
    TemperatureRatioManifest,
    WahlSystematics,
    ZpModel,
    average,
    build_fragment_yield,
    charge_dispersion,
    charge_model,
    charge_model_label,
    charge_polarization,
    column,
    coverage,
    curve_labels,
    effective_charge_distribution,
    energetics,
    factorized_yield,
    fragment_energy_per_nucleon,
    fragmentation_domain,
    heavy_excitation_fraction,
    integer_column,
    is_self_complementary,
    is_wahl_applicable,
    level_density_label,
    level_density_parameter,
    mean_charge_distribution,
    neutron_induced_fission,
    neutron_separation_energy,
    read_charge_distribution,
    read_delimited_table,
    read_kinetic_energy_dispersion,
    read_mass_energy_yield,
    read_mass_excess_table,
    read_segmented_curve,
    read_shell_correction_table,
    read_temperature_ratio_manifest,
    split_weight,
    spontaneous_fission,
    staged_manifest,
    system_label,
    system_record,
    temperature_ratio_path,
    total_average,
    uniform_kinetic_energy_dispersion,
    wahl_1988,
    weights_from,
    yield_over
using Measurements: value
using SpecialFunctions: gamma_inc
using TOML: TOML

export PartitionMethod,
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
export centre_of_mass_spectrum,
    laboratory_spectrum,
    chain_spectrum,
    chain_centre_of_mass_spectrum,
    SpectrumWeighting,
    MultiplicityWeighting,
    ExcitationWeighting,
    pair_weights,
    SpectrumAverage,
    PerFissionAverage,
    PerNeutronAverage,
    total_weight
export EmissionSequence,
    EmissionChain,
    emission_chain,
    is_solved,
    neutron_multiplicity,
    post_neutron_mass,
    sequence_average
export Configuration, load_configuration
export RunResult,
    run_pipeline,
    write_results,
    run_parameters,
    run_identifier,
    run_metadata,
    RUN_IDENTIFIER_ABBREVIATIONS

"""
    PACKAGE_ROOT

The directory the package was loaded from. The default input location is resolved relative to
it, so a configuration loads unchanged from any working directory. No output is written
relative to it; the output location is chosen by the caller, `scripts/run.jl`.
"""
const PACKAGE_ROOT = dirname(@__DIR__)

include("partition.jl")
include("evaporation.jl")
include("sequences.jl")
include("spectrum.jl")
include("configuration.jl")
include("pipeline.jl")

end
