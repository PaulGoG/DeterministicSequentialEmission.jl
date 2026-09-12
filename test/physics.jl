# The configuration of the archived reference run, and helpers shared by the test files.

const CF252 = spontaneous_fission(Nuclide(98, 252))
const CF252_HEAVY_MASSES = 126:174
const CF252_TKE = 130.0:2.0:230.0
const CF252_CHARGES_PER_MASS = 5

mass_table() = read_mass_table(joinpath(INPUT, "AME2020.ANA"))

cf252_polarization() =
    read_charge_polarization(joinpath(INPUT, "Cf-252", "DeltaZ_rms_A.CF52"))

cf252_domain() = fragmentation_domain(
    CF252,
    cf252_polarization(),
    CF252_HEAVY_MASSES;
    charges_per_mass = CF252_CHARGES_PER_MASS,
)

"""
    read_reference(name, columns) -> Dict

Read a reference file into a dictionary keyed by its leading integer columns.
"""
function read_reference(name::AbstractString, columns::Vector{Symbol}, keys::Int)
    table = read_table(joinpath(REFERENCE, name), TableSpec(columns; skip = 1))
    key_columns = [column(table, columns[i], Int) for i in 1:keys]
    values = column(table, columns[end], Float64)
    reference = Dict{NTuple{keys, Int}, Float64}()
    for row in eachindex(values)
        reference[ntuple(i -> key_columns[i][row], keys)] = values[row]
    end
    return reference
end
