# The configuration of the archived reference run, and helpers shared by the test files.

const CF252 = spontaneous_fission(Nuclide(98, 252))
const CF252_HEAVY_MASSES = 126:174
const CF252_TKE = 130.0:2.0:230.0
const CF252_CHARGES_PER_MASS = 5

# The atomic mass evaluation ships with the package, so everything pinned against the archived
# run — Q, the level density parameter, both TXE partitions, the emission sequences — runs on a
# bare clone. The charge distribution systematics does not ship, because it is a third-party
# table this package has no licence to redistribute; the tests that need it are skipped when it
# is absent, exactly as the companion packages do. `data/README.md` says where to obtain it.
const MASS_EXCESS_FILE = joinpath(DATA, "reference", "mass_excess_ame2020.dat")
const CF252_CHARGE_DISTRIBUTION_FILE =
    joinpath(DATA, "Cf252_sf", "charge_distribution_vs_A.dat")
const U235_CHARGE_DISTRIBUTION_FILE =
    joinpath(DATA, "U235_nth", "charge_distribution_vs_A.dat")

const SHELL_CORRECTION_FILE =
    joinpath(DATA, "reference", "shell_corrections_gilbert_cameron.dat")

const CHARGE_DISTRIBUTION_AVAILABLE = isfile(CF252_CHARGE_DISTRIBUTION_FILE)
const U235_POLARIZATION_AVAILABLE = isfile(U235_CHARGE_DISTRIBUTION_FILE)
const SHELL_CORRECTIONS_AVAILABLE = isfile(SHELL_CORRECTION_FILE)

const CF252_YIELD_FILE = joinpath(DATA, "Cf252_sf", "Y_vs_A_TKE.dat")
const YIELD_AVAILABLE = isfile(CF252_YIELD_FILE)

CHARGE_DISTRIBUTION_AVAILABLE || @warn """
    charge distribution systematics not present; the tests that build a fragmentation domain \
    from it are skipped. The reference tests are unaffected — they take their fragmentations \
    from the reference files. See data/README.md.""" file = CF252_CHARGE_DISTRIBUTION_FILE

mass_table() = read_mass_excess_table(MASS_EXCESS_FILE)

cf252_charge_distribution() = read_charge_distribution(CF252_CHARGE_DISTRIBUTION_FILE)

cf252_domain() = fragmentation_domain(
    CF252,
    cf252_charge_distribution(),
    CF252_HEAVY_MASSES;
    charges_per_mass = CF252_CHARGES_PER_MASS,
)

"""
    fragmentation_at(A_H, Z_H) -> Fragmentation

The ²⁵²Cf fragmentation with this heavy fragment, built without the charge distribution.

The reference tests pin `Q`, `a`, `E*_L`, `E*_H` and the emission sequences, none of which
depends on the isobaric charge probability — it is a yield weight, and these are per
configuration. Taking the fragmentation straight from `(A_H, Z_H)` keeps those tests running
without the unshipped systematics, and isolates them from the domain construction, which is
tested separately.
"""
function fragmentation_at(A_H::Integer, Z_H::Integer)
    heavy = Nuclide(Z_H, A_H)
    return Fragmentation(heavy, complementary_fragment(CF252, heavy), 1.0)
end

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
