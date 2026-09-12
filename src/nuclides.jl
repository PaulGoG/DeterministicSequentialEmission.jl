"""
    Nuclide(Z, A)

A nuclide identified by its proton number `Z` and mass number `A`.
"""
struct Nuclide
    Z::Int
    A::Int

    function Nuclide(Z::Integer, A::Integer)
        Z >= 0 || throw(ArgumentError("proton number must be non-negative, got Z = $Z"))
        A >= Z || throw(
            ArgumentError(
                "mass number must be at least the proton number, got A = $A, Z = $Z",
            ),
        )
        return new(Int(Z), Int(A))
    end
end

"""
    neutron_number(nuclide) -> Int

Neutron number `N = A - Z`.
"""
neutron_number(nuclide::Nuclide) = nuclide.A - nuclide.Z

Base.show(io::IO, nuclide::Nuclide) =
    print(io, "Nuclide(Z=", nuclide.Z, ", A=", nuclide.A, ")")

const NEUTRON = Nuclide(0, 1)

"""
    FissioningSystem(; compound, target, incident_energy, spontaneous)

The fissioning system, holding both the target and the compound nucleus explicitly.

For spontaneous fission the two coincide and `incident_energy` is zero. For neutron-induced
fission the compound nucleus is the target plus one neutron, and `incident_energy` is the
kinetic energy of the incident neutron in MeV.

Keeping both nuclides as separate fields is deliberate: the mass number that defines the
fragmentation domain is that of the compound nucleus, while validation of a requested fragment
range and the labelling of a fission case refer to the target. Deriving one from the other at
use sites is what makes the distinction easy to lose.

# Examples

```jldoctest
julia> spontaneous_fission(Nuclide(98, 252))
FissioningSystem(252Cf(SF))

julia> neutron_induced_fission(Nuclide(94, 239), 2.53e-8)
FissioningSystem(239Pu(n,f) at E = 2.53e-8 MeV)
```
"""
struct FissioningSystem
    compound::Nuclide
    target::Nuclide
    incident_energy::Float64
    spontaneous::Bool

    function FissioningSystem(
        compound::Nuclide,
        target::Nuclide,
        incident_energy::Real,
        spontaneous::Bool,
    )
        incident_energy >= 0 || throw(
            ArgumentError("incident energy must be non-negative, got $incident_energy MeV"),
        )
        if spontaneous
            compound == target || throw(
                ArgumentError(
                    "spontaneous fission requires the compound and target nuclides to coincide, " *
                    "got compound $compound and target $target",
                ),
            )
            iszero(incident_energy) || throw(
                ArgumentError(
                    "spontaneous fission requires zero incident energy, got $incident_energy MeV",
                ),
            )
        end
        return new(compound, target, Float64(incident_energy), spontaneous)
    end
end

"""
    spontaneous_fission(nuclide) -> FissioningSystem

The spontaneously fissioning system formed by `nuclide`.
"""
spontaneous_fission(nuclide::Nuclide) = FissioningSystem(nuclide, nuclide, 0.0, true)

"""
    neutron_induced_fission(target, incident_energy) -> FissioningSystem

The system formed by `target` capturing a neutron of kinetic energy `incident_energy` in MeV.
The compound nucleus carries one more mass unit than the target.
"""
function neutron_induced_fission(target::Nuclide, incident_energy::Real)
    compound = Nuclide(target.Z, target.A + 1)
    return FissioningSystem(compound, target, incident_energy, false)
end

"""
    symmetric_mass(system) -> Float64

The mass number of a symmetric split of the compound nucleus, `A₀ / 2`. Not generally an
integer; used to validate requested heavy-fragment ranges.
"""
symmetric_mass(system::FissioningSystem) = system.compound.A / 2

const ELEMENT_SYMBOLS = (
    "n",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
)

"""
    element_symbol(Z) -> String

Chemical symbol for proton number `Z`, or `"Z=\$Z"` beyond the tabulated range.
"""
function element_symbol(Z::Integer)
    return 0 <= Z <= length(ELEMENT_SYMBOLS) - 1 ? ELEMENT_SYMBOLS[Z + 1] : string("Z=", Z)
end

function Base.show(io::IO, system::FissioningSystem)
    label = string(system.target.A, element_symbol(system.target.Z))
    if system.spontaneous
        print(io, "FissioningSystem(", label, "(SF))")
    else
        print(
            io,
            "FissioningSystem(",
            label,
            "(n,f) at E = ",
            system.incident_energy,
            " MeV)",
        )
    end
    return nothing
end
