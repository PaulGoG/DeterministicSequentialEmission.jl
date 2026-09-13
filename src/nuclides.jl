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
    CHANNELS

The entrance channels a fissioning system may be declared with: spontaneous fission, and
thermal, resonance-region and fast neutron-induced fission.

The channel is spelled as the field spells it and not as a reaction code writes it. It is what
distinguishes a thermal run of a system from a resonance run of the same system, which `n,f`
alone cannot, and it is what names the system in a path — see [`system_label`](@ref).
"""
const CHANNELS = ("sf", "nth", "nres", "nfast")

"""
    CHANNEL_REACTION

The reaction code each entrance channel stands for. The code is derived rather than configured,
because a key that can only be redundant or wrong does not belong in a configuration file.
"""
const CHANNEL_REACTION =
    Dict("sf" => "0,f", "nth" => "n,f", "nres" => "n,f", "nfast" => "n,f")

"""
    REACTIONS

Neutrons absorbed by the target in each reaction, which is what derives the compound mass
`A₀ = A_target + REACTIONS[reaction]`.
"""
const REACTIONS = Dict("0,f" => 0, "n,f" => 1)

"""
    FissioningSystem(compound, target, incident_energy, channel)

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
FissioningSystem(252Cf(sf))

julia> neutron_induced_fission(Nuclide(94, 239), 2.53e-8, "nth")
FissioningSystem(239Pu(nth,f) at E = 2.53e-8 MeV)
```
"""
struct FissioningSystem
    compound::Nuclide
    target::Nuclide
    incident_energy::Float64
    channel::String

    function FissioningSystem(
        compound::Nuclide,
        target::Nuclide,
        incident_energy::Real,
        channel::AbstractString,
    )
        channel in CHANNELS || throw(ArgumentError("the entrance channel must be one of \
                 $(join(map(repr, CHANNELS), ", ")), got $(repr(channel))"))
        incident_energy >= 0 || throw(
            ArgumentError("incident energy must be non-negative, got $incident_energy MeV"),
        )
        if channel == "sf"
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
        return new(compound, target, Float64(incident_energy), String(channel))
    end
end

"""
    is_spontaneous(system) -> Bool

Whether the system fissions spontaneously, which is the `"sf"` entrance channel and no other.
"""
is_spontaneous(system::FissioningSystem) = system.channel == "sf"

"""
    reaction(system) -> String

The reaction code of the system's entrance channel, `"0,f"` or `"n,f"`.
"""
reaction(system::FissioningSystem) = CHANNEL_REACTION[system.channel]

"""
    spontaneous_fission(nuclide) -> FissioningSystem

The spontaneously fissioning system formed by `nuclide`.
"""
spontaneous_fission(nuclide::Nuclide) = FissioningSystem(nuclide, nuclide, 0.0, "sf")

"""
    neutron_induced_fission(target, incident_energy, channel) -> FissioningSystem

The system formed by `target` capturing a neutron of kinetic energy `incident_energy` in MeV,
on the entrance channel `channel` — one of `"nth"`, `"nres"` or `"nfast"`.

The compound nucleus carries `REACTIONS["n,f"]` more mass units than the target. The channel
still has to agree with the incident energy: the thermal, resonance and fast boundaries are
conventions rather than thresholds this package can enforce.
"""
function neutron_induced_fission(
    target::Nuclide,
    incident_energy::Real,
    channel::AbstractString,
)
    channel == "sf" && throw(
        ArgumentError("\"sf\" is spontaneous fission; use `spontaneous_fission` instead"),
    )
    compound = Nuclide(target.Z, target.A + REACTIONS["n,f"])
    return FissioningSystem(compound, target, incident_energy, channel)
end

"""
    system_label(system) -> String

The one token that names a fissioning system, `<ElementSymbol><A><_channel>`: `Cf252_sf`,
`U235_nth`, `U235_nres`. It is the directory a run writes to, the directory its input is read
from, and the name of the configuration file that declares it.

The typeset form for a figure or a caption is [`system_notation`](@ref); one name may not mean
both.

# Examples

```jldoctest
julia> system_label(spontaneous_fission(Nuclide(98, 252)))
"Cf252_sf"

julia> system_label(neutron_induced_fission(Nuclide(92, 235), 5.8013e-4, "nres"))
"U235_nres"
```
"""
system_label(system::FissioningSystem) =
    string(element_symbol(system.target.Z), system.target.A, "_", system.channel)

"""
    system_notation(system) -> String

The typeset form of a fissioning system for a figure or a caption: `²⁵²Cf(sf)`, `²³⁵U(nth,f)`.

The channel appears as the plain field spelling inside the reaction parentheses, because that is
how the fission literature writes it and it survives a terminal, a PDF and a grep equally.

# Examples

```jldoctest
julia> system_notation(spontaneous_fission(Nuclide(98, 252)))
"²⁵²Cf(sf)"
```
"""
function system_notation(system::FissioningSystem)
    superscripts = Dict(
        '0' => '⁰',
        '1' => '¹',
        '2' => '²',
        '3' => '³',
        '4' => '⁴',
        '5' => '⁵',
        '6' => '⁶',
        '7' => '⁷',
        '8' => '⁸',
        '9' => '⁹',
    )
    mass = map(digit -> superscripts[digit], string(system.target.A))
    inside = is_spontaneous(system) ? "sf" : string(system.channel, ",f")
    return string(mass, element_symbol(system.target.Z), "(", inside, ")")
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
    if is_spontaneous(system)
        print(io, "FissioningSystem(", label, "(sf))")
    else
        print(
            io,
            "FissioningSystem(",
            label,
            "(",
            system.channel,
            ",f) at E = ",
            system.incident_energy,
            " MeV)",
        )
    end
    return nothing
end
