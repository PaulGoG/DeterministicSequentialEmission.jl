"""
    WahlReactionParameters

`Zₚ`-model parameters determined for one fission reaction by least squares, rather than estimated
for it from a trend across reactions.

Fields are Table A of A. C. Wahl, *Atomic Data and Nuclear Data Tables* **39**, 1 (1988), p. 7,
which fits the model to the evaluated fractional independent yields of that reaction. `assumed`
names the parameters the evaluation could not determine and fixed by assumption, which the table
marks with parentheses; they carry no uncertainty and should not be read as measurements.

`channels` are the entrance channels the fit covers, which is what keeps a reaction reaching the
same compound nucleus at a different excitation energy from claiming these parameters.

`shell` is where the narrow dispersion `σ_Z(50)` applies and `symmetric` where the even-odd
factors are set to one, both as footnotes a and c of that table give them. Each covers the `Z = 50`
closure and its complement by mass, so both legs of the range appear where the table gives two.
"""
struct WahlReactionParameters
    label::String
    channels::Vector{String}
    Z_F::Int
    A_F::Int
    ΔZ140::Float64
    ΔZSL::Float64
    σZ::Float64
    σZ50::Float64
    F_Z::Float64
    F_N::Float64
    ΔA_Z::Float64
    ΔZmax::Float64
    shell::Vector{UnitRange{Int}}
    symmetric::UnitRange{Int}
    assumed::Vector{Symbol}

    # Spelled out rather than left to the default constructor, so the fields are built from
    # exactly these types and inference never has to consider a conversion that cannot happen.
    function WahlReactionParameters(
        label::String,
        channels::Vector{String},
        Z_F::Int,
        A_F::Int,
        ΔZ140::Float64,
        ΔZSL::Float64,
        σZ::Float64,
        σZ50::Float64,
        F_Z::Float64,
        F_N::Float64,
        ΔA_Z::Float64,
        ΔZmax::Float64,
        shell::Vector{UnitRange{Int}},
        symmetric::UnitRange{Int},
        assumed::Vector{Symbol},
    )
        return new(
            label,
            channels,
            Z_F,
            A_F,
            ΔZ140,
            ΔZSL,
            σZ,
            σZ50,
            F_Z,
            F_N,
            ΔA_Z,
            ΔZmax,
            shell,
            symmetric,
            assumed,
        )
    end
end

function Base.show(io::IO, parameters::WahlReactionParameters)
    print(
        io,
        "WahlReactionParameters(",
        parameters.label,
        ", Z_F = ",
        parameters.Z_F,
        ", A_F = ",
        parameters.A_F,
        ")",
    )
    return nothing
end

"""
    WAHL_PER_REACTION

Table A of the 1988 evaluation, keyed by the charge and mass of the fissioning nucleus.

The four reactions it covers are the four that evaluation was made for. Values with their
uncertainties, as published:

| Parameter | U235T | U233T | PU239T | CF252S |
|---|---|---|---|---|
| `ΔZ(A' = 140)` | −0.511(5) | −0.519(8) | −0.544(9) | −0.420(20) |
| `∂ΔZ/∂A'` | −0.008(1) | −0.015(2) | −0.015(2) | −0.015(6) |
| `σ_Z` | 0.531(4) | 0.555(6) | 0.564(6) | 0.589(13) |
| `σ_Z(50)` | 0.33(2) | 0.36(3) | (0.35) | (0.35) |
| `F_Z` | 1.27(1) | 1.27(2) | 1.14(2) | 1.05(4) |
| `F_N` | 1.07(1) | 1.07(2) | 1.05(2) | (1.00) |
| `ΔA'_Z` | 0.9(2) | (0.9) | (0.9) | (0.7) |
| `ΔZ_max` | 0.7(1) | (0.7) | (0.7) | (0.5) |

Parenthesised values were assumed rather than determined, and are listed in each entry's
`assumed` field.

Note that `σ_Z` here is a single average over the mass range rather than a sloped function: the
1988 evaluation takes `σ_Z(A) = σ_Z` outside the shell region and `σ_Z(50)` inside it. It is the
model's Gaussian width, so [`wahl_charge_distribution`](@ref) puts it on the lattice footing of
eq. (8) before returning a [`ChargeDistribution`](@ref).
"""
const WAHL_PER_REACTION = Dict(
    (92, 236) => WahlReactionParameters(
        "U235T",
        ["nth", "nres"],
        92,
        236,
        -0.511,
        -0.008,
        0.531,
        0.33,
        1.27,
        1.07,
        0.9,
        0.7,
        [105:109, 125:129],
        105:129,
        Symbol[],
    ),
    (92, 234) => WahlReactionParameters(
        "U233T",
        ["nth", "nres"],
        92,
        234,
        -0.519,
        -0.015,
        0.555,
        0.36,
        1.27,
        1.07,
        0.9,
        0.7,
        [104:108, 123:128],
        104:128,
        [:ΔA_Z, :ΔZmax],
    ),
    (94, 240) => WahlReactionParameters(
        "PU239T",
        ["nth", "nres"],
        94,
        240,
        -0.544,
        -0.015,
        0.564,
        0.35,
        1.14,
        1.05,
        0.9,
        0.7,
        [109:113, 124:129],
        109:129,
        [:σZ50, :ΔA_Z, :ΔZmax],
    ),
    (98, 252) => WahlReactionParameters(
        "CF252S",
        ["sf"],
        98,
        252,
        -0.420,
        -0.015,
        0.589,
        0.35,
        1.05,
        1.00,
        0.7,
        0.5,
        [119:129],
        119:129,
        [:σZ50, :F_N, :ΔA_Z, :ΔZmax],
    ),
)

"""
    wahl_reaction_parameters(system) -> Union{WahlReactionParameters,Nothing}

The per-reaction `Zₚ` parameters for a fissioning system, or `nothing` where the 1988 evaluation
covers no such reaction.

Matched on the fissioning nucleus **and** the entrance channel. The nucleus alone is not enough:
²⁴⁰Pu(sf) and ²³⁹Pu(n_th,f) reach the same compound nucleus 6.53 MeV apart, and over that gap
eq. (17) moves `σ_Z` by 0.071 — twelve times the uncertainty Table A quotes on the fitted value.
A reaction that far from the one a fit was made for is better served by the systematics, which at
least takes the excitation energy as an argument.

The channels each entry admits are the ones its fit covers. `nres` joins `nth` because the two
differ by the incident energy alone: at the 580 eV resonance ²³⁵U's compound excitation is
6.546097 MeV against 6.545517 MeV at thermal, which moves no parameter in the sixth decimal.
"""
function wahl_reaction_parameters(system::FissioningSystem)
    parameters = get(WAHL_PER_REACTION, (system.compound.Z, system.compound.A), nothing)
    parameters === nothing && return nothing
    return system.channel in parameters.channels ? parameters : nothing
end

"""
    WahlSystematics(parameters, PE) -> WahlSystematics

Carry per-reaction parameters through the region functions of the `Zₚ` model.

The regions and their boundaries are those of LA-13928 — the later and fuller statement of the
same model — while the parameters within them are the ones the 1988 evaluation determined for
this reaction. Parameters that Table A does not carry, the near-symmetry slope `SL50` and the two
wing slopes, are taken from the systematics, which is the only source for them.

`σ_Z` has no slope here: the 1988 evaluation uses one average width outside the shell region, and
footnotes a and c of Table A say where that region and the near-symmetry region are, so those
ranges are carried through and take precedence over eq. (10)'s boundaries.

The footnotes place their ranges in fission-product mass while this package indexes pre-neutron
fragment mass. The two differ by `ν(A')`, which at the `Z = 50` closure these ranges sit on is at
the minimum of its sawtooth — 0.2 to 0.7 neutrons across the four reactions — so the difference is
below the one-mass-unit resolution of the ranges themselves and no shift is applied. Shifting them
by one unit changes the mean absolute difference of `σ_Z` from the four evaluated tables by 6 %,
which is the size of that ambiguity; ignoring the footnotes altogether costs 27 %, and for ²⁵²Cf a
factor of 3.4.
"""
function WahlSystematics(parameters::WahlReactionParameters, PE::Real)
    systematics = WahlSystematics(parameters.Z_F, parameters.A_F, PE)
    A_max = systematics.Bb
    B4 =
        (
            parameters.ΔZmax - parameters.ΔZ140 +
            A_max * systematics.SL50 +
            140 * parameters.ΔZSL
        ) / (systematics.SL50 + parameters.ΔZSL)
    return WahlSystematics(
        parameters.Z_F,
        parameters.A_F,
        Float64(PE),
        parameters.σZ,                 # σZ140: one average width, no slope
        parameters.ΔZ140,
        0.0,                           # σZSL
        parameters.ΔZSL,
        systematics.SL50,              # not in Table A
        parameters.σZ50,
        parameters.ΔZmax,
        systematics.σZSLW,             # not in Table A
        systematics.ΔZSLW,             # not in Table A
        parameters.F_Z,
        parameters.F_N,
        systematics.FZSLW,             # not in Table A
        systematics.FNSLW,             # not in Table A
        systematics.B2,
        B4,
        systematics.B5,
        systematics.B6,
        systematics.Ba,
        systematics.Bb,
        parameters.shell,
        parameters.symmetric,
    )
end

"""
    build_charge_distribution(table, system, masses; even_odd, charges) -> ChargeDistribution

The best isobaric charge distribution available for a fissioning system, resolved in order.

Three layers, each falling through to the next where it cannot answer:

1. **Model with this reaction's own data.** [`WAHL_PER_REACTION`](@ref) carries the `Zₚ`
   parameters the 1988 evaluation determined for this reaction by least squares against its
   measured fractional independent yields. Best where it exists; it exists for four reactions.
2. **The systematics.** LA-13928 eq. (17) estimates the same parameters from a trend across
   reactions, for any nucleus inside [`WAHL_VALIDITY`](@ref). Fig. 18 of that report puts the
   cost of the estimate at a reduced `χ²` of 7.9 against 2.9 for a per-reaction fit.
3. **The conventional means.** `ΔZ = −0.5`, `σ_Z = 0.6` at every mass, which carry no mass
   dependence at all. They are the yield-weighted averages of the tables the first two layers
   approximate, which is what justifies them and also what limits them.

The layer used is recorded in the returned distribution's `source`, with its citation, so a run's
metadata says which one answered rather than leaving it to be inferred.

An evaluated table read from a file is better than any of these and is not a layer here: supply
one through `charge_distribution_file` and this function is not consulted.
"""
function build_charge_distribution(
    table::MassExcessTable,
    system::FissioningSystem,
    masses;
    even_odd::Bool = false,
    charges::Integer = 11,
)
    excitation = compound_nucleus_excitation(table, system)
    parameters = wahl_reaction_parameters(system)
    if parameters !== nothing && excitation !== nothing
        PE = value(excitation)
        if 0 <= PE <= WAHL_VALIDITY.excitation
            return wahl_charge_distribution(
                WahlSystematics(parameters, PE),
                collect(masses);
                even_odd = even_odd,
                charges = charges,
                provenance = "Zₚ model, per-reaction parameters for $(parameters.label), " *
                             "Wahl, At. Data Nucl. Data Tables 39, 1 (1988), Table A",
            )
        end
    end
    if is_wahl_applicable(table, system)
        distribution = wahl_charge_distribution(
            table,
            system,
            collect(masses);
            even_odd = even_odd,
            charges = charges,
        )
        distribution === nothing || return distribution
    end
    return mean_charge_distribution()
end
