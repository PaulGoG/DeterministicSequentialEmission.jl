"""
    LevelDensityModel

Systematics for the energy-independent level density parameter `a(A,Z)`.

The sequential equations close on `Ē_r = a T²`, the Fermi-gas relation, which admits no
explicit energy dependence in `a`. That is why phenomenological systematics are used here
rather than a model such as the superfluid one, whose `a(E*)` would make the closure
circular.

Values are plain `Float64` rather than `Measurement`s. `a` is consumed by the TXE partition
and by the sequential solver, both downstream of the point where uncertainty propagation
stops, so attaching a mass-excess uncertainty here would carry a number that is never
reported and would be misleading if it were.
"""
abstract type LevelDensityModel end

"""
    ShellCorrections

Gilbert-Cameron shell corrections `S(Z)` and `S(N)`, indexed by nucleon number.
"""
struct ShellCorrections
    proton::Dict{Int, Float64}
    neutron::Dict{Int, Float64}
    source::String

    function ShellCorrections(
        proton::Dict{Int, Float64},
        neutron::Dict{Int, Float64},
        source::AbstractString,
    )
        return new(proton, neutron, String(source))
    end
end

"""
    SHELL_CORRECTION_SPEC

Layout of the shell-correction file: one header row, then `n S(N) S(Z)`.

**The neutron correction comes first.** The file's own header reads `n S(N) S(Z)`, and the
earlier coursework implementation of this model read it in that order. The prototype ported
for this package declares `["n", "S_Z", "S_N"]`, i.e. the two swapped. Because the
Gilbert-Cameron expression looks `S(Z)` up at `n = Z` and `S(N)` at `n = N`, a swap is not
absorbed by the sum: it evaluates the neutron correction at the proton number and vice versa.

The archived reference run used the back-shifted Fermi gas, so no published result depends on
which order was used, but a Gilbert-Cameron run of the prototype does.
"""
const SHELL_CORRECTION_SPEC = TableSpec([:n, :neutron, :proton]; skip = 1)

"""
    read_shell_corrections(path; spec = SHELL_CORRECTION_SPEC) -> ShellCorrections

Read the Gilbert-Cameron shell corrections.
"""
function read_shell_corrections(
    path::AbstractString;
    spec::TableSpec = SHELL_CORRECTION_SPEC,
)
    table = read_table(path, spec)
    nucleons = column(table, :n, Int)
    neutron = column(table, :neutron, Float64)
    proton = column(table, :proton, Float64)

    protons = Dict{Int, Float64}()
    neutrons = Dict{Int, Float64}()
    for row in eachindex(nucleons)
        protons[nucleons[row]] = proton[row]
        neutrons[nucleons[row]] = neutron[row]
    end
    return ShellCorrections(protons, neutrons, String(path))
end

function Base.show(io::IO, corrections::ShellCorrections)
    print(
        io,
        "ShellCorrections(",
        length(corrections.proton),
        " entries from ",
        basename(corrections.source),
        ")",
    )
    return nothing
end

"""
    GilbertCameron(corrections)

Gilbert-Cameron systematics with shell corrections,

```
a(A,Z) = A [0.00917 (S(Z) + S(N)) + 0.142]
```

with `S(Z)` taken at the proton number and `S(N)` at the neutron number.
"""
struct GilbertCameron <: LevelDensityModel
    corrections::ShellCorrections
end

const GILBERT_CAMERON_SHELL_WEIGHT = 0.00917
const GILBERT_CAMERON_OFFSET = 0.142

"""
    BackShiftedFermiGas(table)

Von Egidy-Bucurescu systematics for the back-shifted Fermi gas,

```
a(A,Z) = (p₁ + p₂ δW) A^p₃
```

where the shell correction `δW = δW₀ + P_d` is the difference between a liquid-drop mass and
the experimental one, plus a deuteron pairing term built from the mass excesses of the
neighbours `(A+2, Z+1)` and `(A−2, Z−1)`.

Preferred over Gilbert-Cameron for fission fragments: across the fragment range it stays close
to the superfluid model at the excitation energies fragments actually have, departing mainly
near `A ≈ 130` where the `N = 82`, `Z = 50` shells close.

Verified against the source. The systematics is that of the 2009 update, eq. (19):

```
a = (p₁ + p₂ S′) A^p₃,   p₁ = 0.199(7), p₂ = 0.0096(4), p₃ = 0.869(7)
```

The form is additive in `S′`. A variant `p₁(1 + p₂ S′) A^p₃` circulates and is **not** what the
paper gives; the two differ by a factor `p₁` in the weight on `S′`.

The pieces, each read from the papers rather than inferred:

- `S(Z,N) = M_exp − M_LD` (2005, eq. 7), which in binding energies is `E_b^LD − E_b^exp`.
- `M_LD` from Pearson's Weizsäcker-type formula (2005, eq. 9), with `a_vol = −15.65`,
  `a_sf = 17.63`, `a_sym = 27.72`, `a_ss = −25.60` MeV and `r₀ = 1.233` fm, so the Coulomb
  coefficient is `3e²/5r₀ = 0.864/1.233` MeV with `e² = 1.44` MeV fm.
- `S′ = S + 0.5 P_a′` (2009, eq. 19) with
  `P_a′ = ½[M(A+2,Z+1) − 2M(A,Z) + M(A−2,Z−1)]` (2009, eq. 12), so the pairing enters as
  `¼[M(A+2,Z+1) − 2M(A,Z) + M(A−2,Z−1)]`.

Two traps avoided. `P_a′` carries no `(−1)^Z` alternation: it differs by `(−1)^(Z+1)` from the
`P_a` of the 2005 paper and of Audi's tables, as the 2009 paper states where it defines it, and
it is `P_a′` that enters eq. (19). And the 2005 paper uses a different functional form
altogether, `a/A = p₁ + p₂ S′ + p₃ A`, with a parity-branched pairing correction; this follows
the 2009 update throughout rather than mixing the two.

References held in the workspace: von Egidy & Bucurescu, Phys. Rev. C 72, 044311 (2005);
erratum Phys. Rev. C 73, 049901 (2006); Phys. Rev. C 80, 054310 (2009).

"""
struct BackShiftedFermiGas <: LevelDensityModel
    table::MassTable
    proton_excess::Float64
    neutron_excess::Float64

    function BackShiftedFermiGas(table::MassTable)
        proton = mass_excess(table, Nuclide(1, 1))
        neutron = mass_excess(table, NEUTRON)
        (proton === nothing || neutron === nothing) && throw(
            ArgumentError(
                "the back-shifted Fermi gas systematics needs the mass excesses of ¹H and the " *
                "neutron, and $(table.source) supplies at least one of them not at all",
            ),
        )
        return new(table, value(proton), value(neutron))
    end
end

# Liquid-drop coefficients entering δW₀, in MeV.
const LDM_VOLUME = 15.65
const LDM_SURFACE = 17.63
const LDM_COULOMB = 0.864 / 1.233
const LDM_SYMMETRY_VOLUME = 27.72
const LDM_SYMMETRY_SURFACE = 25.6

# Coefficients of the level density parameter itself.
const BSFG_OFFSET = 1.99e-1
const BSFG_SHELL_WEIGHT = 9.6e-3
const BSFG_MASS_EXPONENT = 8.69e-1

"""
    liquid_drop_energy(model, nuclide) -> Float64

Liquid-drop binding energy in MeV, with a mass-dependent symmetry term.
"""
function liquid_drop_energy(::BackShiftedFermiGas, nuclide::Nuclide)
    A = nuclide.A
    Z = nuclide.Z
    symmetry = A * (LDM_SYMMETRY_VOLUME - LDM_SYMMETRY_SURFACE * A^(-1 / 3))
    η = (A - 2Z) / A
    return LDM_VOLUME * A - LDM_SURFACE * A^(2 / 3) - LDM_COULOMB * Z^2 * A^(-1 / 3) -
           symmetry * η^2
end

"""
    shell_correction(model, nuclide) -> Union{Float64,Nothing}

`δW = δW₀ + P_d` in MeV, the difference between the liquid-drop and experimental binding
energies plus the deuteron pairing term

```
P_d = [Δ(A+2, Z+1) − 2Δ(A, Z) + Δ(A−2, Z−1)] / 4
```

`nothing` when a required mass excess is unavailable, which is what makes the parameter
undefined for the most exotic fragments of a deterministic domain.
"""
function shell_correction(model::BackShiftedFermiGas, nuclide::Nuclide)
    Δ = mass_excess(model.table, nuclide)
    Δ === nothing && return nothing
    Δ_up = mass_excess(model.table, Nuclide(nuclide.Z + 1, nuclide.A + 2))
    Δ_up === nothing && return nothing
    nuclide.Z >= 1 && nuclide.A >= 2 || return nothing
    Δ_down = mass_excess(model.table, Nuclide(nuclide.Z - 1, nuclide.A - 2))
    Δ_down === nothing && return nothing

    experimental =
        nuclide.Z * model.proton_excess + neutron_number(nuclide) * model.neutron_excess -
        value(Δ)
    δW₀ = liquid_drop_energy(model, nuclide) - experimental
    pairing = (value(Δ_up) - 2 * value(Δ) + value(Δ_down)) / 4
    return δW₀ + pairing
end

"""
    level_density_parameter(model, nuclide) -> Union{Float64,Nothing}

The energy-independent level density parameter `a(A,Z)` in MeV⁻¹, or `nothing` when the
systematics cannot be evaluated or returns a non-positive value.

A non-positive `a` is not a physical level density parameter and would give an imaginary
temperature through `Ē_r = a T²`, so it is reported as unavailable rather than propagated.

# Examples

```jldoctest
julia> table = read_mass_table(joinpath(pkgdir(DeterministicSequentialEmission), "input_data", "AME2020.ANA"));

julia> a = level_density_parameter(BackShiftedFermiGas(table), Nuclide(52, 134));

julia> round(a; digits = 3)
7.707
```
"""
function level_density_parameter(model::BackShiftedFermiGas, nuclide::Nuclide)
    δW = shell_correction(model, nuclide)
    δW === nothing && return nothing
    a = (BSFG_OFFSET + BSFG_SHELL_WEIGHT * δW) * nuclide.A^BSFG_MASS_EXPONENT
    return a > 0 ? a : nothing
end

function level_density_parameter(model::GilbertCameron, nuclide::Nuclide)
    S_proton = get(model.corrections.proton, nuclide.Z, nothing)
    S_proton === nothing && return nothing
    S_neutron = get(model.corrections.neutron, neutron_number(nuclide), nothing)
    S_neutron === nothing && return nothing
    a =
        nuclide.A *
        (GILBERT_CAMERON_SHELL_WEIGHT * (S_proton + S_neutron) + GILBERT_CAMERON_OFFSET)
    return a > 0 ? a : nothing
end
