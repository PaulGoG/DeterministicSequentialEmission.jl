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
    nucleons = integer_column(table, :n)
    neutron = column(table, :neutron, Float64)
    proton = column(table, :proton, Float64)

    # `S(Z)` is tabulated to Z = 98 and `S(N)` to N = 150, but the file is rectangular: the
    # proton column is padded with an exact 0.00 from 99 on. No tabulated `S(Z)` is zero, so
    # the padding is unambiguous — and storing it would silently hand `GilbertCameron` a
    # correction of zero for Z ≥ 99 instead of letting it decline.
    last_proton = something(findlast(!iszero, proton), 0)

    protons = Dict{Int, Float64}()
    neutrons = Dict{Int, Float64}()
    for row in eachindex(nucleons)
        row <= last_proton && (protons[nucleons[row]] = proton[row])
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
a(A,Z) = A [0.00917 (S(Z) + S(N)) + c],   c = 0.142 undeformed, 0.120 deformed
```

with `S(Z)` taken at the proton number and `S(N)` at the neutron number, and the two shell
corrections read from the paper's Table III.

**The offset is not a universal constant**, and which branch applies is decided by
[`is_deformed_gilbert_cameron`](@ref). Gilbert & Cameron give two parallel correlation lines —
eq. (20), `a/A = 0.00917S + 0.142` for undeformed nuclei, and eq. (21), `a/A = 0.00917S + 0.120`
for deformed ones (p. 1457, restated p. 1475) — together with the regions each covers
(p. 1459). The shell weight `0.00917` is identical in both; only the intercept moves.

`deformed_branch` selects between the two, and **there is a real argument on each side**, so it
is a choice rather than a default with an obvious answer.

Taking both branches, as this does unless told otherwise, is what the paper prescribes, and it
matters: about a third of the yield-weighted ²⁵²Cf fragment population lies inside the first
deformed region, where the undeformed constant overestimates `a` by some 22 % — see
[`is_deformed_gilbert_cameron`](@ref). On those fragments the deformed branch also moves
Gilbert-Cameron towards the back-shifted Fermi gas rather than away from it.

Setting `deformed_branch = false` applies eq. (20) everywhere. That is the convention the
published results of this method were obtained under — Tudora, *EPJA* 54, 87 (2018) labels its
comparison "G-C systm. (spherical nuclei)" — and it is defensible on its own terms: the data
behind eq. (21) is thin, so the deformation boundary is poorly determined, which is part of why
the systematics was superseded by the back-shifted Fermi gas, which needs no such division. It is
also the setting that matches `FissionTemperatureRatio.jl`, which implements eq. (20) alone.

**Use `deformed_branch = false` when the level density parameter has to agree with that package**
— when an `R_T(A_H)` it produced is being inverted here, both sides must form `a_L/a_H` the same
way. In practice both default to the back-shifted Fermi gas, so the interface is unaffected
unless Gilbert-Cameron is selected on both sides.

Nothing archived is affected either way: both reference runs used the back-shifted Fermi gas.

Source: A. Gilbert and A. G. W. Cameron, *Canadian Journal of Physics* **43**, 1446 (1965);
eqs. (20) and (21) p. 1457, the deformation regions p. 1459, Table III pp. 1453-1455.
"""
struct GilbertCameron <: LevelDensityModel
    corrections::ShellCorrections
    deformed_branch::Bool

    function GilbertCameron(corrections::ShellCorrections; deformed_branch::Bool = true)
        return new(corrections, deformed_branch)
    end
end

const GILBERT_CAMERON_SHELL_WEIGHT = 0.00917
const GILBERT_CAMERON_OFFSET_UNDEFORMED = 0.142
const GILBERT_CAMERON_OFFSET_DEFORMED = 0.120

"""
    is_deformed_gilbert_cameron(nuclide) -> Bool

Whether Gilbert-Cameron treats this nuclide as deformed, and so takes the `0.120` branch.

The paper's own definition, p. 1459: "we shall consider the following to be regions of
deformation: `54 ≤ Z ≤ 78, 86 ≤ N ≤ 122`" and "`86 ≤ Z ≤ 122, 130 ≤ N ≤ 182`", with both
bounds of a pair required together. Its further prescription for light nuclei — `11 ≤ Z ≤ 19`
follows the deformed line too, `20 ≤ Z ≤ 29` the undeformed one (p. 1464) — is included for
completeness; it is outside the fragment range.

This is not a corner case for fission. Roughly a quarter of the yield-weighted fragment
population of ²³⁵U(n,f) falls inside the first region — the upper half of the heavy fragment
peak, `A ≳ 141` — where the undeformed constant overestimates `a` by about 23 %. Since
`Ē_r = a T²`, that is a ~10 % deficit in `T` and, through `⟨ε⟩ = 2T`, in the mean emitted
neutron energy, concentrated entirely on the heavy fragment.
"""
function is_deformed_gilbert_cameron(nuclide::Nuclide)
    Z = nuclide.Z
    N = neutron_number(nuclide)
    (54 <= Z <= 78 && 86 <= N <= 122) && return true
    (86 <= Z <= 122 && 130 <= N <= 182) && return true
    return 11 <= Z <= 19
end

"""
    BackShiftedFermiGas(table)

Von Egidy-Bucurescu systematics for the back-shifted Fermi gas,

```
a(A,Z) = (p₁ + p₂ δW) A^p₃
```

where the shell correction `δW = δW₀ + ½P_a′` is the difference between the experimental mass
and a liquid-drop one, plus a pairing term built from the mass excesses of the neighbours
`(A+2, Z+1)` and `(A−2, Z−1)`. The subtraction runs `M_exp − M_LD`, in that order; the 2009
paper's own prose states it the other way round on p. 054310-3 while its eq. (7) gives this
order, and eq. (7) is what the fit was made against.

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

Only the parameter `a` is taken from this systematics, and that is a partial use of it. The
papers fit `a` jointly with the back-shift `E1 = −0.381 + 0.5 P_a′` (2009, eq. 20) for the level
density `ρ(U) ∝ exp(2√(a(U − E1)))`, whose thermodynamic closure is `U − E1 = a T²`. The
sequential equations here close on the un-shifted `Ē_r = a T²`, so `a` enters as an effective
parameter rather than as half of the fitted pair it belongs to.

**This is the source model's convention, not a porting slip**, and the papers that adopt the
systematics are silent about it rather than arguing for it. Tudora, *EPJA* 54, 87 (2018) p. 3
writes the regime as "`E = aT²` where `a` is the level density parameter and `T` the nuclear
temperature" and its eq. (1) carries no shift; *EPJA* 58, 126 (2022) p. 7 does the same. Across
those papers `E1` is never named, used, or discussed — the back-shifted Fermi gas appears only
as the name of the systematics `a` is borrowed from, chosen as a non-energy-dependent stand-in
for the superfluid model on the strength of the shape of `a(A)`.

The omission is not uniform, and does not average away over a cascade. Over the fragment domain
eq. (20) gives `E1` between −1.7 and +1.9 MeV, mean magnitude 0.7 MeV, ordered by parity:
about +0.7 MeV even-even, −0.3 MeV odd-`A`, −1.4 MeV odd-odd. Every emitted neutron flips the
parity of `N`, so the sign alternates down the chain rather than cancelling. Dropping the shift
moves `T` by 1–3 % at the first emission, where `Ē_r` is 10–20 MeV, by 10–12 % near
`Ē_r ≈ 5` MeV, and by 18–40 % at the last sequence. Since `⟨ε⟩ = 2T` for a constant evaporation
cross section, the last figure is a few hundred keV in `⟨ε⟩` — the scale that decides whether
one more neutron is emitted.

The spin-cutoff parameter of the same papers has no consumer here, since `ρ(U, J, π)` is never
evaluated — only the `a` of the energy closure is.

References: T. von Egidy and D. Bucurescu, *Phys. Rev. C* **72**, 044311 (2005); erratum
*Phys. Rev. C* **73**, 049901 (2006); *Phys. Rev. C* **80**, 054310 (2009).

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

`S′ = S + ½P_a′` in MeV, the shell correction plus half the pairing term of the 2009
systematics,

```
½P_a′ = [Δ(A+2, Z+1) − 2Δ(A, Z) + Δ(A−2, Z−1)] / 4
```

The quantity is `P_a′` of the 2009 paper's eq. (12), not the `P_d` of the 2005 paper or of
Audi's tables — those are a different combination, and `P_a′` is the one eq. (19) consumes. It
carries no `(−1)^Z` alternation.

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
julia> table = read_mass_table(joinpath(pkgdir(DeterministicSequentialEmission), "data", "mass_excess", "ame2020.dat"));

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
    offset =
        model.deformed_branch && is_deformed_gilbert_cameron(nuclide) ?
        GILBERT_CAMERON_OFFSET_DEFORMED : GILBERT_CAMERON_OFFSET_UNDEFORMED
    a = nuclide.A * (GILBERT_CAMERON_SHELL_WEIGHT * (S_proton + S_neutron) + offset)
    return a > 0 ? a : nothing
end
