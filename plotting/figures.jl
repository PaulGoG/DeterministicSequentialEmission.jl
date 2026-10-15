# Validation figures: the model against measurement.
#
#   julia plotting/figures.jl [runs directory] [data directory] [figures directory]
#
# Nothing is solved here. The model side is read from what scripts/run.jl wrote under
# data/sims/, so a figure is redrawn from tables rather than by running the model again, and
# the numbers on the page are the numbers in the files. Run the pipeline first:
#
#   julia scripts/run.jl config/Cf252_sf.toml
#
# The measured side is retrieved from EXFOR with ExforFissionData.jl and laid out as
#
#   <data>/<system>/nu_vs_A/*.dat                            A nu nu_uncertainty
#   <data>/<system>/eps_vs_A/*.dat                           A eps eps_uncertainty
#   <data>/<system>/spectrum_maxwellian_ratio_vs_E/*.dat     E spectrum_maxwellian_ratio [...]
#   <data>/<system>/spectrum_cf252_ratio_vs_E/*.dat          E spectrum_cf252_ratio [...]
#   <data>/<system>/Y_vs_A_p/*.dat                           A_p Y Y_uncertainty
#   <data>/<system>/Y_vs_Z_A_p/*.dat                         Z A_p Y Y_uncertainty
#   <data>/<system>/E_K_p_vs_A_p/*.dat                       A_p E_K_p E_K_p_uncertainty
#   <data>/<system>/E_gamma_tot/*.dat                        E_gamma_tot [E_gamma_tot_uncertainty]
#   <data>/Cf252_sf/E_gamma_tot_vs_TKE/*.dat                 TKE E_gamma_tot [...]
#   <data>/Cf252_sf/E_gamma_vs_A/*.dat, E_gamma_tot_vs_A/    A E_gamma [...], A E_gamma_tot [...]
#
# which is the layout that package writes, with a retrieval.toml beside the files that is read
# for the orientation of a ratio, for datasets a later one supersedes, and for the unit and the
# qualifiers of the post-neutron tables. The evaluated product yields under Y_vs_A_p/ and
# Y_vs_Z_A_p/ are staged by scripts/stage_product_yields.jl as <library>.dat with a
# provenance.toml beside them, next to the retrieved sets. The retrieval writes its yields in
# fragments per fission and the evaluations are in per cent; each retrieved table is converted
# from the unit its record states, and a unit the script does not know stops it. The pre-neutron
# mean total kinetic energy a run was solved with is read from the input its configuration.toml
# names, under the same data directory. Columns are taken by position, not by header text. Both
# directories default to data/ under the project and neither is version-controlled.
#
# Each figure is written to results/, or to the figures directory given, as PDF and PNG, beside
# a TOML sidecar that names, for every system drawn, the run directory, the commit recorded in its
# metadata, and the measured, evaluated and input files it read.

include(joinpath(@__DIR__, "activate.jl"))
include(joinpath(@__DIR__, "theme.jl"))

using Printf
using TOML

set_theme!(PUBLICATION_THEME)

const ROOT = dirname(@__DIR__)
const RUNS = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "data", "sims")
const MEASURED = length(ARGS) >= 2 ? ARGS[2] : joinpath(ROOT, "data")
const FIGURES = length(ARGS) >= 3 ? ARGS[3] : joinpath(ROOT, "results")

# LaTeXStrings reaches this environment only through Makie; its string type is taken from the
# macro.
const MathText = typeof(L"")

"""
One fissioning system: the token that names its directories, the typeset notation it is drawn
under, the mass and charge of the compound nucleus, and the evaluated prompt multiplicity it is
compared against: ENDF/B-VIII.0, MF1 MT456 at 0.0253 eV for the induced systems, and for
252-Cf(sf) the constant of the decay sublibrary (MAT 3644), Axton's total 3.7676 less the
delayed 0.0086. The model computes prompt neutrons only, so the total would not compare.

`label` is the system token, as `system_label` returns it; `notation` is the typeset form, as
`system_notation` returns it. They are written out here rather than called, because this
environment draws from tables and does not load the package.
"""
struct Case
    label::String
    notation::MathText
    compound_mass::Int
    compound_charge::Int
    evaluated::Float64
end

const CASES = [
    Case("Cf252_sf", L"$^{252}$Cf(sf)", 252, 98, 3.7590),
    Case("U235_nth", L"$^{235}$U(nth,f)", 236, 92, 2.4140),
    Case("Pu239_nth", L"$^{239}$Pu(nth,f)", 240, 94, 2.8680),
    Case("U233_nth", L"$^{233}$U(nth,f)", 234, 92, 2.4785),
]

# The system every spectrum ratio is taken against.
const REFERENCE = first(CASES)

# Fraction of a run's total `Y(A)` a mass must hold to be drawn: 0.01 % of the fissions per mass
# when the total is 200. The measurement papers stop their ν(A) plots there (Göök 2014,
# doi:10.1103/PhysRevC.90.064611, Fig. 10; Al-Adili 2020, doi:10.1103/PhysRevC.102.064610,
# Fig. 19). On the shipped runs the windows are 81–171 (²⁵²Cf), 76–160 (²³⁵U), 81–159 (²³⁹Pu)
# and 79–155 (²³³U).
const YIELD_WINDOW_FRACTION = 5.0e-5

# Error bars wider than this, in neutrons, do not set the ordinate limits of the ν(A) grid.
const MULTIPLICITY_BAR_LIMIT = 0.5

# Error bars wider than this, in MeV, do not set the ordinate limits of the ⟨ε⟩(A) grid.
const NEUTRON_ENERGY_BAR_LIMIT = 0.3

# Energies drawn in the spectra. Below the window the measured ratio is dominated by detector
# thresholds, above it by counting statistics in the tail.
const ENERGY_WINDOW = (0.08, 11.0)

# Measured mean laboratory energy of the 252-Cf prompt fission neutron spectrum [MeV].
const MEASURED_MEAN_ENERGY = 2.13

# Maxwellian temperature taken when a run's configuration records none [MeV].
const FALLBACK_MAXWELLIAN_TEMPERATURE = 1.42

# Margin on either side of the data range, as a fraction of that range; the band reserved
# beyond the data for one line of in-axis annotation, as a fraction of the data range; and the
# fraction of the axis width, from the annotated side, that those lines cover.
const MARGIN = 0.05
const ANNOTATION_BAND = 0.11
const ANNOTATION_EXTENT = 0.45

# The per-system comparison figures share one grid: two columns, each panel with its legend
# above it, 1200 wide and 500 per row. The legends of one figure all have the same number of
# rows, so the axes of one grid row are the same height.
const GRID_COLUMNS = 2
const GRID_WIDTH = 1200
const GRID_ROW_HEIGHT = 500
const LEGEND_COLUMNS = 2

# Height of one legend row at the base font size, by which a grid row grows where its legends
# need more than two rows.
const LEGEND_ROW_HEIGHT = 34

const SPECTRUM_TICKS =
    ([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0], ["0.1", "0.2", "0.5", "1", "2", "5", "10"])

# Sources of the isobaric charge distribution. Colour and line style are fixed per source.
const CHARGE_SOURCES = (
    (
        file = "charge_distribution_vs_A.dat",
        label = "Used by the run",
        linestyle = :solid,
        colour = BLACK,
    ),
    (
        file = "charge_distribution_vs_A_wahl_per_reaction.dat",
        label = L"$Z_p$ model, fitted parameters",
        linestyle = :dashdot,
        colour = REDDISH_PURPLE,
    ),
    (
        file = "charge_distribution_vs_A_wahl_systematics.dat",
        label = L"$Z_p$ model, systematics",
        linestyle = :dash,
        colour = VERMILION,
    ),
    (
        file = "charge_distribution_vs_A_mean.dat",
        label = "Mean values",
        linestyle = :dot,
        colour = SKY_BLUE,
    ),
)

# Evaluated fission product yield libraries the product figures draw, each `<file>` under
# `<data>/<system>/Y_vs_A_p/` and `<data>/<system>/Y_vs_Z_A_p/`, staged by
# scripts/stage_product_yields.jl. Evaluations are lines, not markers: they are not
# measurements. The isotopic distributions take them as markers, the one given here, so that
# the model's line stays legible between the points of the two libraries.
const EVALUATED_YIELDS = (
    (
        file = "ENDF-B-VIII.0.dat",
        label = "ENDF/B-VIII.0",
        linestyle = :dash,
        colour = VERMILION,
        marker = :circle,
    ),
    (
        file = "JEFF-3.3.dat",
        label = "JEFF-3.3",
        linestyle = :dashdot,
        colour = BLUE,
        marker = :rect,
    ),
)

# Units the retrieval writes its post-neutron tables in, by observable, with the factor that
# takes each to the unit the figures draw: per cent per fission for a yield, 200 over both
# products as the evaluations and the model's tables are, and MeV for a kinetic energy. The
# retrieval writes yields in fragments per fission, 2 over both products.
const PRODUCT_UNITS = Dict(
    "Y_vs_A_p" => Dict("PART/FIS" => 100.0, "PC/FIS" => 1.0),
    "Y_vs_Z_A_p" => Dict("PART/FIS" => 100.0, "PC/FIS" => 1.0),
    "E_K_p_vs_A_p" => Dict("MEV" => 1.0),
)

# Retrieved post-neutron datasets left out of the figures, by accession, with the reason. A set
# its record flags as superseded is left out as well, and the one superseding it is drawn.
const EXCLUDED = Dict(
    "31738002" => "its heavy-product yields alone sum to 182 % per fission, where one side \
        sums to 100 %, and it states no uncertainty",
)

# Archive qualifiers that mark a drawn set, with open markers: SPA, averaged over a neutron
# spectrum other than the one it is filed under, and DERIV, derived from other data rather
# than measured.
const MARKED_QUALIFIERS = ("SPA", "DERIV")

# The pooled measured families of the product figures: the chain yields, muted beneath the
# named sets of the mass yield, and the independent yields of the isotopic distributions.
const CHAIN_COLOUR = RGBf(0.62, 0.62, 0.62)
const INDEPENDENT_COLOUR = BLUISH_GREEN
const INDEPENDENT_MARKER = :diamond
const INDEPENDENT_LABEL = "Measured, EXFOR independent yields"

# Measured kinetic-energy series kept by name in a panel that holds more than the product palette
# can tell apart; the others are pooled into one.
const NAMED_KINETIC_ENERGY_SETS = 4

# Marker sizes of the measured series of the product mass yield and kinetic energy. Up to some
# 700 points share a panel there, a mass unit some 5 px wide, and markers of the base size 14
# merge into bands; the chain yields are smaller still, beneath the named sets.
const DENSE_MARKER_SIZE = 10
const CHAIN_MARKER_SIZE = 7

# Lower limit of the logarithmic yield axes [%], below which the yield figures draw no product;
# and the product yield [%] below which a mean kinetic energy, taken over a handful of
# configurations, is not drawn.
const YIELD_FLOOR = 1.0e-3
const MEAN_FLOOR = 1.0e-2

# Neutron shell closures marked in the neutron-number distribution.
const NEUTRON_SHELLS = (50, 82)

# Elements drawn per side in the isotopic distributions.
const ELEMENTS_PER_SIDE = 4

# Chemical symbols by atomic number, as far as the compound nuclei drawn.
const ELEMENT_SYMBOLS = split(
    """
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As
    Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu
    Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np
    Pu Am Cm Bk Cf
    """,
)

"""
The directory of the run to draw for one system, or `nothing` when none has been written.

A run writes to `data/sims/<system>/<run identifier>`, so a system may hold several. The
most recently written one is drawn, runs moved aside with a `#N` suffix are skipped, and
which that was is stated on the terminal rather than left to be inferred from the figure.
"""
function run_directory(case::Case)
    system = joinpath(RUNS, case.label)
    isdir(system) || return nothing
    runs = filter(
        path -> isdir(path) && !occursin('#', basename(path)),
        readdir(system; join = true),
    )
    isempty(runs) && return nothing
    return last(sort(runs; by = mtime))
end

"The totals a run wrote, by name."
function read_totals(directory::AbstractString)
    path = joinpath(directory, "totals.dat")
    totals = Dict{String, Float64}()
    isfile(path) || return totals
    for (index, line) in enumerate(eachline(path))
        index == 1 && continue
        fields = split(strip(line))
        length(fields) >= 2 || continue
        totals[fields[1]] = parse(Float64, fields[2])
    end
    return totals
end

"""
    yield_window(directory) -> Tuple{Float64, Float64}

The fragment masses a run's figures are drawn over: the least and the greatest mass whose
`Y(A)`, read from the run's `Y_vs_A.dat`, holds at least `YIELD_WINDOW_FRACTION` of the table's
total. Throws an `ArgumentError` naming the directory when no mass qualifies.
"""
function yield_window(directory::AbstractString)
    A, Y, _ = read_columns(joinpath(directory, "Y_vs_A.dat"))
    keep = Y .>= YIELD_WINDOW_FRACTION * sum(Y)
    any(keep) || throw(
        ArgumentError("no mass of $directory holds $(YIELD_WINDOW_FRACTION) of the yield"),
    )
    return (minimum(A[keep]), maximum(A[keep]))
end

"The commit recorded in a run's `metadata.toml`, under `[run]` or at top level."
function gitcommit(directory::AbstractString)
    path = joinpath(directory, "metadata.toml")
    metadata = isfile(path) ? TOML.parsefile(path) : Dict{String, Any}()
    run = get(metadata, "run", Dict{String, Any}())
    commit = run isa AbstractDict ? get(run, "gitcommit", nothing) : nothing
    commit === nothing && (commit = get(metadata, "gitcommit", nothing))
    if commit === nothing
        @warn "no gitcommit recorded for the run" path
        return "unrecorded"
    end
    return string(commit)
end

"The Maxwellian temperature [MeV] a run's configuration takes the spectrum ratio against."
function maxwellian_temperature(directory::AbstractString)
    path = joinpath(directory, "configuration.toml")
    if isfile(path)
        spectrum = get(TOML.parsefile(path), "spectrum", nothing)
        if spectrum isa AbstractDict
            temperature = get(spectrum, "maxwellian_temperature", nothing)
            temperature isa Real && return Float64(temperature)
        end
    end
    @warn "no spectrum.maxwellian_temperature in the run configuration; taking \
        $(FALLBACK_MAXWELLIAN_TEMPERATURE) MeV" path
    return FALLBACK_MAXWELLIAN_TEMPERATURE
end

"""
A measured dataset: its author–year label, the files it was read from — several where one
measurement was tabulated in parts — and the points inside the drawn window with their
uncertainties, zero where none is reported.
"""
struct MeasuredSet
    label::String
    paths::Vector{String}
    x::Vector{Float64}
    y::Vector{Float64}
    σ::Vector{Float64}
end

# Files are named <accession>_<author>_<year>.dat by the retrieval.
accession(name::AbstractString) = String(first(split(basename(name), "_")))

function author_year(name::AbstractString)
    parts = split(replace(basename(name), ".dat" => ""), "_")
    length(parts) >= 3 || return String(name)
    author = replace(parts[2], "." => ". ")
    # Transliterated names restored; EXFOR writes one author as both "Liu Chao" and
    # "Chao Liu", and Liu is the family name.
    author = replace(
        author,
        "Goeoek" => "Göök",
        "Budtz-jorgensen" => "Budtz-Jørgensen",
        "Al-adili" => "Al-Adili",
        "Guenther" => "Günther",
        "LiuChao" => "Liu",
        "ChaoLiu" => "Liu",
    )
    return strip(last(split(author, ". "))) * " " * parts[3]
end

# The measured datasets each figure shows, by EXFOR accession: a fixed selection per system, for
# legibility. The data directory holds the whole retrieval. A system or observable not listed
# shows every dataset retrieved. The post-neutron observables are curated by `load_products`
# instead.
const SHOWN = Dict(
    ("Cf252_sf", "nu_vs_A") =>
        ["23118006", "23175008", "23268005", "23764010", "41425014"],
    ("U235_nth", "nu_vs_A") =>
        ["21095002", "22464004", "23764004", "30909004", "41516012"],
    ("Pu239_nth", "nu_vs_A") => ["14369002", "23012008", "41397002"],
    ("U233_nth", "nu_vs_A") => ["14369004", "22660005", "23595005", "41397004"],
    ("Cf252_sf", "spectrum_maxwellian_ratio_vs_E") =>
        ["14278002", "30775002", "40418008", "41158003"],
)

"""
The accepted datasets of one observable of one system, keyed by accession, from the retrieval
record beside the files; empty where there is none.
"""
function retrieval_record(system::AbstractString, observable::AbstractString)
    path = joinpath(MEASURED, system, observable, "retrieval.toml")
    record = Dict{String, Any}()
    isfile(path) || return record
    for entry in get(TOML.parsefile(path), "accepted", Any[])
        entry isa AbstractDict && haskey(entry, "file") || continue
        record[accession(entry["file"])] = entry
    end
    return record
end

"""
Files of one EXFOR entry — the same first five digits of the accession — under one author–year
label are parts of one measurement, tabulated in subentries by energy range or detector, and
are joined under that label in order of abscissa. Files of different entries under one label
are distinct measurements, and each label carries its accession.
"""
function join_complementary(sets::Vector{MeasuredSet})
    entry(set::MeasuredSet) = first(accession(first(set.paths)), 5)
    joined = MeasuredSet[]
    for label in unique(set.label for set in sets)
        parts = filter(set -> set.label == label, sets)
        entries = unique(entry.(parts))
        for number in entries
            members = filter(part -> entry(part) == number, parts)
            members = members[sortperm(members; by = part -> first(part.x))]
            name = length(entries) == 1 ? label : "$label ($number)"
            push!(
                joined,
                MeasuredSet(
                    name,
                    reduce(vcat, (part.paths for part in members)),
                    reduce(vcat, (part.x for part in members)),
                    reduce(vcat, (part.y for part in members)),
                    reduce(vcat, (part.σ for part in members)),
                ),
            )
        end
    end
    return joined
end

"""
The measured datasets of one observable of one system, restricted to `window` and to `SHOWN`
where the pair is listed, less the accessions in `skip`; those in `invert` are taken as their
reciprocal, the uncertainty propagated as `σ/y²`. Parts of one measurement are joined by
[`join_complementary`](@ref).
"""
function load_measured(
    system::AbstractString,
    observable::AbstractString,
    window;
    skip::Vector{String} = String[],
    invert::Vector{String} = String[],
)
    directory = joinpath(MEASURED, system, observable)
    sets = MeasuredSet[]
    isdir(directory) || return sets
    shown = get(SHOWN, (system, observable), nothing)
    for name in sort(filter(endswith(".dat"), readdir(directory)))
        tag = accession(name)
        shown === nothing || tag in shown || continue
        tag in skip && continue
        path = joinpath(directory, name)
        x, y, σ = read_columns(path)
        if tag in invert
            σ = σ ./ y .^ 2
            y = 1 ./ y
        end
        keep = first(window) .<= x .<= last(window)
        any(keep) || continue
        push!(sets, MeasuredSet(author_year(name), [path], x[keep], y[keep], σ[keep]))
    end
    return join_complementary(sets)
end

"""
A retrieved post-neutron dataset in the unit the figures draw: its accession, author–year label,
the quantity modifier of its reaction code (`SEC`, `CHN` or `IND`), whether a qualifier in
`MARKED_QUALIFIERS` marks it, its file, and its rows — the charge where the table is resolved in
it, the product mass, the value and its uncertainty, zero where none is stated.
"""
struct ProductSet
    accession::String
    label::String
    kind::String
    marked::Bool
    path::String
    Z::Vector{Int}
    A_p::Vector{Float64}
    y::Vector{Float64}
    σ::Vector{Float64}
end

"""
    unit_factor(observable, unit) -> Float64

The factor that takes a retrieved table of `observable`, written in `unit`, to the unit the
figures draw it in. An observable or a unit that `PRODUCT_UNITS` does not list throws an
`ArgumentError`: a table is converted from the unit its record states, never by assumption.

```julia
unit_factor("Y_vs_A_p", "PART/FIS")   # 100.0: fragments per fission to per cent
```
"""
function unit_factor(observable::AbstractString, unit::AbstractString)
    units = get(PRODUCT_UNITS, observable, nothing)
    units === nothing &&
        throw(ArgumentError("no unit conversion is defined for $observable"))
    factor = get(units, unit, nothing)
    known = join(sort!(collect(keys(units))), ", ")
    factor === nothing &&
        throw(ArgumentError("unknown unit \"$unit\" for $observable; known: $known"))
    return factor
end

"The quantity modifier of an EXFOR reaction code: `CHN` in `94-PU-239(N,F)MASS,CHN,FY,,SPA`."
function quantity_modifier(reaction_code::AbstractString)
    fields = split(last(split(reaction_code, ')')), ',')
    length(fields) >= 2 ||
        throw(ArgumentError("no quantity modifier in the reaction code $reaction_code"))
    return String(fields[2])
end

"""
    load_products(system, observable) -> (sets, excluded)

The retrieved datasets of one post-neutron observable of one system — `Y_vs_A_p`, `Y_vs_Z_A_p`
or `E_K_p_vs_A_p` — each converted by [`unit_factor`](@ref) from the unit its record states as
`unit_written`, beside the accessions left out, each with the reason: a set in `EXCLUDED`, and a
set the record flags as superseded. The evaluations staged in the same directory are not read
here. A table the record does not name, or names without a unit, throws an `ArgumentError`.
"""
function load_products(system::AbstractString, observable::AbstractString)
    directory = joinpath(MEASURED, system, observable)
    sets = ProductSet[]
    excluded = Pair{String, String}[]
    isdir(directory) || return sets, excluded
    path = joinpath(directory, "retrieval.toml")
    record = isfile(path) ? TOML.parsefile(path) : Dict{String, Any}()
    entries = Dict{String, Any}(
        entry["file"] => entry for entry in get(record, "accepted", Any[]) if
        entry isa AbstractDict && haskey(entry, "file")
    )
    evaluations = Set(library.file for library in EVALUATED_YIELDS)
    offset = observable == "Y_vs_Z_A_p" ? 1 : 0
    for name in sort(filter(endswith(".dat"), readdir(directory)))
        name in evaluations && continue
        file = joinpath(directory, name)
        entry = get(entries, name, nothing)
        entry === nothing && throw(
            ArgumentError(
                "$file is not named in the retrieval record beside it, so its unit is unknown",
            ),
        )
        tag = accession(name)
        qualifiers = get(entry, "qualifiers", String[])
        if haskey(EXCLUDED, tag)
            push!(excluded, tag => EXCLUDED[tag])
            continue
        end
        index = findfirst(startswith("superseded"), qualifiers)
        if index !== nothing
            by = match(r"by (\w+)", qualifiers[index])
            push!(excluded, tag => "superseded by $(by === nothing ? "?" : by[1])")
            continue
        end
        unit = get(entry, "unit_written", nothing)
        unit isa AbstractString ||
            throw(ArgumentError("no unit_written for $name in $path"))
        factor = unit_factor(observable, unit)
        Z = Int[]
        A_p = Float64[]
        y = Float64[]
        σ = Float64[]
        for row in read_rows(file)
            length(row) >= 2 + offset && isfinite(row[2 + offset]) || continue
            offset == 1 && push!(Z, round(Int, row[1]))
            push!(A_p, row[1 + offset])
            push!(y, factor * row[2 + offset])
            u = length(row) >= 3 + offset ? row[3 + offset] : NaN
            push!(σ, isfinite(u) ? factor * u : 0.0)
        end
        codes = [String(strip(first(split(q, ':')))) for q in qualifiers]
        push!(
            sets,
            ProductSet(
                tag,
                author_year(name),
                quantity_modifier(entry["reaction_code"]),
                any(in(MARKED_QUALIFIERS), codes),
                file,
                Z,
                A_p,
                y,
                σ,
            ),
        )
    end
    for (tag, reason) in excluded
        println(
            "  ",
            rpad(system, 10),
            " ",
            rpad(observable, 13),
            " not drawn: ",
            tag,
            ", ",
            reason,
        )
    end
    return sets, excluded
end

"""
One marker series of a product figure: its legend label, the files drawn, the points, and which
of them belong to a set that a qualifier in `MARKED_QUALIFIERS` marks.
"""
struct ProductSeries
    label::String
    paths::Vector{String}
    x::Vector{Float64}
    y::Vector{Float64}
    σ::Vector{Float64}
    marked::Vector{Bool}
end

"""
    pooled_series(label, sets, window) -> ProductSeries

The points of `sets` whose product mass lies in `window`, pooled under one `label`; a set with no
point in the window does not enter, and the paths name the sets that do.
"""
function pooled_series(label::AbstractString, sets, window)
    paths = String[]
    x = Float64[]
    y = Float64[]
    σ = Float64[]
    marked = Bool[]
    for set in sets
        inside = first(window) .<= set.A_p .<= last(window)
        any(inside) || continue
        push!(paths, set.path)
        append!(x, set.A_p[inside])
        append!(y, set.y[inside])
        append!(σ, set.σ[inside])
        append!(marked, fill(set.marked, count(inside)))
    end
    return ProductSeries(label, paths, x, y, σ, marked)
end

"""
    named_series(sets, window) -> Vector{ProductSeries}

One series per author–year label of `sets`, in `window`. Tables of one EXFOR entry under one
label are one series, as [`join_complementary`](@ref) joins them; tables of different entries
under one label each carry their entry.
"""
function named_series(sets, window)
    entry(set) = first(set.accession, 5)
    series = ProductSeries[]
    for label in unique(set.label for set in sets)
        members = filter(set -> set.label == label, sets)
        entries = unique(entry.(members))
        for number in entries
            name = length(entries) == 1 ? label : "$label ($number)"
            s = pooled_series(name, filter(m -> entry(m) == number, members), window)
            isempty(s.x) || push!(series, s)
        end
    end
    return series
end

"""
    draw_series!(axis, series, colour, marker; markersize = 14, floor = -Inf)

Error bars, cut off at `floor` for a logarithmic axis, then markers, of one product series:
filled for the points of unmarked sets, open for those of marked ones. Returns the plot that
stands for the series in a legend: the filled markers unless every point is marked.
"""
function draw_series!(
    axis::Axis,
    series::ProductSeries,
    colour,
    marker;
    markersize = 14,
    floor::Real = -Inf,
)
    bars = series.σ .> 0
    any(bars) && rangebars!(
        axis,
        series.x[bars],
        max.(series.y[bars] .- series.σ[bars], floor),
        series.y[bars] .+ series.σ[bars];
        color = (colour, 0.5),
        linewidth = 1.5,
    )
    style = (; marker, markersize, strokecolor = darker(colour))
    handle = nothing
    filled = .!series.marked
    if any(filled)
        handle = scatter!(
            axis,
            series.x[filled],
            series.y[filled];
            color = colour,
            strokewidth = 1.5,
            style...,
        )
    end
    if any(series.marked)
        open = scatter!(
            axis,
            series.x[series.marked],
            series.y[series.marked];
            color = :transparent,
            strokewidth = 2,
            style...,
        )
        handle === nothing && (handle = open)
    end
    return handle
end

# The legend entry that explains the open markers of the product figures: an open marker
# beside the qualifiers it stands for, grey where several families are drawn.
const MARKED_LABEL = "SPA or DERIV"
marked_key(; markersize = 14, marker = :circle, colour = GUIDE_LABEL_COLOUR) =
    MarkerElement(;
        marker = marker,
        color = :transparent,
        strokecolor = darker(colour),
        strokewidth = 2,
        markersize = markersize,
    )

function available(case::Case)
    directory = run_directory(case)
    if directory !== nothing && isfile(joinpath(directory, "nu_vs_A.dat"))
        println("  ", rpad(case.label, 10), " ", basename(directory))
        return true
    end
    @warn """
        no run to draw; solve it first with `julia scripts/run.jl \
        config/$(case.label).toml`.""" expected = joinpath(RUNS, case.label)
    return false
end

"Axis limits spanning the finite `values`, widened by `MARGIN` of their range on either side."
function padded_limits(values)
    finite = filter(isfinite, collect(Float64, values))
    isempty(finite) && throw(ArgumentError("no finite values to set axis limits from"))
    low, high = extrema(finite)
    span = high > low ? high - low : max(abs(high), 1.0)
    return (low - MARGIN * span, high + MARGIN * span)
end

"As `padded_limits`, with the margin taken in decades for a logarithmic axis."
function padded_log_limits(values)
    positive = filter(v -> isfinite(v) && v > 0, collect(Float64, values))
    isempty(positive) &&
        throw(ArgumentError("no positive values to set logarithmic axis limits from"))
    low, high = padded_limits(log10.(positive))
    return (exp10(low), exp10(high))
end

"""
Ordinate limits spanning a model curve and every measured value, and also the error bars, on
the side `bars` (`:upper` or `:lower`), of the points whose abscissa lies in `cleared`, the
horizontal extent of an annotation set beyond the data on that side. Elsewhere error bars may
reach into the margins and the annotation band. Bars wider than `bar_limit` are ignored when the
limits are set, so that a few poorly measured points at the edges do not stretch the axis; the
points themselves still count.
"""
function ordinate_limits(
    model::AbstractVector{<:Real},
    sets::Vector{MeasuredSet};
    cleared::Union{Nothing, Tuple{Real, Real}} = nothing,
    bars::Symbol = :upper,
    bar_limit::Real = Inf,
)
    values = collect(Float64, model)
    sign =
        bars === :upper ? 1 :
        bars === :lower ? -1 :
        throw(ArgumentError("bars must be :upper or :lower, got :$bars"))
    for set in sets
        append!(values, set.y)
        cleared === nothing && continue
        under = (first(cleared) .<= set.x .<= last(cleared)) .& (set.σ .<= bar_limit)
        append!(values, set.y[under] .+ sign .* set.σ[under])
    end
    return padded_limits(values)
end

panel(slot; kwargs...) = Axis(slot; xlabelpadding = 8, ylabelpadding = 8, kwargs...)

"Error bars, then markers, of one measured dataset; returns the marker plot for the legend."
function draw_measured!(axis::Axis, set::MeasuredSet)
    colour, marker = dataset_style(set.label)
    any(>(0), set.σ) &&
        errorbars!(axis, set.x, set.y, set.σ; color = (colour, 0.5), linewidth = 1.5)
    return scatter!(
        axis,
        set.x,
        set.y;
        color = colour,
        marker = marker,
        markersize = 14,
        strokecolor = darker(colour),
        strokewidth = 1.5,
    )
end

"""
    label_guide!(slot, value, label, limits; side)

Label a reference line where it meets the axis frame: above the frame for a vertical guide at
`x = value` (`side = :top`), right of it for a horizontal guide at `y = value` (`side = :right`).
A transparent overlay axis in the same layout cell carries the label as its only tick, so the
layout reserves room for it and it cannot collide with the data.
"""
function label_guide!(slot, value::Real, label, limits::Tuple{Real, Real}; side::Symbol)
    common = (; backgroundcolor = :transparent, xgridvisible = false, ygridvisible = false)
    if side === :top
        overlay = Axis(
            slot;
            common...,
            xaxisposition = :top,
            xticks = ([value], [label]),
            xticksvisible = false,
            xticklabelsize = ANNOTATION_SIZE,
            xticklabelcolor = GUIDE_LABEL_COLOUR,
            xticklabelpad = 6,
            xlabelvisible = false,
        )
        hideydecorations!(overlay)
        xlims!(overlay, limits...)
    elseif side === :right
        overlay = Axis(
            slot;
            common...,
            yaxisposition = :right,
            yticks = ([value], [label]),
            yticksvisible = false,
            yticklabelsize = ANNOTATION_SIZE,
            yticklabelcolor = GUIDE_LABEL_COLOUR,
            yticklabelpad = 8,
            ylabelvisible = false,
        )
        hidexdecorations!(overlay)
        ylims!(overlay, limits...)
    else
        throw(ArgumentError("side must be :top or :right, got :$side"))
    end
    hidespines!(overlay)
    return overlay
end

"""
The provenance entry of one system: the run drawn, its commit, the measured files used — of the
`MeasuredSet`s or `ProductSeries` drawn — the evaluated ones where a figure draws evaluations, the
inputs of the run a guide is read from, and the retrieved sets left out, each with the reason.
"""
function system_record(
    directory::AbstractString,
    sets::AbstractVector;
    evaluated::Vector{String} = String[],
    inputs::Vector{String} = String[],
    excluded::Vector{Pair{String, String}} = Pair{String, String}[],
)
    measured = String[]
    for set in sets, path in set.paths
        push!(measured, relpath(path, MEASURED))
    end
    record = Dict{String, Any}(
        "run" => basename(directory),
        "gitcommit" => gitcommit(directory),
        "measured" => measured,
    )
    isempty(evaluated) ||
        (record["evaluated"] = [relpath(path, MEASURED) for path in evaluated])
    isempty(inputs) || (record["input"] = [relpath(path, MEASURED) for path in inputs])
    isempty(excluded) ||
        (record["excluded"] = ["$tag: $reason" for (tag, reason) in excluded])
    return record
end

function export_figure(figure::Figure, name::AbstractString, systems::Dict{String, Any})
    path = joinpath(FIGURES, name)
    save(path * ".pdf", figure)
    save(path * ".png", figure; px_per_unit = 4)
    open(path * ".toml", "w") do io
        TOML.print(
            io,
            Dict{String, Any}("figure" => name, "systems" => systems);
            sorted = true,
        )
    end
    println("wrote ", path, ".pdf")
    return nothing
end

"""
Stacked lines of in-axis annotation in one corner, `ANNOTATION_LINE` apart, read top to bottom:
`:upper_left`, `:upper_right` or `:lower_left`.
"""
function annotate!(axis::Axis, lines::AbstractVector; corner::Symbol)
    # Lines read top to bottom in either corner: the first line is the outermost one in the
    # upper corner and the innermost one in the lower corner.
    if corner === :upper_left
        position, align = (0.03, 0.96), (:left, :top)
        step = index -> -(index - 1) * ANNOTATION_LINE
    elseif corner === :upper_right
        position, align = (0.97, 0.96), (:right, :top)
        step = index -> -(index - 1) * ANNOTATION_LINE
    elseif corner === :lower_left
        position, align = (0.03, 0.04), (:left, :bottom)
        step = index -> (length(lines) - index) * ANNOTATION_LINE
    else
        throw(
            ArgumentError(
                "corner must be :upper_left, :upper_right or :lower_left, got :$corner",
            ),
        )
    end
    for (index, line) in pairs(lines)
        text!(
            axis,
            position...;
            text = line,
            space = :relative,
            align = align,
            offset = (0, step(index)),
            fontsize = ANNOTATION_SIZE,
        )
    end
    return nothing
end

"A legend above its panel: `rows` banks of `LEGEND_COLUMNS` entries, frameless, flush."
function panel_legend!(cell::GridLayout, handles, labels, rows::Integer)
    Legend(
        cell[1, 1],
        handles,
        labels;
        nbanks = rows,
        tellwidth = false,
        tellheight = true,
        padding = (0, 0, 0, 0),
        colgap = 16,
        patchsize = (35, 20),
    )
    rowgap!(cell, 8)
    return nothing
end

"The rows every legend of a grid takes: enough for the panel with the most entries."
legend_rows(counts) = maximum(cld(count, LEGEND_COLUMNS) for count in counts; init = 1)

"""
    grid_slots(figure, index; legend = true) -> (legend position, axis position)

The layout positions of the legend and the axis of panel `index` of a per-system grid. With
legends, each grid row is two rows of the figure, the legends above the axes, so that the legends
of one grid row share one height and its axes align, whatever the number of entries of each
legend; the legend position is `nothing` without them.
"""
function grid_slots(figure::Figure, index::Integer; legend::Bool = true)
    row, column = fldmod1(index, GRID_COLUMNS)
    legend || return nothing, figure[row, column]
    return figure[2row - 1, column], figure[2row, column]
end

"""
    grid_legend!(position, handles, labels, rows)

The legend of one panel of a grid laid out by [`grid_slots`](@ref): `rows` banks of
`LEGEND_COLUMNS` entries, frameless, at the foot of its cell, flush with the axis below it.
"""
function grid_legend!(position, handles, labels, rows::Integer)
    Legend(
        position,
        handles,
        labels;
        nbanks = rows,
        tellwidth = false,
        tellheight = true,
        valign = :bottom,
        padding = (0, 0, 0, 0),
        colgap = 16,
        patchsize = (35, 20),
    )
    return nothing
end

"""
    grid_gaps!(figure, rows; legend = true)

The gaps of a grid of `rows` rows laid out by [`grid_slots`](@ref): 8 between a legend and its
axis, 12 between grid rows and between the columns.
"""
function grid_gaps!(figure::Figure, rows::Integer; legend::Bool = true)
    if legend
        for gap in 1:(2rows - 1)
            rowgap!(figure.layout, gap, isodd(gap) ? 8 : 12)
        end
    else
        rows > 1 && rowgap!(figure.layout, 12)
    end
    colgap!(figure.layout, 12)
    return nothing
end

"One panel of a per-system grid: the system, its run, the model curve and the measured sets."
struct SystemPanel
    case::Case
    directory::String
    x::Vector{Float64}
    y::Vector{Float64}
    sets::Vector{MeasuredSet}
    window::Tuple{Float64, Float64}
end

"""
The panels of one quantity against fragment mass: for every system with a run that wrote
`table`, the model curve and the measured sets of `observable` inside the run's own mass range
— measured points are compared only where the run solved the model — within its yield window.
"""
function mass_panels(
    drawable::Vector{Case},
    table::AbstractString,
    observable::AbstractString,
)
    panels = SystemPanel[]
    for case in drawable
        directory = run_directory(case)
        path = joinpath(directory, table)
        if !isfile(path)
            @warn "no $table to draw" directory
            continue
        end
        A, value, _ = read_columns(path)
        yielded = yield_window(directory)
        window = (max(first(yielded), minimum(A)), min(last(yielded), maximum(A)))
        sets = load_measured(case.label, observable, window)
        push!(panels, SystemPanel(case, directory, A, value, sets, window))
    end
    return panels
end

"""
Draw a quantity against fragment mass for every panel, on the shared grid: the measured sets as
markers, the model as the black line, the A₀/2 guide labelled above the frame, each panel's
legend above it, and `annotate(panel)` as lines in the upper left, above the data. Each panel
takes its own mass range, the run's yield window (`yield_window`), labelled on every panel since
the ranges differ, and its own ordinate range;
error bars wider than `bar_limit` do not set the ordinate range.
"""
function figure_mass_grid(
    name::AbstractString,
    panels::Vector{SystemPanel},
    ylabel,
    annotate::Function;
    bar_limit::Real = Inf,
)
    isempty(panels) && return nothing

    rows = cld(length(panels), GRID_COLUMNS)
    banks = legend_rows(length(p.sets) + 1 for p in panels)
    figure = Figure(size = (GRID_WIDTH, GRID_ROW_HEIGHT * rows))
    provenance = Dict{String, Any}()
    for (index, p) in pairs(panels)
        row, column = fldmod1(index, GRID_COLUMNS)
        cell = GridLayout(figure[row, column])
        has_legend = !isempty(p.sets)
        slot = cell[has_legend ? 2 : 1, 1]
        symmetric = p.case.compound_mass / 2
        xlimits = padded_limits([first(p.window), last(p.window), symmetric])
        annotated = (
            first(xlimits),
            first(xlimits) + ANNOTATION_EXTENT * (last(xlimits) - first(xlimits)),
        )
        axis = panel(
            slot;
            ylabel = ylabel,
            xlabel = L"Fragment mass $A$",
            ylabelvisible = column == 1,
            xticks = mass_ticks(xlimits),
        )

        vlines!(axis, [symmetric]; color = GUIDE_COLOUR, linestyle = :dash, linewidth = 1.5)
        begin_panel!()
        handles = Any[]
        labels = Any[]
        for set in p.sets
            push!(handles, draw_measured!(axis, set))
            push!(labels, set.label)
        end
        push!(handles, lines!(axis, p.x, p.y; color = :black, linewidth = 3))
        push!(labels, "This model")

        lines = annotate(p)
        low, high = ordinate_limits(p.y, p.sets; cleared = annotated, bar_limit = bar_limit)
        ylimits = (low, high + length(lines) * ANNOTATION_BAND * (high - low))
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
        label_guide!(slot, symmetric, L"A_0/2", xlimits; side = :top)
        # The upper left is clear in every panel: light fragments far from symmetry emit least
        # and their neutrons are the slowest.
        annotate!(axis, lines; corner = :upper_left)

        has_legend && panel_legend!(cell, handles, labels, banks)
        provenance[p.case.label] = system_record(p.directory, p.sets)
    end
    rowgap!(figure.layout, 12)
    colgap!(figure.layout, 12)

    export_figure(figure, name, provenance)
    return nothing
end

"""
Draw ν(A) for every system with a run against the measured datasets, each panel annotated with
the system, ⟨ν⟩ and the evaluated value.
"""
function figure_multiplicity(drawable::Vector{Case})
    panels = mass_panels(drawable, "nu_vs_A.dat", "nu_vs_A")
    function annotate(p::SystemPanel)
        modelled = get(read_totals(p.directory), "nu_per_fission", NaN)
        evaluated = p.case.evaluated
        deviation = 100 * (modelled - evaluated) / evaluated
        @printf(
            "  %-10s ⟨ν⟩ = %.4f   evaluated %.4f   %+.1f %%\n",
            p.case.label,
            modelled,
            evaluated,
            deviation
        )
        # Three short lines, so that none reaches the A₀/2 guide, which lies at 0.42 of the
        # axis width for ²³³U.
        return MathText[
            p.case.notation,
            MathText(
                @sprintf(
                    "\$\\langle \\nu \\rangle = %.2f\$ (\$%+.1f\$ %%)",
                    modelled,
                    deviation
                )
            ),
            MathText(@sprintf("Evaluated %.2f", evaluated)),
        ]
    end
    figure_mass_grid(
        "multiplicity_vs_experiment",
        panels,
        L"\nu(A)",
        annotate;
        bar_limit = MULTIPLICITY_BAR_LIMIT,
    )
    return nothing
end

"""
Draw the mean centre-of-mass energy of the prompt neutrons against fragment mass, ⟨ε⟩(A), for
every system with a run against the measured datasets, which are centre-of-mass means as the
retrieval records them.
"""
function figure_mean_neutron_energy(drawable::Vector{Case})
    panels = mass_panels(drawable, "eps_vs_A.dat", "eps_vs_A")
    figure_mass_grid(
        "mean_neutron_energy_vs_experiment",
        panels,
        L"\langle \varepsilon \rangle(A)\ [\mathrm{MeV}]",
        p -> MathText[p.case.notation];
        bar_limit = NEUTRON_ENERGY_BAR_LIMIT,
    )
    return nothing
end

maxwellian(E, T) = 2 / sqrt(π) * T^(-3 / 2) * sqrt(E) * exp(-E / T)

"The laboratory spectrum a run wrote, `(E, N)`, or `nothing`."
function read_spectrum(directory::AbstractString)
    path = joinpath(directory, "spectrum_vs_E.dat")
    isfile(path) || return nothing
    E, N, _ = read_columns(path)
    return (E, N)
end

"Linear interpolation in log–log of a positive tabulated function; `NaN` outside the table."
function loglog_interpolate(xs::Vector{Float64}, ys::Vector{Float64}, x::Real)
    (x < first(xs) || x > last(xs)) && return NaN
    index = searchsortedlast(xs, x)
    index == length(xs) && return ys[end]
    x₁, x₂ = log(xs[index]), log(xs[index + 1])
    y₁, y₂ = log(ys[index]), log(ys[index + 1])
    return exp(y₁ + (y₂ - y₁) * (log(x) - x₁) / (x₂ - x₁))
end

"""
Draw the prompt fission neutron spectra on the shared grid: ²⁵²Cf(sf) as the ratio of its
spectrum to a Maxwellian at the run's `T_M`, then every other system with a run as the ratio of
its spectrum per fission, `ν̄ N(E)`, to that of ²⁵²Cf(sf), which is what the measured ratios
are, the model taking its own multiplicities. The model is the black line and the measured
ratios are markers; a set tabulated as ²⁵²Cf over the system is inverted and a set the record
flags as superseded is skipped. Each panel carries the system and the model's ⟨E⟩ in the lower
left, below the data; the ²⁵²Cf panel adds the measured ⟨E⟩ and labels its unity guide with
`T_M`, the others add the multiplicity ratio the curve carries.
"""
function figure_spectra(drawable::Vector{Case})
    reference = findfirst(case -> case.label == REFERENCE.label, drawable)
    if reference === nothing
        @warn "no $(REFERENCE.label) run; the spectra are drawn against it" expected =
            joinpath(RUNS, REFERENCE.label)
        return nothing
    end
    reference_directory = run_directory(drawable[reference])
    reference_spectrum = read_spectrum(reference_directory)
    if reference_spectrum === nothing
        @warn "no spectrum to draw" expected =
            joinpath(reference_directory, "spectrum_vs_E.dat")
        return nothing
    end
    low, high = ENERGY_WINDOW

    # One entry per panel: the system, its run, the model ratio, the measured sets, the guide
    # label of the unity line, if any, and the annotation.
    panels = []
    let case = drawable[reference], directory = reference_directory
        E, N = reference_spectrum
        T_M = maxwellian_temperature(directory)
        mean_E = get(read_totals(directory), "mean_E", NaN)
        @printf("  ⟨E⟩ = %.4f MeV against a Maxwellian at T_M = %.2f MeV\n", mean_E, T_M)
        inside = low .<= E .<= high
        ratio = N[inside] ./ maxwellian.(E[inside], T_M)
        sets = load_measured(case.label, "spectrum_maxwellian_ratio_vs_E", ENERGY_WINDOW)
        push!(
            panels,
            (;
                case,
                directory,
                E = E[inside],
                ratio,
                sets,
                ylabel = L"N(E)\,/\,M(E,\,T_\mathrm{M})",
                guide = MathText(@sprintf("\$T_\\mathrm{M} = %.2f\$ MeV", T_M)),
                lines = MathText[
                    case.notation,
                    MathText(
                        @sprintf(
                            "\$\\langle E \\rangle = %.3f\$ MeV, measured %.2f MeV",
                            mean_E,
                            MEASURED_MEAN_ENERGY
                        )
                    ),
                ],
            ),
        )
    end
    for case in drawable
        case.label == REFERENCE.label && continue
        directory = run_directory(case)
        spectrum = read_spectrum(directory)
        if spectrum === nothing
            @warn "no spectrum to draw" expected = joinpath(directory, "spectrum_vs_E.dat")
            continue
        end
        E, N = spectrum
        E_ref, N_ref = reference_spectrum
        inside = low .<= E .<= high
        denominator =
            E == E_ref ? N_ref[inside] :
            loglog_interpolate.(Ref(E_ref), Ref(N_ref), E[inside])
        # The measured ratios are of the spectra per fission, ν̄ N(E), so the model's are too,
        # with its own multiplicities: nothing is scaled to the measurement.
        totals = read_totals(directory)
        mean_E = get(totals, "mean_E", NaN)
        multiplicity_ratio =
            get(totals, "nu_per_fission", NaN) /
            get(read_totals(reference_directory), "nu_per_fission", NaN)
        ratio = multiplicity_ratio .* N[inside] ./ denominator
        @printf(
            "  %-10s ⟨E⟩ = %.4f MeV, ν̄/ν̄(%s) = %.3f, drawn as the per-fission ratio\n",
            case.label,
            mean_E,
            REFERENCE.label,
            multiplicity_ratio
        )
        record = retrieval_record(case.label, "spectrum_cf252_ratio_vs_E")
        skip = [
            tag for (tag, entry) in record if
            any(startswith("superseded"), get(entry, "qualifiers", String[]))
        ]
        invert = [
            tag for (tag, entry) in record if
            get(entry, "ratio_orientation", "") == "reference_over_system"
        ]
        sets = load_measured(
            case.label,
            "spectrum_cf252_ratio_vs_E",
            ENERGY_WINDOW;
            skip = skip,
            invert = invert,
        )
        push!(
            panels,
            (;
                case,
                directory,
                E = E[inside],
                ratio,
                sets,
                ylabel = L"\bar{\nu}\,N(E)\,/\,\bar{\nu}_\mathrm{Cf}\,N_\mathrm{Cf}(E)",
                guide = nothing,
                lines = MathText[
                    case.notation,
                    MathText(@sprintf("\$\\langle E \\rangle = %.3f\$ MeV", mean_E)),
                    MathText(
                        @sprintf(
                            "\$\\bar{\\nu}/\\bar{\\nu}_\\mathrm{Cf} = %.3f\$",
                            multiplicity_ratio
                        )
                    ),
                ],
            ),
        )
    end

    xlimits = padded_log_limits(
        reduce(vcat, [vcat(p.E, (set.x for set in p.sets)...) for p in panels]),
    )
    # The annotation sits in the lower left, below the data: at the lowest energies every ratio
    # drawn is at or above unity.
    annotated = (
        first(xlimits),
        exp10(
            log10(first(xlimits)) +
            ANNOTATION_EXTENT * (log10(last(xlimits)) - log10(first(xlimits))),
        ),
    )

    rows = cld(length(panels), GRID_COLUMNS)
    banks = legend_rows(length(p.sets) + 1 for p in panels)
    figure = Figure(size = (GRID_WIDTH, GRID_ROW_HEIGHT * rows))
    provenance = Dict{String, Any}()
    for (index, p) in pairs(panels)
        row, column = fldmod1(index, GRID_COLUMNS)
        bottom = row == rows
        cell = GridLayout(figure[row, column])
        has_legend = !isempty(p.sets)
        slot = cell[has_legend ? 2 : 1, 1]
        axis = panel(
            slot;
            xscale = log10,
            xticks = SPECTRUM_TICKS,
            xlabel = L"Laboratory neutron energy $E$ [MeV]",
            ylabel = p.ylabel,
            xlabelvisible = bottom,
            xticklabelsvisible = bottom,
        )

        # Unity is a reference only for the ratio to the Maxwellian; a ratio of two spectra
        # per fission has no reason to be one.
        guided = p.guide !== nothing
        guided &&
            hlines!(axis, [1.0]; color = GUIDE_COLOUR, linestyle = :dash, linewidth = 1.5)
        begin_panel!()
        handles = Any[]
        labels = Any[]
        for set in p.sets
            push!(handles, draw_measured!(axis, set))
            push!(labels, set.label)
        end
        push!(handles, lines!(axis, p.E, p.ratio; color = :black, linewidth = 3))
        push!(labels, "This model")

        low_y, high_y = ordinate_limits(
            guided ? vcat(p.ratio, 1.0) : p.ratio,
            p.sets;
            cleared = annotated,
            bars = :lower,
        )
        ylimits = (low_y - length(p.lines) * ANNOTATION_BAND * (high_y - low_y), high_y)
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
        guided && label_guide!(slot, 1.0, p.guide, ylimits; side = :right)
        annotate!(axis, p.lines; corner = :lower_left)

        has_legend && panel_legend!(cell, handles, labels, banks)
        provenance[p.case.label] = system_record(p.directory, p.sets)
    end
    rowgap!(figure.layout, 12)
    colgap!(figure.layout, 12)

    export_figure(figure, "prompt_neutron_spectrum", provenance)
    return nothing
end

"""
Draw the post-neutron product mass yield `Y(A_p)` of every system whose run tabulated it, on
the shared grid: the model as the black line, the run's own pre-neutron `Y(A)` as a dotted
guide, the evaluated independent yields summed over charge as coloured lines with their
uncertainty as a band, and the retrieved post-neutron mass yields as markers — each `SEC` set
named, the chain yields `CHN` pooled into one muted series, which β decay, conserving `A`,
makes post-neutron mass yields up to the delayed-neutron transfer between chains. A set that a
qualifier in `MARKED_QUALIFIERS` marks is drawn with open markers. The annotation band above the
data carries the system; each panel takes its own mass range, the run's yield window
(`yield_window`).
"""
function figure_product_yield(drawable::Vector{Case})
    panels = []
    for case in drawable
        directory = run_directory(case)
        path = joinpath(directory, "Y_vs_A_p.dat")
        isfile(path) || continue
        window = yield_window(directory)
        A_p, Y_p, _ = read_columns(path)
        A, Y, _ = read_columns(joinpath(directory, "Y_vs_A.dat"))
        evaluated = []
        for library in EVALUATED_YIELDS
            file = joinpath(MEASURED, case.label, "Y_vs_A_p", library.file)
            isfile(file) || continue
            x, y, σ = read_columns(file)
            keep = first(window) .<= x .<= last(window)
            push!(evaluated, (library, file, x[keep], y[keep], σ[keep]))
        end
        products, excluded = load_products(case.label, "Y_vs_A_p")
        for set in products
            set.kind in ("SEC", "CHN") || throw(
                ArgumentError(
                    "$(set.path): a post-neutron mass yield coded $(set.kind), not SEC or CHN",
                ),
            )
        end
        named = named_series(filter(set -> set.kind == "SEC", products), window)
        pool = pooled_series("", filter(set -> set.kind == "CHN", products), window)
        n = length(pool.paths)
        chain = ProductSeries(
            "Chain yields ($n set$(n == 1 ? "" : "s"))",
            pool.paths,
            pool.x,
            pool.y,
            pool.σ,
            pool.marked,
        )
        push!(
            panels,
            (; case, directory, A_p, Y_p, A, Y, evaluated, named, chain, excluded, window),
        )
    end
    isempty(panels) && return nothing

    # Entries per legend: the guide, the evaluations, the named sets, the chain yields where
    # any are drawn, the key to the open markers where any set is marked, and the model.
    entries(p) =
        2 +
        length(p.evaluated) +
        length(p.named) +
        !isempty(p.chain.x) +
        any(s -> any(s.marked), [p.named; p.chain])
    rows = cld(length(panels), GRID_COLUMNS)
    banks = legend_rows(entries(p) for p in panels)
    # Each legend row beyond the two of the house grid lengthens the row by its own height, so
    # that the axes keep the height of the other grid figures.
    height = GRID_ROW_HEIGHT + LEGEND_ROW_HEIGHT * max(banks - 2, 0)
    # The last tick label may reach past the right spine; the right padding holds it.
    figure = Figure(size = (GRID_WIDTH, height * rows), figure_padding = (10, 24, 10, 10))
    provenance = Dict{String, Any}()
    for (index, p) in pairs(panels)
        row, column = fldmod1(index, GRID_COLUMNS)
        legend_slot, slot = grid_slots(figure, index)
        xlimits = padded_limits([first(p.window), last(p.window)])
        axis = panel(
            slot;
            ylabel = L"Y(A_p)\ [\%]",
            xlabel = L"Product mass $A_p$",
            ylabelvisible = column == 1,
            xticks = mass_ticks(xlimits),
            yticks = 0:2:12,
        )
        begin_panel!()
        handles = Any[]
        labels = Any[]
        # The pre-neutron yield first, beneath everything: it is the guide the shift is read
        # against.
        push!(
            handles,
            lines!(axis, p.A, p.Y; color = GUIDE_COLOUR, linestyle = :dot, linewidth = 1.5),
        )
        push!(labels, L"Pre-neutron $Y(A)$")
        ordinates = vcat(p.Y_p, p.Y)
        for (library, _, x, y, σ) in p.evaluated
            any(>(0), σ) && band!(axis, x, y .- σ, y .+ σ; color = (library.colour, 0.2))
            push!(
                handles,
                lines!(
                    axis,
                    x,
                    y;
                    color = library.colour,
                    linestyle = library.linestyle,
                    linewidth = 2.5,
                ),
            )
            push!(labels, library.label)
            append!(ordinates, y .+ σ)
        end
        # The chain yields beneath the named sets, muted: they are many, and partial chains.
        measured = ProductSeries[]
        if !isempty(p.chain.x)
            draw_series!(
                axis,
                p.chain,
                CHAIN_COLOUR,
                :circle;
                markersize = CHAIN_MARKER_SIZE,
            )
            push!(measured, p.chain)
        end
        for s in p.named
            colour, marker = dataset_style(PRODUCT_PALETTE, s.label)
            push!(
                handles,
                draw_series!(axis, s, colour, marker; markersize = DENSE_MARKER_SIZE),
            )
            push!(labels, s.label)
            push!(measured, s)
        end
        if !isempty(p.chain.x)
            push!(
                handles,
                MarkerElement(;
                    marker = :circle,
                    markersize = CHAIN_MARKER_SIZE,
                    color = CHAIN_COLOUR,
                    strokecolor = darker(CHAIN_COLOUR),
                    strokewidth = 1.5,
                ),
            )
            push!(labels, p.chain.label)
        end
        if any(s -> any(s.marked), measured)
            push!(handles, marked_key(; markersize = DENSE_MARKER_SIZE))
            push!(labels, MARKED_LABEL)
        end
        push!(handles, lines!(axis, p.A_p, p.Y_p; color = :black, linewidth = 3))
        push!(labels, "This model")
        for s in measured
            append!(ordinates, s.y)
        end

        # One line of annotation above all the data: both humps reach the top, so no corner is
        # free, and the band is reserved over the whole width.
        low, high = padded_limits(ordinates)
        ylimits = (max(low, 0.0), high + ANNOTATION_BAND * (high - low))
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
        annotate!(axis, MathText[p.case.notation]; corner = :upper_left)

        grid_legend!(legend_slot, handles, labels, banks)
        provenance[p.case.label] = system_record(
            p.directory,
            measured;
            evaluated = String[e[2] for e in p.evaluated],
            excluded = p.excluded,
        )
    end
    grid_gaps!(figure, rows)

    export_figure(figure, "product_mass_yield", provenance)
    return nothing
end

"""
    read_independent_yields(path) -> Dict{Tuple{Int, Int}, NTuple{2, Float64}}

The independent yields of a table `Z A_p Y [Y_uncertainty]` — a run's `Y_vs_Z_A_p.dat` or a
staged evaluation — keyed by `(Z, A_p)`, each with its uncertainty, zero where none is stated.
Columns are taken by position.
"""
function read_independent_yields(path::AbstractString)
    yields = Dict{Tuple{Int, Int}, NTuple{2, Float64}}()
    for row in read_rows(path)
        length(row) >= 3 || continue
        σ = length(row) >= 4 && isfinite(row[4]) ? row[4] : 0.0
        yields[(round(Int, row[1]), round(Int, row[2]))] = (row[3], σ)
    end
    return yields
end

"""
    grouped_yield(yields, key) -> (x, Y, σ)

The independent yields summed over the products that share `key(Z, A_p)`, in order of the key,
the uncertainties added in quadrature: the evaluations state no correlations.
"""
function grouped_yield(yields::AbstractDict, key::Function)
    sums = Dict{Int, NTuple{2, Float64}}()
    for ((Z, A_p), (Y, σ)) in yields
        k = key(Z, A_p)
        total, variance = get(sums, k, (0.0, 0.0))
        sums[k] = (total + Y, variance + σ^2)
    end
    ordered = sort!(collect(keys(sums)))
    return (
        Float64.(ordered),
        [first(sums[k]) for k in ordered],
        [sqrt(last(sums[k])) for k in ordered],
    )
end

"The elemental yield `Y(Z) = Σ_{A_p} Y(Z, A_p)`."
elemental_yield(yields::AbstractDict) = grouped_yield(yields, (Z, _) -> Z)

"The yield against product neutron number, `Y(N_p) = Σ_{A_p − Z = N_p} Y(Z, A_p)`."
neutron_number_yield(yields::AbstractDict) = grouped_yield(yields, (Z, A_p) -> A_p - Z)

"The isotopic distribution `Y(A_p)` of the element of charge `Z`."
isotopic_yield(yields::AbstractDict, Z::Integer) =
    grouped_yield(filter(entry -> first(first(entry)) == Z, yields), (_, A_p) -> A_p)

"""
The yield of the heavy products against product mass: those of charge above `Z₀/2`, and half of
those at `Z₀/2`. The charge of a product is that of its fragment, which places it on the heavy
side of its split as its mass does in all but the most nearly symmetric splits.
"""
function heavy_product_yield(yields::AbstractDict, compound_charge::Integer)
    heavy = Dict{Int, Float64}()
    for ((Z, A_p), (Y, _)) in yields
        share = 2Z > compound_charge ? 1.0 : 2Z == compound_charge ? 0.5 : 0.0
        heavy[A_p] = get(heavy, A_p, 0.0) + share * Y
    end
    return heavy
end

"""
The global proton even–odd effect of an elemental yield, `δ_p = (Y_e − Y_o)/(Y_e + Y_o)`, `Y_e`
and `Y_o` the yields summed over the even and the odd charges.
"""
function even_odd_effect(Z::Vector{Float64}, Y::Vector{Float64})
    even = sum(Y[i] for i in eachindex(Z) if iseven(round(Int, Z[i])); init = 0.0)
    odd = sum(Y[i] for i in eachindex(Z) if isodd(round(Int, Z[i])); init = 0.0)
    return (even - odd) / (even + odd)
end

"""
The evaluated independent yields staged for one system, `(library, path)` for every library of
`EVALUATED_YIELDS` whose table exists under `<data>/<system>/<observable>/`.
"""
function evaluated_files(system::AbstractString, observable::AbstractString)
    files = Tuple{eltype(EVALUATED_YIELDS), String}[]
    for library in EVALUATED_YIELDS
        path = joinpath(MEASURED, system, observable, library.file)
        if isfile(path)
            push!(files, (library, path))
        else
            @warn "no $(library.label) independent yields to draw" expected = path
        end
    end
    return files
end

"""
Decade ticks inside `limits`, labelled as plain decimals — 0.001, 0.01, 0.1, 1, 10 — for a
logarithmic axis spanning a few decades. Given the range of the data rather than of the axis,
they leave the annotation band above the data unlabelled.
"""
function decade_ticks(limits::Tuple{Real, Real})
    exponents = ceil(Int, log10(first(limits))):floor(Int, log10(last(limits)))
    labels = [k < 0 ? "0." * "0"^(-k - 1) * "1" : "1" * "0"^k for k in exponents]
    return (exp10.(exponents), labels)
end

"""
Logarithmic ordinate limits from `floor` up to `top`, with the fraction `reserved` of the axis
height left clear above the data for annotation.
"""
function log_limits(floor::Real, top::Real, reserved::Real)
    low, high = log10(floor), log10(top)
    return (Float64(floor), exp10(low + (high - low) / (1 - reserved)))
end

"""
Linear ordinate limits spanning `values`, with a band of `lines` lines of annotation reserved
above the data (`corner = :upper_left` or `:upper_right`) or below it (`:lower_left`).
"""
function banded_limits(values, lines::Integer, corner::Symbol)
    low, high = padded_limits(values)
    band = lines * ANNOTATION_BAND * (high - low)
    corner in (:upper_left, :upper_right) && return (low, high + band)
    corner === :lower_left && return (low - band, high)
    throw(
        ArgumentError(
            "corner must be :upper_left, :upper_right or :lower_left, got :$corner",
        ),
    )
end

"""
Integer mass ticks inside `limits`: every 10 over a full yield window, every 4 over a range of
some tens of masses, every 2 over a short one, so that the three-digit labels stay apart.
"""
function mass_ticks(limits::Tuple{Real, Real})
    span = last(limits) - first(limits)
    step = span > 40 ? 10 : span > 8 ? 4 : 2
    return (ceil(Int, first(limits) / step) * step):step:floor(Int, last(limits))
end

"""
    input_kinetic_energy_by_mass(directory, case) -> (A, TKE, path) or nothing

The pre-neutron mean total kinetic energy `⟨TKE⟩(A)` of the yield a run was solved with, read
from the input its `configuration.toml` names under `[yield]`: the mean over the measured matrix
`mass_energy_yield_file`, or the measured `mean_kinetic_energy_file`. Where the run imposed the
identity of the two fragments of a split, `symmetrize`, the complement is averaged in as the
package does it, cell by cell for a matrix and mass by mass for a mean. `nothing` where the run
was solved without a yield or its input is not under the data directory.
"""
function input_kinetic_energy_by_mass(directory::AbstractString, case::Case)
    path = joinpath(directory, "configuration.toml")
    isfile(path) || return nothing
    section = get(TOML.parsefile(path), "yield", nothing)
    section isa AbstractDict || return nothing
    symmetrize = get(section, "symmetrize", true)
    A₀ = case.compound_mass
    if haskey(section, "mass_energy_yield_file")
        file = joinpath(MEASURED, section["mass_energy_yield_file"])
        isfile(file) || return nothing
        cells = Dict{Tuple{Int, Float64}, Float64}()
        for row in read_rows(file)
            length(row) >= 3 && isfinite(row[3]) || continue
            cells[(round(Int, row[1]), row[2])] = row[3]
        end
        if symmetrize
            paired = Dict{Tuple{Int, Float64}, Float64}()
            for ((A, TKE), Y) in cells
                partner = get(cells, (A₀ - A, TKE), nothing)
                paired[(A, TKE)] =
                    paired[(A₀ - A, TKE)] = partner === nothing ? Y : (Y + partner) / 2
            end
            cells = paired
        end
        sums = Dict{Int, NTuple{2, Float64}}()
        for ((A, TKE), Y) in cells
            weighted, norm = get(sums, A, (0.0, 0.0))
            sums[A] = (weighted + Y * TKE, norm + Y)
        end
        masses = sort!([A for (A, (_, norm)) in sums if norm > 0])
        return (Float64.(masses), [first(sums[A]) / last(sums[A]) for A in masses], file)
    elseif haskey(section, "mean_kinetic_energy_file")
        file = joinpath(MEASURED, section["mean_kinetic_energy_file"])
        isfile(file) || return nothing
        A, TKE, _ = read_columns(file)
        means = Dict(round(Int, a) => t for (a, t) in zip(A, TKE))
        if symmetrize
            paired = Dict{Int, Float64}()
            for (a, t) in means
                partner = get(means, A₀ - a, nothing)
                paired[a] = paired[A₀ - a] = partner === nothing ? t : (t + partner) / 2
            end
            means = paired
        end
        masses = sort!(collect(keys(means)))
        return (Float64.(masses), [means[a] for a in masses], file)
    end
    return nothing
end

"""
    figure_system_grid(name, panels, xlimits, draw!; xlabel, ylabel, legend = true, axis...)

Lay out one panel per system on the shared grid of the comparison figures — two columns, the
legend of each panel above it unless `legend` is false, one abscissa range, its label and tick
labels on the bottom row only — and export it. `draw!(axis, panel, slot)` draws one panel and
returns `handles`, `labels`, `ylimits`, the `lines` of annotation and their `corner`, and the
provenance `record`; `slot` is the layout cell of the axis, for guides labelled outside the
frame. The legends of one grid take the same number of rows, so the axes of one row are the
same height.
"""
function figure_system_grid(
    name::AbstractString,
    panels::AbstractVector,
    xlimits::Tuple{Real, Real},
    draw!::Function;
    xlabel,
    ylabel,
    legend::Bool = true,
    axis_options...,
)
    isempty(panels) && return nothing
    rows = cld(length(panels), GRID_COLUMNS)
    figure = Figure(size = (GRID_WIDTH, GRID_ROW_HEIGHT * rows))
    provenance = Dict{String, Any}()
    legends = []
    for (index, p) in pairs(panels)
        row, column = fldmod1(index, GRID_COLUMNS)
        bottom = row == rows
        legend_slot, slot = grid_slots(figure, index; legend = legend)
        axis = panel(
            slot;
            xlabel = xlabel,
            ylabel = ylabel,
            ylabelvisible = column == 1,
            xlabelvisible = bottom,
            xticklabelsvisible = bottom,
            axis_options...,
        )
        begin_panel!()
        drawn = draw!(axis, p, slot)
        xlims!(axis, xlimits...)
        ylims!(axis, drawn.ylimits...)
        annotate!(axis, drawn.lines; corner = drawn.corner)
        push!(legends, (legend_slot, drawn.handles, drawn.labels))
        provenance[p.case.label] = drawn.record
    end
    if legend
        banks = legend_rows(length(handles) for (_, handles, _) in legends)
        for (position, handles, labels) in legends
            grid_legend!(position, handles, labels, banks)
        end
        # Legend rows beyond the two of the house grid lengthen each grid row by their height,
        # so that the axes keep the height of the other grid figures.
        banks > 2 && resize!(
            figure,
            GRID_WIDTH,
            (GRID_ROW_HEIGHT + LEGEND_ROW_HEIGHT * (banks - 2)) * rows,
        )
    end
    grid_gaps!(figure, rows; legend = legend)

    export_figure(figure, name, provenance)
    return nothing
end

"""
The run of one system and the independent yields it wrote, beside the evaluated ones staged for
it, or `nothing` where the run wrote none.
"""
function independent_yield_panel(case::Case)
    directory = run_directory(case)
    path = joinpath(directory, "Y_vs_Z_A_p.dat")
    if !isfile(path)
        @warn "no Y_vs_Z_A_p.dat to draw" directory
        return nothing
    end
    evaluated = [
        (library, file, read_independent_yields(file)) for
        (library, file) in evaluated_files(case.label, "Y_vs_Z_A_p")
    ]
    return (; case, directory, model = read_independent_yields(path), evaluated)
end

evaluated_paths(p) = String[file for (_, file, _) in p.evaluated]

"""
Draw the elemental yield `Y(Z)` of every system whose run tabulated its independent yields, on
the shared grid and a logarithmic axis: the model as the black line, the evaluated independent
yields summed over mass as coloured lines. Each legend entry carries the global proton even–odd
effect `δ_p` of its curve.
"""
function figure_charge_yield(drawable::Vector{Case})
    panels = []
    for case in drawable
        p = independent_yield_panel(case)
        p === nothing && continue
        model = elemental_yield(p.model)
        evaluated =
            [(library, elemental_yield(yields)) for (library, _, yields) in p.evaluated]
        push!(panels, (; p..., curves = [(nothing, model); evaluated]))
    end
    isempty(panels) && return nothing

    charges = Float64[]
    top = 0.0
    for p in panels, (_, (Z, Y, _)) in p.curves
        append!(charges, Z[Y .>= YIELD_FLOOR])
        top = max(top, maximum(Y))
    end
    xlimits = padded_limits(charges)
    # Three lines of annotation take a third of a grid panel's height.
    ylimits = log_limits(YIELD_FLOOR, top, 0.36)

    function draw!(axis::Axis, p, _)
        handles = Any[]
        labels = Any[]
        effects = String[]
        # The evaluations first, the model on top of them.
        for (library, (Z, Y, _)) in [p.curves[2:end]; p.curves[1:1]]
            keep = Y .> 0
            δ = @sprintf("%.0f %%", 100 * even_odd_effect(Z, Y))
            if library === nothing
                push!(
                    handles,
                    lines!(axis, Z[keep], Y[keep]; color = :black, linewidth = 3),
                )
                push!(labels, "This model")
                pushfirst!(effects, "model $δ")
            else
                push!(
                    handles,
                    lines!(
                        axis,
                        Z[keep],
                        Y[keep];
                        color = library.colour,
                        linestyle = library.linestyle,
                        linewidth = 2.5,
                    ),
                )
                push!(labels, library.label)
                push!(effects, "$(library.label) $δ")
            end
        end
        # A string without mathematics is set as text, not as a formula.
        lines = Any[
            p.case.notation,
            MathText("Even–odd effect \$\\delta_p\$: $(first(effects))"),
            join(effects[2:end], ", "),
        ]
        return (;
            handles,
            labels,
            ylimits,
            lines,
            corner = :upper_left,
            record = system_record(
                p.directory,
                MeasuredSet[];
                evaluated = evaluated_paths(p),
            ),
        )
    end
    figure_system_grid(
        "product_charge_yield",
        panels,
        xlimits,
        draw!;
        xlabel = L"Product charge $Z$",
        ylabel = L"Y(Z)\ [\%]",
        yscale = log10,
        yticks = decade_ticks((YIELD_FLOOR, top)),
        xticks = 20:10:80,
    )
    return nothing
end

"""
Draw the yield against product neutron number, `Y(N_p)`, of every system whose run tabulated its
independent yields, on the shared grid: the model as the black line, the evaluated independent
yields regrouped by `N_p = A_p − Z` as coloured lines with their uncertainty as a band, and the
closed neutron shells as labelled guides.
"""
function figure_neutron_number_yield(drawable::Vector{Case})
    panels = []
    for case in drawable
        p = independent_yield_panel(case)
        p === nothing && continue
        model = neutron_number_yield(p.model)
        evaluated = [
            (library, neutron_number_yield(yields)) for (library, _, yields) in p.evaluated
        ]
        push!(panels, (; p..., model_curve = model, evaluated_curves = evaluated))
    end
    isempty(panels) && return nothing

    neutrons = Float64[]
    for p in panels, (N, Y, _) in [p.model_curve; last.(p.evaluated_curves)]
        append!(neutrons, N[Y .>= MEAN_FLOOR])
    end
    xlimits = padded_limits(neutrons)

    function draw!(axis::Axis, p, slot)
        for shell in NEUTRON_SHELLS
            first(xlimits) < shell < last(xlimits) || continue
            vlines!(axis, [shell]; color = GUIDE_COLOUR, linestyle = :dash, linewidth = 1.5)
            label_guide!(slot, shell, MathText("\$N = $(shell)\$"), xlimits; side = :top)
        end
        handles = Any[]
        labels = Any[]
        ordinates = copy(p.model_curve[2])
        for (library, (N, Y, σ)) in p.evaluated_curves
            keep = first(xlimits) .<= N .<= last(xlimits)
            N, Y, σ = N[keep], Y[keep], σ[keep]
            any(>(0), σ) && band!(axis, N, Y .- σ, Y .+ σ; color = (library.colour, 0.2))
            push!(
                handles,
                lines!(
                    axis,
                    N,
                    Y;
                    color = library.colour,
                    linestyle = library.linestyle,
                    linewidth = 2.5,
                ),
            )
            push!(labels, library.label)
            append!(ordinates, Y .+ σ)
        end
        N, Y, _ = p.model_curve
        keep = first(xlimits) .<= N .<= last(xlimits)
        push!(handles, lines!(axis, N[keep], Y[keep]; color = :black, linewidth = 3))
        push!(labels, "This model")
        # The guide at N = 50 crosses the upper left; the upper right, beyond the heavy hump,
        # is clear.
        low, high = banded_limits(ordinates, 1, :upper_right)
        return (;
            handles,
            labels,
            ylimits = (max(low, 0.0), high),
            lines = MathText[p.case.notation],
            corner = :upper_right,
            record = system_record(
                p.directory,
                MeasuredSet[];
                evaluated = evaluated_paths(p),
            ),
        )
    end
    figure_system_grid(
        "product_neutron_number_yield",
        panels,
        xlimits,
        draw!;
        xlabel = L"Product neutron number $N_p$",
        ylabel = L"Y(N_p)\ [\%]",
        xticks = 40:10:110,
    )
    return nothing
end

"""
    kinetic_energy_series(sets, window) -> Vector{ProductSeries}

The retrieved product kinetic energies of one system in `window`, one series per author–year
label. Where there are more series than `PRODUCT_PALETTE` has colour and marker pairs, the
`NAMED_KINETIC_ENERGY_SETS` with the most points keep their names and the rest are pooled into
one, labelled with their number.
"""
function kinetic_energy_series(sets::Vector{ProductSet}, window)
    series = named_series(sets, window)
    length(series) <= length(PRODUCT_PALETTE.colours) && return series
    order = sortperm(series; by = s -> length(s.x), rev = true)
    named = series[sort(order[1:NAMED_KINETIC_ENERGY_SETS])]
    rest = series[order[(NAMED_KINETIC_ENERGY_SETS + 1):end]]
    pool = ProductSeries(
        "Other sets ($(length(rest)))",
        reduce(vcat, (s.paths for s in rest)),
        reduce(vcat, (s.x for s in rest)),
        reduce(vcat, (s.y for s in rest)),
        reduce(vcat, (s.σ for s in rest)),
        reduce(vcat, (s.marked for s in rest)),
    )
    return [named; pool]
end

"""
Draw the mean kinetic energy of the products against product mass, `E_K_p(A_p)`, of every system
whose run tabulated it, on the shared grid: the model as the black line, the pre-neutron
`E_K(A) = ⟨TKE⟩(A) (A₀ − A)/A₀` of the yield the run was solved with as a dashed guide, and the
retrieved post-neutron kinetic energies as markers, as [`kinetic_energy_series`](@ref) groups
them, open where a qualifier in `MARKED_QUALIFIERS` marks the set. Masses whose product yield is
below `MEAN_FLOOR` are not drawn.
"""
function figure_product_kinetic_energy(drawable::Vector{Case})
    panels = []
    for case in drawable
        directory = run_directory(case)
        path = joinpath(directory, "E_K_p_vs_A_p.dat")
        if !isfile(path)
            @warn "no E_K_p_vs_A_p.dat to draw" directory
            continue
        end
        A_p, E_K_p, _ = read_columns(path)
        masses, yields, _ = read_columns(joinpath(directory, "Y_vs_A_p.dat"))
        drawn = Set(masses[yields .>= MEAN_FLOOR])
        keep = [a in drawn for a in A_p]
        A_p, E_K_p = A_p[keep], E_K_p[keep]
        window = extrema(A_p)
        guide = input_kinetic_energy_by_mass(directory, case)
        if guide !== nothing
            A, TKE, file = guide
            inside = first(window) .<= A .<= last(window)
            guide = (
                A[inside],
                TKE[inside] .* (case.compound_mass .- A[inside]) ./ case.compound_mass,
                file,
            )
        end
        products, excluded = load_products(case.label, "E_K_p_vs_A_p")
        series = kinetic_energy_series(products, window)
        push!(panels, (; case, directory, A_p, E_K_p, guide, series, excluded))
    end
    isempty(panels) && return nothing
    xlimits = padded_limits(reduce(vcat, [p.A_p for p in panels]))

    function draw!(axis::Axis, p, _)
        handles = Any[]
        labels = Any[]
        ordinates = copy(p.E_K_p)
        inputs = String[]
        if p.guide !== nothing
            A, E_K, file = p.guide
            push!(
                handles,
                lines!(
                    axis,
                    A,
                    E_K;
                    color = GUIDE_COLOUR,
                    linestyle = :dash,
                    linewidth = 1.5,
                ),
            )
            push!(labels, L"Pre-neutron $E_\mathrm{K}(A)$")
            append!(ordinates, E_K)
            push!(inputs, file)
        end
        for s in p.series
            pooled = startswith(s.label, "Other sets")
            colour, marker =
                pooled ? (CHAIN_COLOUR, :circle) : dataset_style(PRODUCT_PALETTE, s.label)
            push!(
                handles,
                draw_series!(axis, s, colour, marker; markersize = DENSE_MARKER_SIZE),
            )
            push!(labels, s.label)
            append!(ordinates, s.y)
        end
        if any(s -> any(s.marked), p.series)
            push!(handles, marked_key(; markersize = DENSE_MARKER_SIZE))
            push!(labels, MARKED_LABEL)
        end
        push!(handles, lines!(axis, p.A_p, p.E_K_p; color = :black, linewidth = 3))
        push!(labels, "This model")
        # The light products are the fast ones: the lower left is clear.
        lines = Any[p.case.notation]
        isempty(p.series) && push!(lines, "No measured set in EXFOR")
        return (;
            handles,
            labels,
            ylimits = banded_limits(ordinates, length(lines), :lower_left),
            lines,
            corner = :lower_left,
            record = system_record(
                p.directory,
                p.series;
                inputs = inputs,
                excluded = p.excluded,
            ),
        )
    end
    figure_system_grid(
        "product_kinetic_energy",
        panels,
        xlimits,
        draw!;
        xlabel = L"Product mass $A_p$",
        ylabel = L"E_{\mathrm{K},p}(A_p)\ [\mathrm{MeV}]",
        xticks = 80:20:180,
    )
    return nothing
end

"""
Draw the mean kinetic energy of the products against their charge, `E_K_p(Z)`, of every system
whose run tabulated it, on the shared grid, the model as the black line. Charges whose elemental
yield is below `MEAN_FLOOR` are not drawn.
"""
function figure_product_kinetic_energy_by_charge(drawable::Vector{Case})
    panels = []
    for case in drawable
        directory = run_directory(case)
        path = joinpath(directory, "E_K_p_vs_Z.dat")
        independent = joinpath(directory, "Y_vs_Z_A_p.dat")
        if !isfile(path) || !isfile(independent)
            @warn "no E_K_p_vs_Z.dat or Y_vs_Z_A_p.dat to draw" directory
            continue
        end
        Z, E_K_p, _ = read_columns(path)
        charges, yields, _ = elemental_yield(read_independent_yields(independent))
        drawn = Set(charges[yields .>= MEAN_FLOOR])
        keep = [z in drawn for z in Z]
        push!(panels, (; case, directory, Z = Z[keep], E_K_p = E_K_p[keep]))
    end
    isempty(panels) && return nothing
    xlimits = padded_limits(reduce(vcat, [p.Z for p in panels]))

    function draw!(axis::Axis, p, _)
        handle = lines!(axis, p.Z, p.E_K_p; color = :black, linewidth = 3)
        return (;
            handles = Any[handle],
            labels = Any["This model"],
            ylimits = banded_limits(p.E_K_p, 1, :lower_left),
            lines = MathText[p.case.notation],
            corner = :lower_left,
            record = system_record(p.directory, MeasuredSet[]),
        )
    end
    figure_system_grid(
        "product_kinetic_energy_by_charge",
        panels,
        xlimits,
        draw!;
        xlabel = L"Product charge $Z$",
        ylabel = L"E_{\mathrm{K},p}(Z)\ [\mathrm{MeV}]",
        # One curve per panel needs no legend.
        legend = false,
        xticks = 20:10:80,
    )
    return nothing
end

"""
Draw the mean post-neutron total kinetic energy against the heavy product mass, `TKE_p(A_p)`, of
every system whose run tabulated it, on the shared grid: the model as the black line and the
pre-neutron `⟨TKE⟩(A)` of the yield the run was solved with as a dashed guide, over the heavy
masses. Masses whose heavy-product yield is below `MEAN_FLOOR` are not drawn. The annotation
states the means of both over the run.
"""
function figure_product_total_kinetic_energy(drawable::Vector{Case})
    panels = []
    for case in drawable
        directory = run_directory(case)
        path = joinpath(directory, "TKE_p_vs_A_p.dat")
        independent = joinpath(directory, "Y_vs_Z_A_p.dat")
        if !isfile(path) || !isfile(independent)
            @warn "no TKE_p_vs_A_p.dat or Y_vs_Z_A_p.dat to draw" directory
            continue
        end
        A_p, TKE_p, _ = read_columns(path)
        heavy =
            heavy_product_yield(read_independent_yields(independent), case.compound_charge)
        keep = [get(heavy, round(Int, a), 0.0) >= MEAN_FLOOR for a in A_p]
        A_p, TKE_p = A_p[keep], TKE_p[keep]
        guide = input_kinetic_energy_by_mass(directory, case)
        if guide !== nothing
            A, TKE, file = guide
            inside = (2 .* A .>= case.compound_mass) .& (A .<= last(A_p))
            guide = (A[inside], TKE[inside], file)
        end
        push!(panels, (; case, directory, A_p, TKE_p, guide))
    end
    isempty(panels) && return nothing
    masses = Float64[]
    for p in panels
        append!(masses, p.A_p)
        p.guide === nothing || append!(masses, p.guide[1])
    end
    xlimits = padded_limits(masses)

    function draw!(axis::Axis, p, _)
        handles = Any[]
        labels = Any[]
        ordinates = copy(p.TKE_p)
        inputs = String[]
        if p.guide !== nothing
            A, TKE, file = p.guide
            push!(
                handles,
                lines!(
                    axis,
                    A,
                    TKE;
                    color = GUIDE_COLOUR,
                    linestyle = :dash,
                    linewidth = 1.5,
                ),
            )
            push!(labels, L"Pre-neutron $\langle \mathrm{TKE} \rangle(A)$")
            append!(ordinates, TKE)
            push!(inputs, file)
        end
        push!(handles, lines!(axis, p.A_p, p.TKE_p; color = :black, linewidth = 3))
        push!(labels, "This model")
        totals = read_totals(p.directory)
        lines = MathText[
            p.case.notation,
            MathText(
                @sprintf(
                    "\$\\langle \\mathrm{TKE}_p \\rangle = %.1f\$ MeV",
                    get(totals, "mean_TKE_p", NaN)
                )
            ),
            MathText(
                @sprintf(
                    "\$\\langle \\mathrm{TKE} \\rangle = %.1f\$ MeV",
                    get(totals, "mean_TKE", NaN)
                )
            ),
        ]
        # Near symmetry the heavy products are the slow ones: the upper left is clear.
        return (;
            handles,
            labels,
            ylimits = banded_limits(ordinates, length(lines), :upper_left),
            lines,
            corner = :upper_left,
            record = system_record(p.directory, MeasuredSet[]; inputs = inputs),
        )
    end
    figure_system_grid(
        "product_total_kinetic_energy",
        panels,
        xlimits,
        draw!;
        xlabel = L"Heavy product mass $A_p$",
        ylabel = L"\mathrm{TKE}_p(A_p)\ [\mathrm{MeV}]",
        xticks = 110:10:180,
    )
    return nothing
end

"""
    independent_series(sets, Z) -> ProductSeries

The retrieved independent yields of charge `Z` at or above `YIELD_FLOOR`, every set pooled into
one series; isomeric states as each record states them.
"""
function independent_series(sets::Vector{ProductSet}, Z::Integer)
    paths = String[]
    x = Float64[]
    y = Float64[]
    σ = Float64[]
    marked = Bool[]
    for set in sets
        rows = (set.Z .== Z) .& (set.y .>= YIELD_FLOOR)
        any(rows) || continue
        push!(paths, set.path)
        append!(x, set.A_p[rows])
        append!(y, set.y[rows])
        append!(σ, set.σ[rows])
        append!(marked, fill(set.marked, count(rows)))
    end
    return ProductSeries(INDEPENDENT_LABEL, paths, x, y, σ, marked)
end

"""
Draw, for one system, the isotopic distributions `Y(A_p)` at fixed `Z` of the
`ELEMENTS_PER_SIDE` most abundant light and heavy elements as small multiples on a logarithmic
axis: the model as the black line, the evaluated independent yields as markers, and the
retrieved independent yields `IND`, all sets pooled, as a third marker family, open where a
qualifier in `MARKED_QUALIFIERS` marks the set. The elements are ranked by the mean of the
evaluated elemental yields, so the panels do not move with the model; by the model's own where
no evaluation is staged.
"""
function figure_independent_yield(case::Case)
    p = independent_yield_panel(case)
    p === nothing && return nothing
    products, excluded = load_products(case.label, "Y_vs_Z_A_p")
    for set in products
        set.kind == "IND" || throw(
            ArgumentError("$(set.path): an independent yield coded $(set.kind), not IND"),
        )
    end
    reference =
        isempty(p.evaluated) ? [p.model] : [yields for (_, _, yields) in p.evaluated]
    abundance = Dict{Int, Float64}()
    for yields in reference, ((Z, _), (Y, _)) in yields
        abundance[Z] = get(abundance, Z, 0.0) + Y / length(reference)
    end
    function most_abundant(side)
        ranked = sort(
            [Z for Z in keys(abundance) if side(Z)];
            by = Z -> abundance[Z],
            rev = true,
        )
        return sort(ranked[1:min(ELEMENTS_PER_SIDE, length(ranked))])
    end
    rows = (
        most_abundant(Z -> 2Z < case.compound_charge),
        most_abundant(Z -> 2Z > case.compound_charge),
    )

    curves = Dict{Int, Any}()
    top = 0.0
    for Z in Iterators.flatten(rows)
        model = isotopic_yield(p.model, Z)
        evaluated =
            [(library, isotopic_yield(yields, Z)) for (library, _, yields) in p.evaluated]
        measured = independent_series(products, Z)
        masses = copy(measured.x)
        top = max(top, maximum(measured.y; init = 0.0))
        for (A_p, Y, _) in [model; last.(evaluated)]
            append!(masses, A_p[Y .>= YIELD_FLOOR])
            top = max(top, maximum(Y; init = 0.0))
        end
        curves[Z] = (; model, evaluated, measured, xlimits = padded_limits(masses))
    end
    # Two lines of annotation above the data, and a margin below the floor so that a marker
    # at the floor is drawn whole.
    ylimits = log_limits(YIELD_FLOOR, top, 0.24)
    ylimits = (
        exp10(log10(first(ylimits)) - MARGIN * log10(last(ylimits) / first(ylimits))),
        last(ylimits),
    )

    columns = maximum(length, rows)
    figure = Figure(size = (GRID_WIDTH, 2 * 350 + 100))
    handles = Dict{String, Any}()
    measured_paths = Set{String}()
    for (row, charges) in pairs(rows), (column, Z) in pairs(charges)
        c = curves[Z]
        axis = panel(
            figure[row + 1, column];
            yscale = log10,
            yticks = decade_ticks((YIELD_FLOOR, top)),
            yticklabelsvisible = column == 1,
            xticks = mass_ticks(c.xlimits),
        )
        for (library, (A_p, Y, _)) in c.evaluated
            keep = Y .>= YIELD_FLOOR
            handles[library.label] = scatter!(
                axis,
                A_p[keep],
                Y[keep];
                color = library.colour,
                marker = library.marker,
                markersize = library.marker === :circle ? 14 : 10,
                strokecolor = darker(library.colour),
                strokewidth = 1.5,
            )
        end
        if !isempty(c.measured.x)
            handle = draw_series!(
                axis,
                c.measured,
                INDEPENDENT_COLOUR,
                INDEPENDENT_MARKER;
                floor = first(ylimits),
            )
            any(.!c.measured.marked) && (handles[c.measured.label] = handle)
            any(c.measured.marked) && (
                handles[MARKED_LABEL] = marked_key(;
                    marker = INDEPENDENT_MARKER,
                    colour = INDEPENDENT_COLOUR,
                )
            )
            union!(measured_paths, c.measured.paths)
        end
        A_p, Y, _ = c.model
        keep = Y .> 0
        handles["This model"] =
            lines!(axis, A_p[keep], Y[keep]; color = :black, linewidth = 3)
        xlims!(axis, c.xlimits...)
        ylims!(axis, ylimits...)
        lines = MathText[MathText("\$\\mathrm{$(ELEMENT_SYMBOLS[Z])}\$, \$Z = $Z\$")]
        row == 1 && column == 1 && pushfirst!(lines, case.notation)
        annotate!(axis, lines; corner = :upper_left)
    end
    # A family drawn only open in this figure is still named in the legend, by its filled marker.
    if any(!isempty(curves[Z].measured.x) for Z in keys(curves)) &&
       !haskey(handles, INDEPENDENT_LABEL)
        handles[INDEPENDENT_LABEL] = MarkerElement(;
            marker = INDEPENDENT_MARKER,
            markersize = 14,
            color = INDEPENDENT_COLOUR,
            strokecolor = darker(INDEPENDENT_COLOUR),
            strokewidth = 1.5,
        )
    end
    order = [
        [library.label for library in EVALUATED_YIELDS]
        INDEPENDENT_LABEL
        MARKED_LABEL
        "This model"
    ]
    shown = [label for label in order if haskey(handles, label)]
    # One row holds the model and the two evaluations; with the measured family the entries
    # outrun the width of the figure and take two.
    Legend(
        figure[1, 1:columns],
        [handles[label] for label in shown],
        shown;
        nbanks = length(shown) > 3 ? 2 : 1,
        tellwidth = false,
        tellheight = true,
        padding = (0, 0, 0, 0),
        colgap = 24,
        patchsize = (35, 20),
    )
    Label(figure[2:3, 0], L"Y(Z, A_p)\ [\%]"; rotation = π / 2, tellheight = false)
    Label(figure[4, 1:columns], L"Product mass $A_p$"; tellwidth = false)
    rowgap!(figure.layout, 10)
    colgap!(figure.layout, 12)

    # The sets drawn, those with a point in a panel shown; the others are retrieved and read.
    drawn = [(; paths = [path]) for path in sort!(collect(measured_paths))]
    record = system_record(
        p.directory,
        drawn;
        evaluated = evaluated_paths(p),
        excluded = excluded,
    )
    record["elements"] = collect(Iterators.flatten(rows))
    export_figure(
        figure,
        "product_independent_yield_$(case.label)",
        Dict{String, Any}(case.label => record),
    )
    return nothing
end

"""
Draw the quantities internal to the sequential treatment against fragment mass, one panel each
in a 2 × 2 grid, every system as a line in its fixed style. None has a measured counterpart, so
the comparison is between systems.
"""
function figure_residual_quantities(drawable::Vector{Case})
    isempty(drawable) && return nothing
    quantities = (
        ("E_star_vs_A.dat", L"E^*(A)\ [\mathrm{MeV}]"),
        ("E_r_vs_A.dat", L"\bar{E}_\mathrm{r}(A)\ [\mathrm{MeV}]"),
        ("T_vs_A.dat", L"\langle T \rangle(A)\ [\mathrm{MeV}]"),
        ("eps_vs_A.dat", L"\langle \varepsilon \rangle(A)\ [\mathrm{MeV}]"),
    )
    directories = Dict(case.label => run_directory(case) for case in drawable)
    tables = Dict{Tuple{Int, String}, Tuple{Vector{Float64}, Vector{Float64}}}()
    for (index, (name, _)) in pairs(quantities), case in drawable
        path = joinpath(directories[case.label], name)
        isfile(path) || continue
        A, value, _ = read_columns(path)
        tables[(index, case.label)] = (A, value)
    end
    isempty(tables) && return nothing
    xlimits = padded_limits(reduce(vcat, [table[1] for table in values(tables)]))

    shown = [
        case for case in drawable if
        any(haskey(tables, (index, case.label)) for index in eachindex(quantities))
    ]
    has_legend = length(shown) > 1
    first_row = has_legend ? 2 : 1

    figure = Figure(size = (1200, 950))
    handles = Dict{String, Any}()
    for (index, (_, ylabel)) in pairs(quantities)
        row, column = fldmod1(index, 2)
        axis = panel(
            figure[first_row + row - 1, column];
            ylabel = ylabel,
            xlabel = L"Fragment mass $A$",
            xlabelvisible = row == 2,
            xticklabelsvisible = row == 2,
            xticks = 80:20:180,
        )
        ordinates = Float64[]
        for case in shown
            haskey(tables, (index, case.label)) || continue
            A, value = tables[(index, case.label)]
            colour, linestyle = SYSTEM_STYLE[case.label]
            curve =
                lines!(axis, A, value; color = colour, linestyle = linestyle, linewidth = 3)
            get!(handles, case.label, curve)
            append!(ordinates, value)
        end
        xlims!(axis, xlimits...)
        isempty(ordinates) || ylims!(axis, padded_limits(ordinates)...)
    end

    if has_legend
        Legend(
            figure[1, 1:2],
            [handles[case.label] for case in shown],
            [case.notation for case in shown];
            tellwidth = false,
            tellheight = true,
            padding = (0, 0, 0, 0),
            colgap = 24,
            patchsize = (50, 20),
        )
    end
    rowgap!(figure.layout, 10)
    colgap!(figure.layout, 12)

    provenance = Dict{String, Any}(
        case.label => system_record(directories[case.label], MeasuredSet[]) for
        case in shown
    )
    export_figure(figure, "residual_quantities", provenance)
    return nothing
end

"""
Draw the isobaric charge distribution each run was weighted by, ΔZ(A) and σ_Z(A) in two columns
with one row per system, against the fitted-parameter and systematics `Zₚ` model and the mean
values. All rows share one heavy-fragment mass range.
"""
function figure_charge_distribution(drawable::Vector{Case})
    isempty(drawable) && return nothing
    directories = [run_directory(case) for case in drawable]
    tables = Dict{Tuple{String, String}, NTuple{3, Vector{Float64}}}()
    for (case, directory) in zip(drawable, directories), source in CHARGE_SOURCES
        path = joinpath(directory, source.file)
        isfile(path) && (tables[(case.label, source.file)] = read_columns(path))
    end
    isempty(tables) && return nothing
    xlimits = padded_limits(reduce(vcat, [table[1] for table in values(tables)]))

    present = [
        source for source in CHARGE_SOURCES if
        any(haskey(tables, (case.label, source.file)) for case in drawable)
    ]
    has_legend = length(present) > 1
    first_row = has_legend ? 2 : 1

    n = length(drawable)
    figure = Figure(size = (1200, 600 + 350 * (n - 1)))
    handles = Dict{String, Any}()
    for (row, case) in pairs(drawable)
        bottom = row == n
        for (column, ylabel) in pairs((L"\Delta Z(A)", L"\sigma_Z(A)"))
            axis = panel(
                figure[first_row + row - 1, column];
                ylabel = ylabel,
                xlabel = L"Heavy fragment mass $A_\mathrm{H}$",
                xlabelvisible = bottom,
                xticklabelsvisible = bottom,
                xticks = 110:10:180,
            )
            ordinates = Float64[]
            for source in CHARGE_SOURCES
                haskey(tables, (case.label, source.file)) || continue
                A, ΔZ, σ_Z = tables[(case.label, source.file)]
                value = column == 1 ? ΔZ : σ_Z
                curve = lines!(
                    axis,
                    A,
                    value;
                    color = source.colour,
                    linestyle = source.linestyle,
                    linewidth = 3,
                )
                get!(handles, source.file, curve)
                append!(ordinates, value)
            end
            xlims!(axis, xlimits...)
            isempty(ordinates) || ylims!(axis, padded_limits(ordinates)...)
            # ΔZ falls with mass, leaving the upper right clear.
            column == 1 && text!(
                axis,
                0.97,
                0.95;
                text = case.notation,
                space = :relative,
                align = (:right, :top),
                fontsize = ANNOTATION_SIZE,
            )
        end
    end

    if has_legend
        Legend(
            figure[1, 1:2],
            [handles[source.file] for source in present],
            [source.label for source in present];
            nbanks = 2,
            tellwidth = false,
            tellheight = true,
            padding = (0, 0, 0, 0),
            colgap = 24,
            patchsize = (50, 20),
        )
    end
    rowgap!(figure.layout, 10)
    colgap!(figure.layout, 12)

    provenance = Dict{String, Any}(
        case.label => system_record(directory, MeasuredSet[]) for
        (case, directory) in zip(drawable, directories)
    )
    export_figure(figure, "charge_distribution", provenance)
    return nothing
end

# Prompt-γ datasets per fission left out of the summary, by accession, with the reason. The
# members of a dependent group are left out where the group's own average is drawn.
const GAMMA_PER_FISSION_SKIP = Dict(
    "23541005" => "the γ energy between 1 and 6 MeV only, not a total",
    "231970023" => "one detector of Oberstedt 2015, whose average 231970093 is drawn",
    "231970033" => "one detector of Oberstedt 2015, whose average 231970093 is drawn",
    "231970053" => "one detector of Oberstedt 2015, whose average 231970093 is drawn",
    "231970063" => "one detector of Oberstedt 2015, whose average 231970093 is drawn",
    "231970073" => "one detector of Oberstedt 2015, whose average 231970093 is drawn",
    "231970083" => "one detector of Oberstedt 2015, whose average 231970093 is drawn",
    "31778002" => "one detector of Gatera 2017, whose average 31778005 is drawn",
    "31778003" => "one detector of Gatera 2017, whose average 31778005 is drawn",
    "31778004" => "one detector of Gatera 2017, whose average 31778005 is drawn",
)

"""
A measured prompt-γ energy per fission: the system, the author–year label, the file, the value,
its uncertainty (zero where none is stated), and whether the record flags it — coded partial in
γ energy (`PAR`) or marked preliminary — which the figure draws as an open marker.
"""
struct GammaValue
    system::String
    label::String
    path::String
    value::Float64
    σ::Float64
    flagged::Bool
end

"""
    gamma_values(case) -> (values, excluded)

The total prompt-γ energies per fission retrieved for one system, each file one value, less those
in `GAMMA_PER_FISSION_SKIP`, which are returned with their reasons.
"""
function gamma_values(case::Case)
    directory = joinpath(MEASURED, case.label, "E_gamma_tot")
    values = GammaValue[]
    excluded = Pair{String, String}[]
    isdir(directory) || return values, excluded
    record = retrieval_record(case.label, "E_gamma_tot")
    path = joinpath(directory, "retrieval.toml")
    warnings = isfile(path) ? get(TOML.parsefile(path), "datasets", Dict()) : Dict()
    preliminary = string(get(warnings, "preliminary_warning", ""))
    for name in sort(filter(endswith(".dat"), readdir(directory)))
        tag = accession(name)
        if haskey(GAMMA_PER_FISSION_SKIP, tag)
            push!(excluded, tag => GAMMA_PER_FISSION_SKIP[tag])
            continue
        end
        rows = read_rows(joinpath(directory, name))
        length(rows) == 1 ||
            throw(ArgumentError("$name holds $(length(rows)) values, not one"))
        row = only(rows)
        σ = length(row) >= 2 && isfinite(row[2]) ? row[2] : 0.0
        entry = get(record, tag, Dict{String, Any}())
        flagged =
            get(entry, "gamma_range_partial", false) === true || occursin(tag, preliminary)
        push!(
            values,
            GammaValue(
                case.label,
                author_year(name),
                joinpath(directory, name),
                row[1],
                σ,
                flagged,
            ),
        )
    end
    return values, excluded
end

"Weighted least-squares slope of `y` on `x`, weights `1/σ²` where every `σ` is positive."
function fitted_slope(x::Vector{Float64}, y::Vector{Float64}, σ::Vector{Float64})
    w = all(>(0), σ) ? 1 ./ σ .^ 2 : ones(length(x))
    W = sum(w)
    x̄ = sum(w .* x) / W
    ȳ = sum(w .* y) / W
    return sum(w .* (x .- x̄) .* (y .- ȳ)) / sum(w .* (x .- x̄) .^ 2)
end

"Linear interpolation of a tabulated curve at `x`, `NaN` outside it."
function interpolate_curve(xs::Vector{Float64}, ys::Vector{Float64}, x::Real)
    (x < first(xs) || x > last(xs)) && return NaN
    i = searchsortedlast(xs, x)
    i == length(xs) && return ys[end]
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return (1 - t) * ys[i] + t * ys[i + 1]
end

"""
Draw the total prompt-γ energy per fission of every system with a run against the measured
values, one row per dataset grouped by system, the model as a black bar across its system's
rows. What the model draws is the energy its neutrons leave, which counts conversion electrons,
X-rays and isomers as well; the measured values are γ rays above a threshold within a time
window, which their records state where the entry or its publication does.
"""
function figure_gamma_energy_per_fission(drawable::Vector{Case})
    groups =
        Tuple{Case, String, Float64, Vector{GammaValue}, Vector{Pair{String, String}}}[]
    for case in drawable
        directory = run_directory(case)
        modelled = get(read_totals(directory), "E_gamma_per_fission", NaN)
        isfinite(modelled) || continue
        values, excluded = gamma_values(case)
        isempty(values) && continue
        push!(groups, (case, directory, modelled, values, excluded))
        @printf("  %-10s E_γ = %.2f MeV per fission\n", case.label, modelled)
    end
    isempty(groups) && return nothing

    rows = sum(length(g[4]) for g in groups)
    figure = Figure(size = (1000, 140 + 46 * (rows + length(groups))))
    axis = panel(
        figure[2, 1];
        xlabel = L"Prompt $\gamma$ energy per fission $E_{\gamma,\mathrm{tot}}$ [MeV]",
        ygridvisible = false,
        xticks = 0:1:20,
    )
    positions = Float64[]
    labels = String[]
    provenance = Dict{String, Any}()
    position = 0.0
    measured_handle = nothing
    flagged_handle = nothing
    model_handle = nothing
    everything = Float64[]
    system_labels = Tuple{MathText, Float64}[]
    for (index, (case, directory, modelled, values, excluded)) in pairs(groups)
        index > 1 && (position += 1.0)
        # Colours are distinct within a system's rows and a label keeps its colour across them.
        begin_panel!()
        first_row = position + 1
        for value in values
            position += 1
            push!(positions, position)
            push!(labels, value.label)
            colour, marker = dataset_style(value.label)
            value.σ > 0 && rangebars!(
                axis,
                [position],
                [value.value - value.σ],
                [value.value + value.σ];
                direction = :x,
                color = (colour, 0.5),
                linewidth = 1.5,
            )
            handle = scatter!(
                axis,
                [value.value],
                [position];
                color = value.flagged ? :transparent : colour,
                marker = marker,
                markersize = 16,
                strokecolor = darker(colour),
                strokewidth = value.flagged ? 2 : 1.5,
            )
            value.flagged ? (flagged_handle = handle) : (measured_handle = handle)
            push!(everything, value.value - value.σ, value.value + value.σ)
        end
        model_handle = linesegments!(
            axis,
            [modelled, modelled],
            [first_row - 0.45, position + 0.45];
            color = :black,
            linewidth = 3,
        )
        push!(everything, modelled)
        push!(system_labels, (case.notation, (first_row + position) / 2))
        provenance[case.label] = Dict{String, Any}(
            "run" => basename(directory),
            "gitcommit" => gitcommit(directory),
            "measured" => [relpath(v.path, MEASURED) for v in values],
            "excluded" => ["$tag: $reason" for (tag, reason) in excluded],
        )
        if index < length(groups)
            hlines!(
                axis,
                [position + 1.0];
                color = GUIDE_COLOUR,
                linestyle = :dash,
                linewidth = 1.0,
            )
        end
    end
    low, high = padded_limits(everything)
    # A band on the right holds the system labels, clear of every marker and bar.
    right = high + 0.3 * (high - low)
    xlims!(axis, low, right)
    # Descending limits put the first system at the top.
    ylims!(axis, last(positions) + 0.8, first(positions) - 0.8)
    axis.yticks = (positions, labels)
    for (notation, middle) in system_labels
        text!(
            axis,
            right - 0.015 * (right - low),
            middle;
            text = notation,
            align = (:right, :center),
            fontsize = ANNOTATION_SIZE,
        )
    end

    handles = Any[model_handle]
    entries = Any["This model"]
    if measured_handle !== nothing
        push!(
            handles,
            MarkerElement(;
                marker = :circle,
                color = GUIDE_LABEL_COLOUR,
                strokecolor = darker(GUIDE_LABEL_COLOUR),
                strokewidth = 1.5,
                markersize = 14,
            ),
        )
        push!(entries, "Measured")
    end
    if flagged_handle !== nothing
        push!(handles, marked_key())
        push!(entries, "Partial γ range or preliminary")
    end
    Legend(
        figure[1, 1],
        handles,
        entries;
        orientation = :horizontal,
        tellwidth = false,
        tellheight = true,
        padding = (0, 0, 0, 0),
        colgap = 20,
        patchsize = (35, 20),
    )
    rowgap!(figure.layout, 8)
    export_figure(figure, "prompt_gamma_energy_per_fission", provenance)
    return nothing
end

"""
Draw the prompt-γ energy of ²⁵²Cf(sf) against measurement, three panels: the total per fission
against TKE with both fitted slopes, the energy of one fragment against its mass, and the total
per fission against the mass of either fragment. Only ²⁵²Cf has these measured.
"""
function figure_gamma_energy_cf252(drawable::Vector{Case})
    index = findfirst(case -> case.label == "Cf252_sf", drawable)
    index === nothing && return nothing
    case = drawable[index]
    directory = run_directory(case)
    specification = (
        (
            "E_gamma_tot_vs_TKE.dat",
            "E_gamma_tot_vs_TKE",
            "Total kinetic energy TKE [MeV]",
            L"E_{\gamma,\mathrm{tot}}\ [\mathrm{MeV}]",
        ),
        (
            "E_gamma_vs_A.dat",
            "E_gamma_vs_A",
            L"Fragment mass $A$",
            L"E_{\gamma}\ [\mathrm{MeV}]",
        ),
        (
            "E_gamma_tot_vs_A.dat",
            "E_gamma_tot_vs_A",
            L"Fragment mass $A$",
            L"E_{\gamma,\mathrm{tot}}\ [\mathrm{MeV}]",
        ),
    )
    available = [isfile(joinpath(directory, table)) for (table, _, _, _) in specification]
    all(available) ||
        (@warn "no prompt-γ tables to draw; re-run the system" directory; return nothing)

    figure = Figure(size = (1500, 600))
    legend_handles = Dict{String, Any}()
    measured_sets = MeasuredSet[]
    model_handle = nothing
    for (column, (table, observable, xlabel, ylabel)) in pairs(specification)
        x, y, _ = read_columns(joinpath(directory, table))
        sets = load_measured(case.label, observable, (-Inf, Inf))
        span = isempty(sets) ? extrema(x) : extrema(reduce(vcat, (set.x for set in sets)))
        xlimits = padded_limits(collect(span))
        inside = first(xlimits) .<= x .<= last(xlimits)
        axis = panel(figure[2, column]; xlabel = xlabel, ylabel = ylabel)
        begin_panel!()
        if observable != "E_gamma_tot_vs_TKE"
            symmetric = case.compound_mass / 2
            first(xlimits) < symmetric < last(xlimits) && vlines!(
                axis,
                [symmetric];
                color = GUIDE_COLOUR,
                linestyle = :dash,
                linewidth = 1.5,
            )
        end
        for set in sets
            legend_handles[set.label] = draw_measured!(axis, set)
            push!(measured_sets, set)
        end
        model_handle = lines!(axis, x[inside], y[inside]; color = :black, linewidth = 3)
        values = vcat(y[inside], (set.y for set in sets)...)
        low, high = padded_limits(values)
        lines = MathText[case.notation]
        if observable == "E_gamma_tot_vs_TKE" && !isempty(sets)
            set = only(sets)
            measured_slope = fitted_slope(set.x, set.y, set.σ)
            modelled = [interpolate_curve(x, y, e) for e in set.x]
            model_slope = fitted_slope(set.x, modelled, zeros(length(set.x)))
            @printf(
                "  Cf252_sf   dE_γ,tot/dTKE: model %.3f, %s %.3f MeV/MeV\n",
                model_slope,
                set.label,
                measured_slope
            )
            push!(
                lines,
                MathText(@sprintf("Slope, measured: \$%.3f\$", measured_slope)),
                MathText(@sprintf("Slope, model: \$%.3f\$", model_slope)),
            )
        end
        high += length(lines) * ANNOTATION_BAND * (high - low)
        xlims!(axis, xlimits...)
        ylims!(axis, low, high)
        # Upper right against TKE, where the energy falls; upper left against mass, where the
        # light fragments far from symmetry leave least and the A₀/2 guide is far.
        corner = observable == "E_gamma_tot_vs_TKE" ? :upper_right : :upper_left
        annotate!(axis, lines; corner = corner)
    end
    labels = sort(collect(keys(legend_handles)))
    Legend(
        figure[1, 1:3],
        vcat([legend_handles[l] for l in labels], [model_handle]),
        vcat(labels, ["This model"]);
        orientation = :horizontal,
        tellwidth = false,
        tellheight = true,
        padding = (0, 0, 0, 0),
        colgap = 20,
        patchsize = (35, 20),
    )
    rowgap!(figure.layout, 8)
    colgap!(figure.layout, 12)
    export_figure(
        figure,
        "prompt_gamma_energy_cf252",
        Dict{String, Any}(case.label => system_record(directory, measured_sets)),
    )
    return nothing
end

function main()
    mkpath(FIGURES)
    isdir(RUNS) || @warn """
        no runs directory; solve the systems first with \
        `julia scripts/run.jl config/<system>.toml`.""" directory = RUNS
    isdir(MEASURED) || @warn """
        no data directory; the figures will carry the model curves alone. Retrieve the \
        measurements with ExforFissionData.jl and lay them out as described at the top of this \
        file.""" directory = MEASURED
    println("drawing from:")
    drawable = filter(available, CASES)
    figure_multiplicity(drawable)
    figure_spectra(drawable)
    figure_mean_neutron_energy(drawable)
    figure_product_yield(drawable)
    figure_charge_yield(drawable)
    foreach(figure_independent_yield, drawable)
    figure_neutron_number_yield(drawable)
    figure_product_kinetic_energy(drawable)
    figure_product_kinetic_energy_by_charge(drawable)
    figure_product_total_kinetic_energy(drawable)
    figure_residual_quantities(drawable)
    figure_charge_distribution(drawable)
    figure_gamma_energy_per_fission(drawable)
    figure_gamma_energy_cf252(drawable)
    return nothing
end

main()
