# Validation figures: the model against measurement.
#
#   julia plotting/figures.jl [runs directory] [data directory]
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
#
# which is the layout that package writes, with a retrieval.toml beside the files that is read
# for the orientation of a ratio and for datasets a later one supersedes. Columns are taken by
# position, not by header text. Both directories default to data/ under the project and neither
# is version-controlled.
#
# Each figure is written to results/ as PDF and PNG, beside a TOML sidecar that names, for every
# system drawn, the run directory, the commit recorded in its metadata, and the measured files.

include(joinpath(@__DIR__, "activate.jl"))
include(joinpath(@__DIR__, "theme.jl"))

using Printf
using TOML

set_theme!(PUBLICATION_THEME)

const ROOT = dirname(@__DIR__)
const FIGURES = joinpath(ROOT, "results")
const RUNS = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "data", "sims")
const MEASURED = length(ARGS) >= 2 ? ARGS[2] : joinpath(ROOT, "data")

# LaTeXStrings reaches this environment only through Makie; its string type is taken from the
# macro.
const MathText = typeof(L"")

"""
One fissioning system: the token that names its directories, the typeset notation it is drawn
under, and the evaluated total multiplicity it is compared against — ENDF/B-VIII.0 for the
induced systems and the Axton standard for 252-Cf.

`label` is the system token, as `system_label` returns it; `notation` is the typeset form, as
`system_notation` returns it. They are written out here rather than called, because this
environment draws from tables and does not load the package.
"""
struct Case
    label::String
    notation::MathText
    compound_mass::Int
    evaluated::Float64
end

const CASES = [
    Case("Cf252_sf", L"$^{252}$Cf(sf)", 252, 3.7573),
    Case("U235_nth", L"$^{235}$U(nth,f)", 236, 2.4140),
    Case("Pu239_nth", L"$^{239}$Pu(nth,f)", 240, 2.8750),
    Case("U233_nth", L"$^{233}$U(nth,f)", 234, 2.4968),
]

# The system every spectrum ratio is taken against.
const REFERENCE = first(CASES)

# Measured masses drawn. Outside this window the yields vanish and the measured uncertainties
# are unbounded.
const MASS_WINDOW = (76.0, 176.0)

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
    author = replace(
        author,
        "Goeoek" => "Göök",
        "Budtz-jorgensen" => "Budtz-Jørgensen",
        "Al-adili" => "Al-Adili",
    )
    return strip(last(split(author, ". "))) * " " * parts[3]
end

# The measured datasets each figure shows, by EXFOR accession: a fixed selection per system, for
# legibility. The data directory holds the whole retrieval. A system or observable not listed
# shows every dataset retrieved.
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
reach into the margins and the annotation band.
"""
function ordinate_limits(
    model::AbstractVector{<:Real},
    sets::Vector{MeasuredSet};
    cleared::Union{Nothing, Tuple{Real, Real}} = nothing,
    bars::Symbol = :upper,
)
    values = collect(Float64, model)
    sign =
        bars === :upper ? 1 :
        bars === :lower ? -1 :
        throw(ArgumentError("bars must be :upper or :lower, got :$bars"))
    for set in sets
        append!(values, set.y)
        cleared === nothing && continue
        under = first(cleared) .<= set.x .<= last(cleared)
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

"The provenance entry of one system: the run drawn, its commit, and the measured files used."
function system_record(directory::AbstractString, sets::Vector{MeasuredSet})
    measured = String[]
    for set in sets, path in set.paths
        push!(measured, relpath(path, MEASURED))
    end
    return Dict{String, Any}(
        "run" => basename(directory),
        "gitcommit" => gitcommit(directory),
        "measured" => measured,
    )
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
Stacked lines of in-axis annotation in one corner, `ANNOTATION_LINE` apart, read top to bottom.
"""
function annotate!(axis::Axis, lines::Vector{MathText}; corner::Symbol)
    # Lines read top to bottom in either corner: the first line is the outermost one in the
    # upper corner and the innermost one in the lower corner.
    if corner === :upper_left
        position, align = (0.03, 0.96), (:left, :top)
        step = index -> -(index - 1) * ANNOTATION_LINE
    elseif corner === :lower_left
        position, align = (0.03, 0.04), (:left, :bottom)
        step = index -> (length(lines) - index) * ANNOTATION_LINE
    else
        throw(ArgumentError("corner must be :upper_left or :lower_left, got :$corner"))
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

"One panel of a per-system grid: the system, its run, the model curve and the measured sets."
struct SystemPanel
    case::Case
    directory::String
    x::Vector{Float64}
    y::Vector{Float64}
    sets::Vector{MeasuredSet}
end

"""
The panels of one quantity against fragment mass: for every system with a run that wrote
`table`, the model curve and the measured sets of `observable` inside the run's own mass range
— measured points are compared only where the run solved the model — within the outer window.
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
        window = (max(first(MASS_WINDOW), minimum(A)), min(last(MASS_WINDOW), maximum(A)))
        sets = load_measured(case.label, observable, window)
        push!(panels, SystemPanel(case, directory, A, value, sets))
    end
    return panels
end

"""
Draw a quantity against fragment mass for every panel, on the shared grid: the measured sets as
markers, the model as the black line, the A₀/2 guide labelled above the frame, each panel's
legend above it, and `annotate(panel)` as lines in the upper left, above the data. All panels
share one mass range; ordinate ranges are each panel's own.
"""
function figure_mass_grid(
    name::AbstractString,
    panels::Vector{SystemPanel},
    ylabel,
    annotate::Function,
)
    isempty(panels) && return nothing
    masses = Float64[]
    for p in panels
        append!(masses, p.x)
        foreach(set -> append!(masses, set.x), p.sets)
        push!(masses, p.case.compound_mass / 2)
    end
    xlimits = padded_limits(masses)
    annotated = (
        first(xlimits),
        first(xlimits) + ANNOTATION_EXTENT * (last(xlimits) - first(xlimits)),
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
            ylabel = ylabel,
            xlabel = L"Fragment mass $A$",
            ylabelvisible = column == 1,
            xlabelvisible = bottom,
            xticklabelsvisible = bottom,
            xticks = 80:20:180,
        )

        symmetric = p.case.compound_mass / 2
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
        low, high = ordinate_limits(p.y, p.sets; cleared = annotated)
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
    figure_mass_grid("multiplicity_vs_experiment", panels, L"\nu(A)", annotate)
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
        p -> MathText[p.case.notation],
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
    figure_residual_quantities(drawable)
    figure_charge_distribution(drawable)
    return nothing
end

main()
