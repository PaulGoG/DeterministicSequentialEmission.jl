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
#   <data>/<system>/spectrum_maxwellian_ratio_vs_E/*.dat     E spectrum_maxwellian_ratio [...]
#
# which is the layout that package writes. Columns are taken by position, not by header text.
# Both directories default to data/ under the project and neither is version-controlled.
#
# Each figure is written to plots/ as PDF and PNG, beside a TOML sidecar that names, for every
# system drawn, the run directory, the commit recorded in its metadata, and the measured files.

include(joinpath(@__DIR__, "activate.jl"))
include(joinpath(@__DIR__, "theme.jl"))

using Printf
using TOML

set_theme!(PUBLICATION_THEME)

const ROOT = dirname(@__DIR__)
const FIGURES = joinpath(ROOT, "plots")
const RESULTS = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "data", "sims")
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

# Measured masses drawn. Outside this window the yields vanish and the measured uncertainties
# are unbounded.
const MASS_WINDOW = (76.0, 176.0)

# Energies drawn in the spectrum. Below the window the measured ratio is dominated by detector
# thresholds, above it by counting statistics in the tail.
const ENERGY_WINDOW = (0.08, 11.0)

# Measured mean laboratory energy of the 252-Cf prompt fission neutron spectrum [MeV].
const MEASURED_MEAN_ENERGY = 2.13

# Maxwellian temperature taken when a run's configuration records none [MeV].
const FALLBACK_MAXWELLIAN_TEMPERATURE = 1.42

# Margin on either side of the data range, as a fraction of that range, and the band reserved
# above the data for two lines of in-axis annotation.
const MARGIN = 0.05
const ANNOTATION_BAND = 0.30

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
    system = joinpath(RESULTS, case.label)
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
A measured dataset: its author–year label, the file it was read from, and the points inside the
drawn window with their uncertainties, zero where none is reported.
"""
struct MeasuredSet
    label::String
    path::String
    x::Vector{Float64}
    y::Vector{Float64}
    σ::Vector{Float64}
end

# Files are named <accession>_<author>_<year>.dat by the retrieval.
function author_year(name::AbstractString)
    parts = split(replace(name, ".dat" => ""), "_")
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

"The measured datasets of one observable of one system, restricted to `window`."
function load_measured(system::AbstractString, observable::AbstractString, window)
    directory = joinpath(MEASURED, system, observable)
    sets = MeasuredSet[]
    isdir(directory) || return sets
    for name in sort(filter(endswith(".dat"), readdir(directory)))
        path = joinpath(directory, name)
        x, y, σ = read_columns(path)
        keep = first(window) .<= x .<= last(window)
        any(keep) || continue
        push!(sets, MeasuredSet(author_year(name), path, x[keep], y[keep], σ[keep]))
    end
    return sets
end

function available(case::Case)
    directory = run_directory(case)
    if directory !== nothing && isfile(joinpath(directory, "nu_vs_A.dat"))
        println("  ", rpad(case.label, 10), " ", basename(directory))
        return true
    end
    @warn """
        no run to draw; solve it first with `julia scripts/run.jl \
        config/$(case.label).toml`.""" expected = joinpath(RESULTS, case.label)
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

"Ordinate limits spanning a model curve and every measured value; error bars may reach into the
margins and the annotation band above."
function ordinate_limits(model::AbstractVector{<:Real}, sets::Vector{MeasuredSet})
    values = collect(Float64, model)
    for set in sets
        append!(values, set.y)
    end
    return padded_limits(values)
end

panel(slot; kwargs...) = Axis(slot; xlabelpadding = 8, ylabelpadding = 8, kwargs...)

# Legends hold at most three entries per row.
legend_rows(entries::Integer) = cld(entries, 3)

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
    return Dict{String, Any}(
        "run" => basename(directory),
        "gitcommit" => gitcommit(directory),
        "measured" => String[relpath(set.path, MEASURED) for set in sets],
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
Draw ν(A) for every system with a run, one stacked panel each, with the measured datasets as
markers and the model as the black line. Each panel carries the A₀/2 guide, the system notation,
and ⟨ν⟩ against the evaluated value; all panels share one mass range.
"""
function figure_multiplicity(drawable::Vector{Case})
    isempty(drawable) && return nothing
    panels = map(drawable) do case
        directory = run_directory(case)
        A, ν, _ = read_columns(joinpath(directory, "nu_vs_A.dat"))
        # Measured points are compared only where the run solved the model: the run's own mass
        # range, within the outer window.
        window =
            (max(first(MASS_WINDOW), minimum(A)), min(last(MASS_WINDOW), maximum(A)))
        sets = load_measured(case.label, "nu_vs_A", window)
        (; case, directory, A, ν, sets)
    end

    masses = Float64[]
    for p in panels
        append!(masses, p.A)
        foreach(set -> append!(masses, set.x), p.sets)
        push!(masses, p.case.compound_mass / 2)
    end
    xlimits = padded_limits(masses)

    n = length(panels)
    figure = Figure(size = (900, 600 + 350 * (n - 1)))
    provenance = Dict{String, Any}()
    for (row, p) in pairs(panels)
        bottom = row == n
        cell = GridLayout(figure[row, 1])
        has_legend = !isempty(p.sets)
        slot = cell[has_legend ? 2 : 1, 1]
        axis = panel(
            slot;
            ylabel = L"\nu(A)",
            xlabel = L"Fragment mass $A$",
            xlabelvisible = bottom,
            xticklabelsvisible = bottom,
            xticks = 80:10:180,
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
        push!(handles, lines!(axis, p.A, p.ν; color = :black, linewidth = 3))
        push!(labels, "This model")

        low, high = ordinate_limits(p.ν, p.sets)
        ylimits = (low, high + ANNOTATION_BAND * (high - low))
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
        label_guide!(slot, symmetric, L"A_0/2", xlimits; side = :top)

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
        takeaway = MathText(
            @sprintf(
                "\$\\langle \\nu \\rangle = %.2f\$, evaluated %.2f (\$%+.1f\$ %%)",
                modelled,
                evaluated,
                deviation
            )
        )
        # The upper left is clear of the sawtooth: light fragments far from symmetry emit least.
        text!(
            axis,
            0.03,
            0.96;
            text = p.case.notation,
            space = :relative,
            align = (:left, :top),
            fontsize = ANNOTATION_SIZE,
        )
        text!(
            axis,
            0.03,
            0.96;
            text = takeaway,
            space = :relative,
            align = (:left, :top),
            offset = (0, -ANNOTATION_LINE),
            fontsize = ANNOTATION_SIZE,
        )

        if has_legend
            Legend(
                cell[1, 1],
                handles,
                labels;
                nbanks = legend_rows(length(handles)),
                tellwidth = false,
                tellheight = true,
                padding = (0, 0, 0, 0),
                colgap = 16,
                patchsize = (35, 20),
            )
            rowgap!(cell, 8)
        end
        provenance[p.case.label] = system_record(p.directory, p.sets)
    end
    n > 1 && rowgap!(figure.layout, 12)

    export_figure(figure, "multiplicity_vs_experiment", provenance)
    return nothing
end

maxwellian(E, T) = 2 / sqrt(π) * T^(-3 / 2) * sqrt(E) * exp(-E / T)

"""
Draw the 252-Cf prompt neutron spectrum as its ratio to a Maxwellian at the run's `T_M`, the
model as the black line and the measured ratios as markers. The unity guide is labelled with
`T_M`, and ⟨E⟩ is compared to the measured value in the lower left.
"""
function figure_spectrum(case::Case)
    directory = run_directory(case)
    path = directory === nothing ? "" : joinpath(directory, "spectrum_vs_E.dat")
    if !isfile(path)
        @warn "no spectrum to draw" expected = path
        return nothing
    end
    E, N, _ = read_columns(path)
    T_M = maxwellian_temperature(directory)
    mean_E = get(read_totals(directory), "mean_E", NaN)
    @printf("  ⟨E⟩ = %.4f MeV against a Maxwellian at T_M = %.2f MeV\n", mean_E, T_M)

    low, high = ENERGY_WINDOW
    inside = low .<= E .<= high
    energies = E[inside]
    ratio = N[inside] ./ maxwellian.(energies, T_M)
    sets = load_measured(case.label, "spectrum_maxwellian_ratio_vs_E", ENERGY_WINDOW)

    xlimits = padded_log_limits(vcat(energies, (set.x for set in sets)...))
    ylimits = ordinate_limits(vcat(ratio, 1.0), sets)

    figure = Figure(size = (900, 600))
    has_legend = !isempty(sets)
    slot = figure[has_legend ? 2 : 1, 1]
    axis = panel(
        slot;
        xscale = log10,
        xticks = SPECTRUM_TICKS,
        xlabel = L"Laboratory neutron energy $E$ [MeV]",
        ylabel = L"N(E)\,/\,M(E,\,T_\mathrm{M})",
    )

    hlines!(axis, [1.0]; color = GUIDE_COLOUR, linestyle = :dash, linewidth = 1.5)
    begin_panel!()
    handles = Any[]
    labels = Any[]
    for set in sets
        push!(handles, draw_measured!(axis, set))
        push!(labels, set.label)
    end
    push!(handles, lines!(axis, energies, ratio; color = :black, linewidth = 3))
    push!(labels, "This model")

    xlims!(axis, xlimits...)
    ylims!(axis, ylimits...)
    guide = MathText(@sprintf("\$T_\\mathrm{M} = %.2f\$ MeV", T_M))
    label_guide!(slot, 1.0, guide, ylimits; side = :right)

    # Below the data at low energy: every measured ratio there stays near unity.
    takeaway = MathText(
        @sprintf(
            "\$\\langle E \\rangle = %.3f\$ MeV, measured %.2f MeV",
            mean_E,
            MEASURED_MEAN_ENERGY
        )
    )
    text!(
        axis,
        0.03,
        0.04;
        text = takeaway,
        space = :relative,
        align = (:left, :bottom),
        fontsize = ANNOTATION_SIZE,
    )
    text!(
        axis,
        0.03,
        0.04;
        text = case.notation,
        space = :relative,
        align = (:left, :bottom),
        offset = (0, ANNOTATION_LINE),
        fontsize = ANNOTATION_SIZE,
    )

    if has_legend
        Legend(
            figure[1, 1],
            handles,
            labels;
            nbanks = legend_rows(length(handles)),
            tellwidth = false,
            tellheight = true,
            padding = (0, 0, 0, 0),
            colgap = 16,
            patchsize = (35, 20),
        )
        rowgap!(figure.layout, 8)
    end

    provenance = Dict{String, Any}(case.label => system_record(directory, sets))
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
    isdir(RESULTS) || @warn """
        no runs directory; solve the systems first with \
        `julia scripts/run.jl config/<system>.toml`.""" directory = RESULTS
    isdir(MEASURED) || @warn """
        no data directory; the figures will carry the model curves alone. Retrieve the \
        measurements with ExforFissionData.jl and lay them out as described at the top of this \
        file.""" directory = MEASURED
    println("drawing from:")
    drawable = filter(available, CASES)
    figure_multiplicity(drawable)
    figure_spectrum(first(CASES))
    figure_residual_quantities(drawable)
    figure_charge_distribution(drawable)
    return nothing
end

main()
