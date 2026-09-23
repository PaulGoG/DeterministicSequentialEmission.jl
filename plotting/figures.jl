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
# Both directories default to beside the project and neither is version-controlled; nothing here
# knows where yours lives.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie, MathTeXEngine
using Printf

const ROOT = dirname(@__DIR__)
const FIGURES = joinpath(ROOT, "plots")
const RESULTS = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "data", "sims")
const MEASURED = length(ARGS) >= 2 ? ARGS[2] : joinpath(ROOT, "data")

"""
One fissioning system: the token that names its directories, the typeset notation it is drawn
under, and the evaluated total multiplicity it is compared against — ENDF/B-VIII.0 for the
induced systems and the Axton standard for 252-Cf.

`label` is the system token, as `system_label` returns it; `notation` is the typeset form, as
`system_notation` returns it. One name may not mean both. They are written out here rather than
called, because this environment draws from tables and does not load the package.
"""
struct Case
    label::String
    notation::Any
    compound_mass::Int
    measured::Float64
end

const CASES = [
    Case("Cf252_sf", L"$^{252}$Cf(sf)", 252, 3.7573),
    Case("U235_nth", L"$^{235}$U(nth,f)", 236, 2.4140),
    Case("Pu239_nth", L"$^{239}$Pu(nth,f)", 240, 2.8750),
    Case("U233_nth", L"$^{233}$U(nth,f)", 234, 2.4968),
]

# One common abscissa, so the panels can be read against each other and the inner tick labels
# dropped. Nothing is plotted outside it.
const MASS_RANGE = (76.0, 176.0)

# The heavy branch alone, which is where a charge distribution is tabulated. Wide enough for the
# widest window any shipped system sweeps.
const CHARGE_MASS_RANGE = (115.0, 176.0)

const THEME = Theme(
    fonts = (; regular = texfont(:text), bold = texfont(:bold), italic = texfont(:italic)),
    fontsize = 9,
    figure_padding = 6,
    Axis = (
        xgridstyle = :dash,
        ygridstyle = :dash,
        xgridcolor = (:grey, 0.12),
        ygridcolor = (:grey, 0.12),
        xminorticksvisible = false,
        yminorticksvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticksize = 3,
        yticksize = 3,
    ),
)

const PALETTE = [
    RGBf(0.835, 0.369, 0.0),
    RGBf(0.0, 0.620, 0.451),
    RGBf(0.800, 0.475, 0.655),
    RGBf(0.337, 0.706, 0.914),
    RGBf(0.941, 0.667, 0.0),
    RGBf(0.494, 0.184, 0.557),
]
const MARKERS = [:circle, :rect, :utriangle, :diamond, :star5, :cross]
# Colour carries the series family and style carries the role, so the figures survive greyscale.
const LINESTYLES = [:solid, :dash, :dashdot, :dot, (:dot, :dense), :solid]

"Read a two- or three-column file: abscissa, value, optional uncertainty."
function read_columns(path)
    x = Float64[]
    y = Float64[]
    σ = Float64[]
    for (index, line) in enumerate(eachline(path))
        index == 1 && continue
        fields = split(strip(line))
        length(fields) >= 2 || continue
        push!(x, parse(Float64, fields[1]))
        push!(y, parse(Float64, fields[2]))
        push!(σ, length(fields) >= 3 ? something(tryparse(Float64, fields[3]), 0.0) : 0.0)
    end
    return x, y, σ
end

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
function read_totals(case::Case)
    directory = run_directory(case)
    directory === nothing && return Dict{String, Float64}()
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

curve(case::Case, name) = read_columns(joinpath(run_directory(case), name))

function measured_sets(system, observable)
    directory = joinpath(MEASURED, system, observable)
    isdir(directory) || return Tuple{String, String}[]
    files = sort(filter(name -> endswith(name, ".dat"), readdir(directory)))
    return [(_author(name), joinpath(directory, name)) for name in files]
end

# Files are named <accession>_<author>_<year>.dat by the retrieval.
function _author(name)
    parts = split(replace(name, ".dat" => ""), "_")
    length(parts) >= 3 || return name
    author = replace(parts[2], "." => ". ")
    author = replace(author, "Goeoek" => "Göök", "Budtz-jorgensen" => "Budtz-Jørgensen")
    return strip(last(split(author, ". "))) * " " * parts[3]
end

function _available(case::Case)
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

function figure_multiplicity()
    set_theme!(THEME)
    drawable = filter(_available, CASES)
    isempty(drawable) && return nothing

    figure = Figure(size = (510, 30 + 140 * length(drawable)))
    summary = String[]

    for (row, case) in pairs(drawable)
        A, ν, _ = curve(case, "nu_vs_A.dat")
        totals = read_totals(case)
        modelled = get(totals, "nu_per_fission", NaN)
        push!(
            summary,
            @sprintf(
                "%-24s ⟨ν⟩ = %.4f   evaluated %.4f   %+.1f %%",
                string(case.notation),
                modelled,
                case.measured,
                100 * (modelled - case.measured) / case.measured
            )
        )

        axis = Axis(
            figure[2 * row, 1];
            ylabel = L"\nu(A)",
            xlabel = row == length(drawable) ? L"Fragment mass $A$" : "",
            xticks = 80:10:170,
            xticklabelsvisible = row == length(drawable),
        )

        handles = []
        labels = Any[]
        for (index, (author, path)) in pairs(measured_sets(case.label, "nu_vs_A"))
            index > length(PALETTE) && break
            A_measured, ν_measured, σ = read_columns(path)
            keep = first(MASS_RANGE) .<= A_measured .<= last(MASS_RANGE)
            handle = scatter!(
                axis,
                A_measured[keep],
                ν_measured[keep];
                color = (PALETTE[index], 0.7),
                marker = MARKERS[index],
                markersize = 4,
            )
            any(>(0), σ[keep]) && errorbars!(
                axis,
                A_measured[keep],
                ν_measured[keep],
                σ[keep];
                color = (PALETTE[index], 0.3),
                linewidth = 0.5,
            )
            push!(handles, handle)
            push!(labels, author)
        end
        model = lines!(axis, A, ν; color = :black, linewidth = 1.6)
        push!(handles, model)
        push!(labels, "This model")

        symmetric = case.compound_mass / 2
        vlines!(axis, [symmetric]; color = (:grey, 0.55), linestyle = :dot, linewidth = 0.8)
        text!(
            axis,
            symmetric + 1.5,
            3.4;
            text = L"A_0/2",
            color = (:grey, 0.75),
            fontsize = 7,
            align = (:left, :top),
        )
        text!(
            axis,
            last(MASS_RANGE) - 2,
            0.1;
            text = case.notation,
            fontsize = 9,
            align = (:right, :bottom),
        )
        xlims!(axis, MASS_RANGE...)
        ylims!(axis, 0.0, 3.5)

        Legend(
            figure[2 * row - 1, 1],
            handles,
            labels;
            orientation = :horizontal,
            nbanks = 1,
            framevisible = false,
            padding = (0, 0, 0, 0),
            colgap = 8,
            patchsize = (12, 2),
            labelsize = 7.5,
        )
    end

    # Gaps are set once the layout is complete: a legend row tight against its panel, and a
    # wider separation between systems.
    for row in 1:length(drawable)
        rowgap!(figure.layout, 2 * row - 1, 2)
        row < length(drawable) && rowgap!(figure.layout, 2 * row, 8)
    end

    path = joinpath(FIGURES, "multiplicity_vs_experiment")
    save(path * ".pdf", figure)
    save(path * ".png", figure; px_per_unit = 3)
    foreach(println, summary)
    println("wrote ", path, ".pdf")
    return nothing
end

maxwellian(E, T) = 2 / sqrt(π) * T^(-3 / 2) * sqrt(E) * exp(-E / T)

# Where the retrieved sets have useful statistics. Below this the ratio is dominated by detector
# thresholds, above it by counting statistics in the tail.
const LOW = 0.08
const HIGH = 11.0
const T_M = 1.42     # the Maxwellian temperature every retrieved set is reported against

function figure_spectrum()
    set_theme!(THEME)
    case = first(CASES)
    directory = run_directory(case)
    path = directory === nothing ? "" : joinpath(directory, "spectrum_vs_E.dat")
    if !isfile(path)
        @warn "no spectrum to draw" expected = path
        return nothing
    end
    E, N, _ = read_columns(path)
    totals = read_totals(case)
    @printf(
        "⟨E⟩ = %.4f MeV against a Maxwellian at T_M = %.2f MeV\n",
        get(totals, "mean_E", NaN),
        T_M
    )

    inside = LOW .<= E .<= HIGH
    ratio = N[inside] ./ maxwellian.(E[inside], T_M)

    figure = Figure(size = (510, 300))
    axis = Axis(
        figure[2, 1];
        xscale = log10,
        xlabel = L"Laboratory neutron energy $E$ [MeV]",
        ylabel = L"N(E)\,/\,M(E,\,T_\mathrm{M})",
    )

    handles = []
    labels = Any[]
    for (index, (author, file)) in
        pairs(measured_sets(case.label, "spectrum_maxwellian_ratio_vs_E"))
        index > length(PALETTE) && break
        energies, observed, σ = read_columns(file)
        keep = LOW .<= energies .<= HIGH
        any(keep) || continue
        handle = scatter!(
            axis,
            energies[keep],
            observed[keep];
            color = (PALETTE[index], 0.6),
            marker = MARKERS[index],
            markersize = 3.5,
        )
        any(>(0), σ[keep]) && errorbars!(
            axis,
            energies[keep],
            observed[keep],
            σ[keep];
            color = (PALETTE[index], 0.25),
            linewidth = 0.5,
        )
        push!(handles, handle)
        push!(labels, author)
    end

    model = lines!(axis, E[inside], ratio; color = :black, linewidth = 1.8)
    push!(handles, model)
    push!(labels, "This model")
    hlines!(axis, [1.0]; color = (:grey, 0.6), linestyle = :dot, linewidth = 0.8)

    text!(
        axis,
        LOW * 1.04,
        0.635;
        text = case.notation,
        fontsize = 9,
        align = (:left, :bottom),
    )
    text!(
        axis,
        HIGH * 0.97,
        0.635;
        text = L"T_\mathrm{M} = 1.42\ \mathrm{MeV}",
        fontsize = 8,
        color = (:grey, 0.8),
        align = (:right, :bottom),
    )
    ticks = [0.1, 0.2, 0.5, 1, 2, 5, 10]
    axis.xticks = (ticks, [t < 1 ? string(t) : string(Int(t)) for t in ticks])
    xlims!(axis, LOW, HIGH)
    ylims!(axis, 0.62, 1.18)

    Legend(
        figure[1, 1],
        handles,
        labels;
        orientation = :horizontal,
        nbanks = 1,
        framevisible = false,
        padding = (0, 0, 0, 0),
        colgap = 10,
        patchsize = (12, 2),
        labelsize = 7.5,
    )
    rowgap!(figure.layout, 1, 4)

    path = joinpath(FIGURES, "prompt_neutron_spectrum")
    save(path * ".pdf", figure)
    save(path * ".png", figure; px_per_unit = 3)
    println("wrote ", path, ".pdf")
    return nothing
end

"""
The quantities internal to the sequential treatment, against fragment mass, for every system a
run wrote. None has an experimental counterpart — they are intermediate to the emission chain,
not observables — so the comparison that means anything is between systems, and between a system
and the shape its own sawtooth implies.

Diagnostic rather than published: not tracked, and regenerated by whoever wants it.
"""
function figure_residual_quantities()
    set_theme!(THEME)
    drawable = filter(_available, CASES)
    isempty(drawable) && return nothing

    panels = (
        ("E_star_vs_A.dat", L"E^*(A)\ [\mathrm{MeV}]"),
        ("E_r_vs_A.dat", L"\bar{E}_r(A)\ [\mathrm{MeV}]"),
        ("T_vs_A.dat", L"\langle T \rangle(A)\ [\mathrm{MeV}]"),
        ("eps_vs_A.dat", L"\langle \varepsilon \rangle(A)\ [\mathrm{MeV}]"),
    )

    figure = Figure(size = (510, 340))
    handles = []
    labels = Any[]
    for (index, (name, ylabel)) in pairs(panels)
        row, column = fldmod1(index, 2)
        axis = Axis(
            figure[row + 1, column];
            ylabel = ylabel,
            xlabel = row == 2 ? L"Fragment mass $A$" : "",
            xticks = 80:20:160,
            xticklabelsvisible = row == 2,
        )
        for (position, case) in pairs(drawable)
            directory = run_directory(case)
            directory === nothing && continue
            path = joinpath(directory, name)
            isfile(path) || continue
            A, quantity, _ = read_columns(path)
            handle = lines!(
                axis,
                A,
                quantity;
                color = PALETTE[position],
                linestyle = LINESTYLES[position],
                linewidth = 1.3,
            )
            if index == 1
                push!(handles, handle)
                push!(labels, case.notation)
            end
        end
        xlims!(axis, MASS_RANGE...)
    end

    Legend(
        figure[1, 1:2],
        handles,
        labels;
        orientation = :horizontal,
        framevisible = false,
        padding = (0, 0, 0, 0),
        colgap = 14,
        patchsize = (16, 2),
        labelsize = 8,
    )
    rowgap!(figure.layout, 1, 4)
    rowgap!(figure.layout, 2, 6)
    colgap!(figure.layout, 1, 14)

    path = joinpath(FIGURES, "residual_quantities")
    save(path * ".pdf", figure)
    save(path * ".png", figure; px_per_unit = 3)
    println("wrote ", path, ".pdf")
    return nothing
end

"""
The isobaric charge distribution each run was weighted by, against the alternatives it could have
taken: the `Zₚ` model with the parameters fitted to this reaction, the same model with its
parameters estimated from the systematics, and the conventional means that carry no mass
dependence.

Every run writes all of them, so this draws tables rather than recomputing anything. Read it as a
diagnostic of the *input*. The gap between the fitted and the estimated parameters is what a
system without a fit of its own inherits; the gap between either and the used table is what the
model still misses; and the flat means are the floor all three improve on.
"""
function figure_charge_distribution()
    set_theme!(THEME)
    drawable = filter(_available, CASES)
    isempty(drawable) && return nothing

    sources = (
        ("charge_distribution_vs_A.dat", "Used by the run", :solid, :black),
        (
            "charge_distribution_vs_A_wahl_per_reaction.dat",
            L"$Z_p$ model, fitted parameters",
            :dashdot,
            PALETTE[2],
        ),
        (
            "charge_distribution_vs_A_wahl_systematics.dat",
            L"$Z_p$ model, systematics",
            :dash,
            PALETTE[1],
        ),
        ("charge_distribution_vs_A_mean.dat", "Mean values", :dot, PALETTE[4]),
    )

    figure = Figure(size = (510, 40 + 105 * length(drawable)))
    # Keyed by file name rather than taken from the first panel, so a source that only some
    # systems have still reaches the legend.
    entries = Dict{String, Any}()
    for (row, case) in pairs(drawable)
        directory = run_directory(case)
        directory === nothing && continue
        bottom = row == length(drawable)
        for (column, (ylabel, index)) in pairs(((L"\Delta Z(A)", 2), (L"\sigma_Z(A)", 3)))
            axis = Axis(
                figure[row + 1, column];
                ylabel = ylabel,
                xlabel = bottom ? L"Heavy fragment mass $A_H$" : "",
                xticklabelsvisible = bottom,
            )
            for (name, label, style, colour) in sources
                path = joinpath(directory, name)
                isfile(path) || continue
                A, second, third = read_columns(path)
                handle = lines!(
                    axis,
                    A,
                    index == 2 ? second : third;
                    color = colour,
                    linestyle = style,
                    linewidth = index == 2 ? 1.4 : 1.2,
                )
                index == 2 && get!(entries, name, (handle, label))
            end
            column == 2 && text!(
                axis,
                0.97,
                0.06;
                text = case.notation,
                fontsize = 8.5,
                align = (:right, :bottom),
                space = :relative,
            )
            # One range for every panel. The systems sweep different heavy-mass windows, and
            # letting each autoscale while only the bottom row carries tick labels would leave
            # the upper panels read against the wrong axis.
            xlims!(axis, CHARGE_MASS_RANGE...)
        end
    end

    present = [entries[name] for (name, _, _, _) in sources if haskey(entries, name)]
    Legend(
        figure[1, 1:2],
        first.(present),
        last.(present);
        orientation = :horizontal,
        framevisible = false,
        padding = (0, 0, 0, 0),
        colgap = 14,
        patchsize = (16, 2),
        labelsize = 7.5,
    )
    rowgap!(figure.layout, 1, 4)
    for gap in 2:(size(figure.layout)[1] - 1)
        rowgap!(figure.layout, gap, 6)
    end
    colgap!(figure.layout, 1, 16)

    path = joinpath(FIGURES, "charge_distribution")
    save(path * ".pdf", figure)
    save(path * ".png", figure; px_per_unit = 3)
    println("wrote ", path, ".pdf")
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
    figure_multiplicity()
    figure_spectrum()
    figure_residual_quantities()
    figure_charge_distribution()
    return nothing
end

main()
