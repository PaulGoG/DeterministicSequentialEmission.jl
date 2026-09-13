# Validation figures: the model against measurement.
#
#   julia --project=plotting plotting/figures.jl [results directory] [experimental directory]
#
# Nothing is solved here. The model side is read from what scripts/run.jl wrote under results/,
# so a figure is redrawn from tables rather than by running the model again, and the numbers on
# the page are the numbers in the files. Run the pipeline first:
#
#   julia --project scripts/run.jl config/Cf252_0f.toml
#
# The measured side is retrieved from EXFOR with ExforFissionData.jl and laid out as
#
#   <experimental>/multiplicity/<system>/*.dat     A ν [σν]
#   <experimental>/spectrum_ratio/<system>/*.dat   E N(E)/M(E,T_M) [σ]
#
# which is the layout that package writes. Both directories default to beside the project and
# neither is version-controlled; nothing here knows where yours lives.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie, MathTeXEngine
using Printf

const ROOT = dirname(@__DIR__)
const FIGURES = joinpath(ROOT, "figures")
const RESULTS = length(ARGS) >= 1 ? ARGS[1] : joinpath(ROOT, "results")
const EXPERIMENTAL = length(ARGS) >= 2 ? ARGS[2] : joinpath(ROOT, "data", "experimental")

"""
One fissioning system: the label it is drawn under, the results directory it is read from, and
the evaluated total multiplicity it is compared against — ENDF/B-VIII.0 for the induced systems
and the Axton standard for 252-Cf.
"""
struct Case
    label::Any
    identifier::String
    compound_mass::Int
    measured::Float64
end

const CASES = [
    Case(L"$^{252}$Cf(sf)", "Cf252_0f", 252, 3.7573),
    Case(L"$^{235}$U(n$_\mathrm{th}$,f)", "U235_nf", 236, 2.4140),
    Case(L"$^{239}$Pu(n$_\mathrm{th}$,f)", "Pu239_nf", 240, 2.8750),
    Case(L"$^{233}$U(n$_\mathrm{th}$,f)", "U233_nf", 234, 2.4968),
]

# One common abscissa, so the panels can be read against each other and the inner tick labels
# dropped. Nothing is plotted outside it.
const MASS_RANGE = (76.0, 176.0)

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

"The totals a run wrote, by name."
function read_totals(case::Case)
    path = joinpath(RESULTS, case.identifier, "totals.dat")
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

curve(case::Case, name) = read_columns(joinpath(RESULTS, case.identifier, name))

function measured_sets(kind, system)
    directory = joinpath(EXPERIMENTAL, kind, system)
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
    path = joinpath(RESULTS, case.identifier, "multiplicity_A.dat")
    isfile(path) && return true
    @warn """
        no run to draw; solve it first with `julia --project scripts/run.jl \
        config/$(case.identifier).toml`.""" expected = path
    return false
end

function figure_multiplicity()
    set_theme!(THEME)
    drawable = filter(_available, CASES)
    isempty(drawable) && return nothing

    figure = Figure(size = (510, 30 + 140 * length(drawable)))
    summary = String[]

    for (row, case) in pairs(drawable)
        A, ν, _ = curve(case, "multiplicity_A.dat")
        totals = read_totals(case)
        modelled = get(totals, "multiplicity_per_fission", NaN)
        push!(
            summary,
            @sprintf(
                "%-24s ⟨ν⟩ = %.4f   evaluated %.4f   %+.1f %%",
                string(case.label),
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
        for (index, (author, path)) in pairs(measured_sets("multiplicity", case.identifier))
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
            text = case.label,
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
    path = joinpath(RESULTS, case.identifier, "spectrum_E.dat")
    if !isfile(path)
        @warn "no spectrum to draw" expected = path
        return nothing
    end
    E, N, _ = read_columns(path)
    totals = read_totals(case)
    @printf(
        "⟨E⟩ = %.4f MeV against a Maxwellian at T_M = %.2f MeV\n",
        get(totals, "mean_neutron_laboratory_energy", NaN),
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
    for (index, (author, file)) in pairs(measured_sets("spectrum_ratio", case.identifier))
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
        text = case.label,
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
        ("excitation_A.dat", L"E^*(A)\ [\mathrm{MeV}]"),
        ("residual_energy_A.dat", L"\bar{E}_r(A)\ [\mathrm{MeV}]"),
        ("temperature_A.dat", L"\langle T \rangle(A)\ [\mathrm{MeV}]"),
        ("neutron_energy_A.dat", L"\langle \varepsilon \rangle(A)\ [\mathrm{MeV}]"),
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
            path = joinpath(RESULTS, case.identifier, name)
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
                push!(labels, case.label)
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

function main()
    mkpath(FIGURES)
    isdir(RESULTS) || @warn """
        no results directory; solve the systems first with \
        `julia --project scripts/run.jl config/<case>.toml`.""" directory = RESULTS
    isdir(EXPERIMENTAL) || @warn """
        no experimental directory; the figures will carry the model curves alone. Retrieve the \
        measurements with ExforFissionData.jl and lay them out as described at the top of this \
        file.""" directory = EXPERIMENTAL
    figure_multiplicity()
    figure_spectrum()
    figure_residual_quantities()
    return nothing
end

main()
