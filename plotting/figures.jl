# Validation figures: the model against measurement.
#
#   julia --project=plotting plotting/figures.jl [experimental directory]
#
# The model side needs the input data described in data/README.md, including the temperature
# ratios under data/partition. The measured side needs fission observables retrieved from EXFOR
# with ExforFissionData.jl, laid out as
#
#   <experimental>/multiplicity/<system>/*.dat     A ν [σν]
#   <experimental>/spectrum_ratio/<system>/*.dat   E N(E)/M(E,T_M) [σ]
#
# which is the layout that package writes. The directory defaults to data/experimental and is
# not version-controlled; nothing here knows where yours lives.

include(joinpath(@__DIR__, "activate.jl"))

using DeterministicSequentialEmission
using CairoMakie, MathTeXEngine
using Measurements: value
using QuadGK: quadgk
using Printf

const ROOT = dirname(@__DIR__)
const DATA = joinpath(ROOT, "data")
const FIGURES = joinpath(ROOT, "figures")
const EXPERIMENTAL =
    length(ARGS) >= 1 ? ARGS[1] : joinpath(DATA, "experimental")

"""
One fissioning system, with everything needed to run it and everything to compare it against.
`sweep` matches the kinetic energy grid of the measured yield, so the whole distribution is
reached rather than every other point of it. `measured` is the evaluated total multiplicity per
fission, ENDF/B-VIII.0 for the induced systems and the Axton standard for 252-Cf.
"""
struct Case
    label::Any
    system::FissioningSystem
    polarization::String
    yield::String
    heavy_masses::UnitRange{Int}
    sweep::StepRangeLen{Float64}
    identifier::String
    measured::Float64
end

const CASES = [
    Case(
        L"$^{252}$Cf(sf)",
        spontaneous_fission(Nuclide(98, 252)),
        "Cf252.dat",
        joinpath("Cf252", "vespa_Y_A_TKE.dat"),
        126:174,
        130.0:1.0:230.0,
        "Cf252_0f",
        3.7573,
    ),
    Case(
        L"$^{235}$U(n$_\mathrm{th}$,f)",
        neutron_induced_fission(Nuclide(92, 235), 2.53e-8),
        "U235.dat",
        joinpath("U235", "Y_A_TKE.dat"),
        118:160,
        100.0:1.0:200.0,
        "U235_nf",
        2.4140,
    ),
    Case(
        L"$^{239}$Pu(n$_\mathrm{th}$,f)",
        neutron_induced_fission(Nuclide(94, 239), 2.53e-8),
        "Pu239.dat",
        joinpath("Pu239", "Y_A_TKE.dat"),
        120:159,
        130.0:1.0:210.0,
        "Pu239_nf",
        2.8750,
    ),
]

# One common abscissa, so the three panels can be read against each other and the inner tick
# labels dropped. Nothing is plotted outside it.
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

"Read a two- or three-column measured file: abscissa, value, optional uncertainty."
function read_measured(path)
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

measured_sets(kind, system) = begin
    directory = joinpath(EXPERIMENTAL, kind, system)
    isdir(directory) || return Tuple{String, String}[]
    files = sort(filter(name -> endswith(name, ".dat"), readdir(directory)))
    [(_author(name), joinpath(directory, name)) for name in files]
end

# Files are named <accession>_<author>_<year>.dat by the retrieval.
function _author(name)
    parts = split(replace(name, ".dat" => ""), "_")
    length(parts) >= 3 || return name
    author = replace(parts[2], "." => ". ")
    author = replace(author, "Goeoek" => "Göök", "Budtz-jorgensen" => "Budtz-Jørgensen")
    return strip(last(split(author, ". "))) * " " * parts[3]
end

"""
Everything one case needs, assembled once: the mass table, the fragmentation domain, the yield
weights, and the TXE partition. Assembling it separately from the runs keeps the two figures on
identical inputs.
"""
function setup(case::Case)
    table = read_mass_table(joinpath(DATA, "mass_excess", "ame2020.dat"))
    polarization =
        read_charge_polarization(joinpath(DATA, "charge_distribution", case.polarization))
    domain = fragmentation_domain(case.system, polarization, case.heavy_masses)
    model = BackShiftedFermiGas(table)
    yields = build_fragment_yields(
        domain,
        polarization,
        read_experimental_yield(joinpath(DATA, "yield", case.yield)),
    )
    # The TXE partition is the one FissionTemperatureRatio.jl extracts for this system from its
    # own measured ν(A): a temperature ratio R_T against heavy fragment mass. Each system is
    # therefore partitioned by its own data rather than by a ratio carried over from another.
    segments = read_segments(
        joinpath(DATA, "partition", case.identifier, "R_T_systematic_trend.csv"),
    )
    method = TemperatureRatio(segments, model)
    return (
        table = table,
        domain = domain,
        model = model,
        cross_section = VariableCrossSection(),
        yields = yields,
        weights = weights_from(yields),
        method = method,
        segments = segment_coverage(method, case.heavy_masses),
        reached = coverage(yields, case.sweep).reached,
    )
end

"Run one case and return ν(A) together with the multiplicity per fragment."
function multiplicity(case::Case, context = setup(case))
    ν = ConfigurationValues{Float64}()
    for entry in context.domain.entries, TKE in case.sweep
        balance = energetics(context.table, case.system, entry, TKE)
        balance === nothing && continue
        split = partition(context.method, entry, value(balance.TXE))
        split === nothing && continue
        for (fragment, E) in ((entry.light, split.light), (entry.heavy, split.heavy))
            chain = emission_chain(
                context.table,
                context.model,
                context.cross_section,
                fragment,
                E,
            )
            ν[(fragment.A, fragment.Z, TKE)] = neutron_multiplicity(chain)
        end
    end
    return average(ν, context.weights; over = (:charge, :kinetic_energy)),
    total_average(ν, context.weights)
end

function figure_multiplicity()
    set_theme!(THEME)
    figure = Figure(size = (510, 560))
    summary = String[]

    for (row, case) in pairs(CASES)
        context = setup(case)
        curve, total = multiplicity(case, context)
        push!(
            summary,
            @sprintf(
                "%-22s ⟨ν⟩_pair = %.3f   measured %.3f   %+.1f %%   R_T covers %.0f %%   \
                 yield reached %.0f %%",
                string(case.label),
                2 * total,
                case.measured,
                100 * (2 * total - case.measured) / case.measured,
                100 * context.segments.covered,
                100 * context.reached
            )
        )

        axis = Axis(
            figure[2 * row, 1];
            ylabel = L"\nu(A)",
            xlabel = row == length(CASES) ? L"Fragment mass $A$" : "",
            xticks = 80:10:170,
            xticklabelsvisible = row == length(CASES),
        )

        handles = []
        labels = Any[]
        for (index, (author, path)) in pairs(measured_sets("multiplicity", case.identifier))
            index > length(PALETTE) && break
            A, ν, σ = read_measured(path)
            keep = first(MASS_RANGE) .<= A .<= last(MASS_RANGE)
            handle = scatter!(
                axis,
                A[keep],
                ν[keep];
                color = (PALETTE[index], 0.7),
                marker = MARKERS[index],
                markersize = 4,
            )
            any(>(0), σ[keep]) && errorbars!(
                axis,
                A[keep],
                ν[keep],
                σ[keep];
                color = (PALETTE[index], 0.3),
                linewidth = 0.5,
            )
            push!(handles, handle)
            push!(labels, author)
        end
        model = lines!(axis, first.(curve), last.(curve); color = :black, linewidth = 1.6)
        push!(handles, model)
        push!(labels, "This model")

        symmetric = case.system.compound.A / 2
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
    # wider separation between the three systems.
    for row in 1:length(CASES)
        rowgap!(figure.layout, 2 * row - 1, 2)
        row < length(CASES) && rowgap!(figure.layout, 2 * row, 8)
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

function figure_spectrum()
    set_theme!(THEME)
    case = first(CASES)
    context = setup(case)
    cross_section = context.cross_section

    records = []
    total_weight = 0.0
    for entry in context.domain.entries, TKE in 130.0:2.0:230.0
        y = get(context.yields.values, (entry.heavy.A, entry.heavy.Z, TKE), nothing)
        (y === nothing || y <= 0) && continue
        balance = energetics(context.table, case.system, entry, TKE)
        balance === nothing && continue
        split = partition(context.method, entry, value(balance.TXE))
        split === nothing && continue
        light =
            emission_chain(context.table, context.model, cross_section, entry.light, split.light)
        heavy =
            emission_chain(context.table, context.model, cross_section, entry.heavy, split.heavy)
        w = pair_weights(MultiplicityWeighting(), light, heavy, split)
        w === nothing && continue
        E_light, E_heavy = fragment_energy_per_nucleon(entry, TKE)
        push!(records, (y, light, heavy, w, E_light, E_heavy))
        total_weight += y
    end

    spectrum(E) =
        sum(
            y * (
                w[1] * chain_spectrum(l, cross_section, El, E) +
                w[2] * chain_spectrum(h, cross_section, Eh, E)
            ) for (y, l, h, w, El, Eh) in records
        ) / total_weight

    # The normalisation is the model's own: the spectrum integrates to one neutron, so the ratio
    # is an absolute statement and nothing is fitted to the data below.
    norm, _ = quadgk(spectrum, 0, Inf; rtol = 1e-8)
    moment, _ = quadgk(E -> E * spectrum(E), 0, Inf; rtol = 1e-8)
    T_M = 1.42     # the Maxwellian temperature every retrieved set is reported against
    @printf(
        "⟨E⟩ = %.4f MeV against a Maxwellian at T_M = %.2f MeV; spectrum integrates to %.6f\n",
        moment / norm,
        T_M,
        norm
    )

    grid = 10 .^ range(log10(LOW), log10(HIGH); length = 260)
    ratio = [spectrum(E) / (norm * maxwellian(E, T_M)) for E in grid]

    figure = Figure(size = (510, 300))
    axis = Axis(
        figure[2, 1];
        xscale = log10,
        xlabel = L"Laboratory neutron energy $E$ [MeV]",
        ylabel = L"N(E)\,/\,M(E,\,T_\mathrm{M})",
    )

    handles = []
    labels = Any[]
    for (index, (author, path)) in pairs(measured_sets("spectrum_ratio", case.identifier))
        index > length(PALETTE) && break
        E, R, σ = read_measured(path)
        keep = LOW .<= E .<= HIGH
        any(keep) || continue
        handle = scatter!(
            axis,
            E[keep],
            R[keep];
            color = (PALETTE[index], 0.6),
            marker = MARKERS[index],
            markersize = 3.5,
        )
        any(>(0), σ[keep]) && errorbars!(
            axis,
            E[keep],
            R[keep],
            σ[keep];
            color = (PALETTE[index], 0.25),
            linewidth = 0.5,
        )
        push!(handles, handle)
        push!(labels, author)
    end

    model_line = lines!(axis, grid, ratio; color = :black, linewidth = 1.8)
    push!(handles, model_line)
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

function main()
    mkpath(FIGURES)
    isdir(EXPERIMENTAL) || @warn """
        no experimental directory; the figures will carry the model curves alone. Retrieve the \
        measurements with ExforFissionData.jl and lay them out as described at the top of this \
        file.""" directory = EXPERIMENTAL
    figure_multiplicity()
    figure_spectrum()
    return nothing
end

main()
