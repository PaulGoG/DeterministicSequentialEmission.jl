# Theme, colour semantics and table reader shared by the validation figures.

using CairoMakie, MathTeXEngine

const PUBLICATION_THEME = Theme(
    fonts = (; regular = texfont(:text), bold = texfont(:bold), italic = texfont(:italic)),
    fontsize = 26,
    figure_padding = 10,
    linewidth = 3,
    markersize = 14,
    Axis = (
        spinewidth = 1.5,
        xticklabelsize = 22,
        yticklabelsize = 22,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xgridcolor = (:grey, 0.12),
        ygridcolor = (:grey, 0.12),
        xminorticksvisible = false,
        yminorticksvisible = false,
        xtickalign = 1,
        ytickalign = 1,
    ),
    Scatter = (strokewidth = 1.5,),
    Legend = (framevisible = false, orientation = :horizontal, titlefont = :bold),
)

# In-axis annotations at 0.8 of the base size, and the vertical step between stacked lines of
# them [px].
const ANNOTATION_SIZE = 21
const ANNOTATION_LINE = 36

# Okabe–Ito palette.
const BLACK = RGBf(0.0, 0.0, 0.0)
const ORANGE = RGBf(0.902, 0.624, 0.0)
const SKY_BLUE = RGBf(0.337, 0.706, 0.914)
const BLUISH_GREEN = RGBf(0.0, 0.620, 0.451)
const YELLOW = RGBf(0.941, 0.894, 0.259)
const BLUE = RGBf(0.0, 0.447, 0.698)
const VERMILION = RGBf(0.835, 0.369, 0.0)
const REDDISH_PURPLE = RGBf(0.800, 0.475, 0.655)

# Reference lines and their labels.
const GUIDE_COLOUR = RGBf(0.5, 0.5, 0.5)
const GUIDE_LABEL_COLOUR = RGBf(0.35, 0.35, 0.35)

"""
    SYSTEM_STYLE

Colour and line style of each fissioning system, keyed by its label, identical in every figure.
The line style repeats the distinction the colour makes, so the curves survive greyscale.
"""
const SYSTEM_STYLE = Dict{String, Tuple{RGBf, Any}}(
    "Cf252_sf" => (ORANGE, :solid),
    "U235_nth" => (BLUE, :dash),
    "Pu239_nth" => (VERMILION, :dashdot),
    "U233_nth" => (BLUISH_GREEN, :dot),
    "U235_nres" => (SKY_BLUE, (:dot, :dense)),
)

"""
    DatasetPalette(colours, markers)

Colour and marker pairs handed to measured datasets in order of first appearance, the slot each
author–year label has been given, and the slots taken in the panel being drawn.
"""
struct DatasetPalette
    colours::Vector{RGBf}
    markers::Vector{Symbol}
    assignment::Dict{String, Int}
    taken::Set{Int}
    exhausted::Base.RefValue{Bool}
end

function DatasetPalette(colours, markers)
    length(colours) == length(markers) || throw(
        DimensionMismatch("$(length(colours)) colours against $(length(markers)) markers"),
    )
    return DatasetPalette(
        collect(RGBf, colours),
        collect(Symbol, markers),
        Dict{String, Int}(),
        Set{Int}(),
        Ref(false),
    )
end

# The measured datasets of the multiplicity, energy and spectrum figures. Black is reserved for
# the model curve.
const DATASET_PALETTE = DatasetPalette(
    (VERMILION, BLUISH_GREEN, REDDISH_PURPLE, SKY_BLUE, ORANGE, BLUE, YELLOW),
    (:circle, :rect, :utriangle, :diamond, :dtriangle, :pentagon, :star5),
)

# The measured datasets of the product figures, which draw the evaluations in vermilion and
# blue: those two are left out, and a dark purple outside the Okabe–Ito set takes the sixth slot.
const DARK_PURPLE = RGBf(0.42, 0.24, 0.6)
const PRODUCT_PALETTE = DatasetPalette(
    (BLUISH_GREEN, REDDISH_PURPLE, SKY_BLUE, ORANGE, YELLOW, DARK_PURPLE),
    (:circle, :rect, :utriangle, :diamond, :dtriangle, :pentagon),
)

"Start a new panel: the slots taken within it are free again."
begin_panel!() =
    (foreach(p -> empty!(p.taken), (DATASET_PALETTE, PRODUCT_PALETTE)); nothing)

"""
    dataset_style([palette], author_year) -> (colour, marker)

The colour and marker of a measured dataset, keyed by its author–year label, from `palette`,
`DATASET_PALETTE` by default. Slots are taken per panel, so no two datasets of one panel share a
pair; a label keeps the slot it was first given wherever that slot is still free, so a dataset
recurring in another panel looks the same. Beyond the length of the palette in one panel the
table is reused from the start, with a single warning.
"""
dataset_style(author_year::String) = dataset_style(DATASET_PALETTE, author_year)

function dataset_style(palette::DatasetPalette, author_year::String)
    slot = get(palette.assignment, author_year, 0)
    if slot == 0 || slot in palette.taken
        free = [s for s in eachindex(palette.colours) if !(s in palette.taken)]
        if isempty(free)
            palette.exhausted[] || @warn "more than $(length(palette.colours)) measured \
                datasets in one panel; colours and markers are reused" dataset =
                author_year
            palette.exhausted[] = true
            slot = mod1(length(palette.taken) + 1, length(palette.colours))
        else
            slot = first(free)
        end
        palette.assignment[author_year] = slot
    end
    push!(palette.taken, slot)
    return (palette.colours[slot], palette.markers[slot])
end

"""
    darker(colour) -> RGBf

The same hue at 70 % of its HSV value, for marker strokes. Scaling every RGB channel by one
factor leaves hue and saturation unchanged.
"""
darker(colour::RGBf) = RGBf(0.7 * colour.r, 0.7 * colour.g, 0.7 * colour.b)

"""
    read_columns(path) -> (x, y, σ)

Read a two- or three-column table with a one-line header: abscissa, value, and an optional
uncertainty, zero where the file carries none or writes `NaN` (not quoted). Columns are taken
by position.
"""
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
        u = length(fields) >= 3 ? something(tryparse(Float64, fields[3]), NaN) : NaN
        push!(σ, isfinite(u) ? u : 0.0)
    end
    return x, y, σ
end

"""
    read_rows(path) -> Vector{Vector{Float64}}

The rows of a whitespace-separated table with a one-line header, every field parsed by position,
`NaN` where a field is not a number.
"""
function read_rows(path)
    rows = Vector{Vector{Float64}}()
    for (index, line) in enumerate(eachline(path))
        index == 1 && continue
        fields = split(strip(line))
        isempty(fields) && continue
        push!(rows, [something(tryparse(Float64, field), NaN) for field in fields])
    end
    return rows
end
