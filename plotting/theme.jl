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

# Colour and marker pairs handed to measured datasets in order of first appearance. Black is
# reserved for the model curve.
const DATASET_COLOURS =
    (VERMILION, BLUISH_GREEN, REDDISH_PURPLE, SKY_BLUE, ORANGE, BLUE, YELLOW)
const DATASET_MARKERS =
    (:circle, :rect, :utriangle, :diamond, :dtriangle, :pentagon, :star5)

const DATASET_ASSIGNMENT = Dict{String, Int}()
const PANEL_SLOTS = Set{Int}()
const DATASET_TABLE_EXHAUSTED = Ref(false)

"Start a new panel: the slots taken within it are free again."
begin_panel!() = (empty!(PANEL_SLOTS); nothing)

"""
    dataset_style(author_year) -> (colour, marker)

The colour and marker of a measured dataset, keyed by its author–year label. Slots are taken
per panel, so no two datasets of one panel share a pair; a label keeps the slot it was first
given wherever that slot is still free, so a dataset recurring in another panel looks the same.
Beyond seven datasets in one panel the table is reused from the start, with a single warning.
"""
function dataset_style(author_year::String)
    slot = get(DATASET_ASSIGNMENT, author_year, 0)
    if slot == 0 || slot in PANEL_SLOTS
        free = [s for s in eachindex(DATASET_COLOURS) if !(s in PANEL_SLOTS)]
        if isempty(free)
            DATASET_TABLE_EXHAUSTED[] ||
                @warn "more than $(length(DATASET_COLOURS)) measured \
                datasets in one panel; colours and markers are reused" dataset =
                    author_year
            DATASET_TABLE_EXHAUSTED[] = true
            slot = mod1(length(PANEL_SLOTS) + 1, length(DATASET_COLOURS))
        else
            slot = first(free)
        end
        DATASET_ASSIGNMENT[author_year] = slot
    end
    push!(PANEL_SLOTS, slot)
    return (DATASET_COLOURS[slot], DATASET_MARKERS[slot])
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
uncertainty, zero where the file carries none. Columns are taken by position.
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
        push!(σ, length(fields) >= 3 ? something(tryparse(Float64, fields[3]), 0.0) : 0.0)
    end
    return x, y, σ
end
