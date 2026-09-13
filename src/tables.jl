"""
    TableSpec(columns; delimiter = ' ', skip = 0, optional = ())

How to read one whitespace- or character-delimited data table.

`columns` names the fields in order; `skip` is the number of leading lines to discard, which
is 1 for a file carrying a header row and 0 for one that does not. Names in `optional` may be
absent from the file: a table is accepted if it provides every required column, in order,
followed by any leading prefix of the optional ones.

That last rule is what lets one reader serve both the legacy data files and the output of
`ExforFissionData.jl`, which emits an uncertainty column only when the underlying measurements
carry uncertainties.

# Examples

```jldoctest
julia> TableSpec([:A, :yield, :uncertainty]; skip = 1, optional = (:uncertainty,))
TableSpec(A, yield, [uncertainty]; delimiter=' ', skip=1)
```
"""
struct TableSpec
    columns::Vector{Symbol}
    delimiter::Char
    skip::Int
    optional::Vector{Symbol}

    function TableSpec(
        columns::Vector{Symbol};
        delimiter::AbstractChar = ' ',
        skip::Integer = 0,
        optional::Union{Vector{Symbol}, Tuple{Vararg{Symbol}}} = (),
    )
        names = copy(columns)
        isempty(names) && throw(ArgumentError("a table spec needs at least one column"))
        allunique(names) ||
            throw(ArgumentError("duplicate column names in table spec: $(names)"))
        optional_names = optional isa Vector{Symbol} ? copy(optional) : Symbol[optional...]
        for name in optional_names
            name in names || throw(
                ArgumentError("optional column $name is not among the spec's columns"),
            )
        end
        # Optional columns must form a suffix, otherwise "a leading prefix of the optional
        # ones" is not a well-defined shape.
        required = length(names) - length(optional_names)
        for (index, name) in enumerate(names)
            if index <= required
                name in optional_names && throw(
                    ArgumentError(
                        "optional columns must come last; $name at position $index is optional " *
                        "but precedes a required column",
                    ),
                )
            end
        end
        skip >= 0 || throw(ArgumentError("skip must be non-negative, got $skip"))
        return new(names, Char(delimiter), Int(skip), optional_names)
    end
end

function Base.show(io::IO, spec::TableSpec)
    required = setdiff(spec.columns, spec.optional)
    print(io, "TableSpec(", join(required, ", "))
    isempty(spec.optional) || print(io, ", [", join(spec.optional, ", "), "]")
    print(io, "; delimiter=", repr(spec.delimiter), ", skip=", spec.skip, ")")
    return nothing
end

required_columns(spec::TableSpec) =
    spec.columns[1:(length(spec.columns) - length(spec.optional))]

"""
    DataTable

A table read from disk: `columns` maps each column name to its parsed values, `path` records
where it came from, and `present` lists the columns the file actually supplied.
"""
struct DataTable
    columns::Dict{Symbol, Vector{String}}
    present::Vector{Symbol}
    path::String
    rows::Int

    # Explicit, so Julia generates no `Any`-accepting fallback whose conversion paths would
    # have to be reasoned about.
    function DataTable(
        columns::Dict{Symbol, Vector{String}},
        present::Vector{Symbol},
        path::AbstractString,
        rows::Integer,
    )
        return new(columns, present, String(path), Int(rows))
    end
end

Base.length(table::DataTable) = table.rows
has_column(table::DataTable, name::Symbol) = name in table.present

function Base.show(io::IO, table::DataTable)
    print(
        io,
        "DataTable(",
        table.rows,
        " rows, ",
        join(table.present, ", "),
        " from ",
        basename(table.path),
        ")",
    )
    return nothing
end

"""
    read_table(path, spec) -> DataTable

Read a delimited table according to `spec`, failing with a message that names the file and the
offending line.

Repeated delimiters are treated as one, so the space-padded fixed-width layouts of the legacy
data files read correctly without a separate code path.
"""
function read_table(path::AbstractString, spec::TableSpec)
    isfile(path) || throw(ArgumentError("data file does not exist: $path"))
    lines = readlines(path)
    length(lines) > spec.skip || throw(
        ArgumentError(
            "$path holds $(length(lines)) lines, fewer than the $(spec.skip) skipped",
        ),
    )

    required = required_columns(spec)
    fields_per_row = 0
    raw = Vector{Vector{SubString{String}}}()

    for (offset, line) in enumerate(@view lines[(spec.skip + 1):end])
        stripped = strip(line)
        isempty(stripped) && continue
        fields = filter(!isempty, split(stripped, spec.delimiter))
        if isempty(raw)
            fields_per_row = length(fields)
            if fields_per_row < length(required) || fields_per_row > length(spec.columns)
                throw(
                    ArgumentError(
                        "$path line $(offset + spec.skip): found $fields_per_row fields; " *
                        "spec requires $(length(required)) ($(join(required, ", "))) and " *
                        "accepts up to $(length(spec.columns))",
                    ),
                )
            end
        elseif length(fields) != fields_per_row
            throw(
                ArgumentError(
                    "$path line $(offset + spec.skip): found $(length(fields)) fields, " *
                    "expected $fields_per_row to match the rest of the file",
                ),
            )
        end
        push!(raw, fields)
    end

    isempty(raw) && throw(ArgumentError("$path holds no data rows"))

    present = spec.columns[1:fields_per_row]
    columns = Dict{Symbol, Vector{String}}(
        name => Vector{String}(undef, length(raw)) for name in present
    )
    for (row, fields) in enumerate(raw), (index, name) in enumerate(present)
        columns[name][row] = String(fields[index])
    end
    return DataTable(columns, collect(present), String(path), length(raw))
end

"""
    column(table, name, T) -> Vector{T}

Parse one column as `T`, naming the file, the column and the row on failure.
"""
function column(table::DataTable, name::Symbol, ::Type{T}) where {T}
    has_column(table, name) || throw(
        ArgumentError(
            "$(table.path) has no column $name; it supplied $(join(table.present, ", "))",
        ),
    )
    values = Vector{T}(undef, table.rows)
    for (row, text) in enumerate(table.columns[name])
        parsed = tryparse(T, text)
        parsed === nothing && throw(
            ArgumentError(
                "$(table.path) row $row, column $name: cannot read $(repr(text)) as $T",
            ),
        )
        values[row] = parsed
    end
    return values
end

column(table::DataTable, name::Symbol) = table.columns[name]

"""
    integer_column(table, name) -> Vector{Int}

Parse a nucleon-number column, accepting an integral value written as a float.

Nucleon numbers are integers, but the files that carry them are not all written by hand:
`data/U235_nth/charge_distribution_vs_A.dat` tabulates its mass column as `118.0`, `119.0`, …, because a
generator wrote it out in floating point. Read strictly as `Int` that file cannot be loaded at
all. A value with a fractional part is still an error, and says so.
"""
function integer_column(table::DataTable, name::Symbol)
    has_column(table, name) || throw(
        ArgumentError(
            "$(table.path) has no column $name; it supplied $(join(table.present, ", "))",
        ),
    )
    values = Vector{Int}(undef, table.rows)
    for (row, text) in enumerate(table.columns[name])
        parsed = tryparse(Int, text)
        if parsed === nothing
            approximate = tryparse(Float64, text)
            (approximate === nothing || !isinteger(approximate)) && throw(
                ArgumentError(
                    "$(table.path) row $row, column $name: cannot read $(repr(text)) as a " *
                    "nucleon number",
                ),
            )
            parsed = Int(approximate)
        end
        values[row] = parsed
    end
    return values
end

"""
    column_or(table, name, T, fallback) -> Vector{T}

Parse an optional column, or return a vector of `fallback` when the file omitted it.
"""
function column_or(table::DataTable, name::Symbol, ::Type{T}, fallback::T) where {T}
    has_column(table, name) || return fill(fallback, table.rows)
    return column(table, name, T)
end
