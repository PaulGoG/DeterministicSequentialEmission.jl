# Entry point: julia scripts/stage_product_yields.jl [data directory] [--library <library>]
#              [--system <system>] [--file <library>=<path>] [--downloads <directory>]
#
# Stages the evaluated independent fission product yields of ENDF/B-VIII.0 and JEFF-3.3 for
# ²⁵²Cf(sf) and for ²³³U, ²³⁵U and ²³⁹Pu at thermal incident energy: the post-neutron product
# yields the model's own are compared with. The ENDF-6 sublibraries come from their
# distributors, ENDF/B-VIII.0 NFY and SFY as the zip archives of the NNDC
# (https://www.nndc.bnl.gov/endf-b8.0/download.html), JEFF-3.3 NFY and SFY as the two ENDF-6
# tapes of the NEA Data Bank record doi:10.82555/nhfqy-hvz09.
#
# From MF8 MT454 of each material the script takes the thermal point, 0.0253 eV, or the single
# point of spontaneous fission; sums the isomeric states of each (Z, A_p), and sums over Z for
# the mass yield, adding uncertainties in quadrature. The ternary particles JEFF-3.3 lists among
# the products are not fragments and are left out. Yields are written in per cent per fission,
# 200 over both products, to data/<system>/Y_vs_A_p/<library>.dat (A_p Y Y_uncertainty)
# and to data/<system>/Y_vs_Z_A_p/<library>.dat (Z A_p Y Y_uncertainty), <library> being
# ENDF-B-VIII.0 or JEFF-3.3: the directories and headers ExforFissionData.jl writes for the
# measured yields against the same abscissae. A provenance.toml in each directory records, per
# file, the source, its SHA-256, the material, the energy point and the handling of isomers and
# ternary particles.
#
# A table is never overwritten: one that differs from the new one moves aside as
# <library>#1.dat, <library>#2.dat, …, and its provenance entry moves with it; an identical one
# is left in place. A download that fails is reported with its URL and HTTP status, and the
# other libraries are staged regardless. The file can then be fetched in a browser and passed
# with --file <library>=<path>, a zip archive or an ENDF-6 tape, repeated as needed. Downloads
# go to a temporary directory unless --downloads names one to keep them in.

include(joinpath(@__DIR__, "..", "activate.jl"))

using DrWatson: datadir, gitdescribe, projectdir
using Dates: Dates
using Downloads: Downloads
using Printf: @printf
using SHA: sha256
using TOML: TOML
using p7zip_jll: p7zip

const usage = """
usage: julia scripts/stage_product_yields.jl [data directory]
           [--library ENDF-B-VIII.0|JEFF-3.3] [--system Cf252_sf|U235_nth|Pu239_nth|U233_nth]
           [--file <library>=<path>] [--downloads <directory>]
--library and --system restrict the staging and may be repeated; --file names a local copy of a
library's archive or tape and may be repeated."""

# The thermal point of an evaluated yield file, in eV.
const thermal_energy = 0.0253
# Products lighter than this are ternary particles. JEFF-3.3 lists light nuclei up to A_p = 21
# beside the binary products, whose lightest mass in the four materials is 50.
const ternary_mass_limit = 40
# An ENDF-6 field holds at most 7 significant digits in exponent form; the per-cent sums are
# written to 10, which keeps every digit the file states and drops the binary residue of the
# scaling.
const significant_digits = 10

# A fissioning system as the data tree names it: the token, the target nuclide, its standard
# ENDF-6 material number, and whether the fission is spontaneous or induced by thermal neutrons.
struct FissioningSystem
    token::String
    Z::Int
    A::Int
    symbol::String
    MAT::Int
    spontaneous::Bool
end

const systems = (
    FissioningSystem("Cf252_sf", 98, 252, "Cf", 9861, true),
    FissioningSystem("U235_nth", 92, 235, "U", 9228, false),
    FissioningSystem("Pu239_nth", 94, 239, "Pu", 9437, false),
    FissioningSystem("U233_nth", 92, 233, "U", 9222, false),
)

target_ZA(system::FissioningSystem) = 1000 * system.Z + system.A
compound_mass(system::FissioningSystem) = system.A + (system.spontaneous ? 0 : 1)
sublibrary(system::FissioningSystem) = system.spontaneous ? "sfy" : "nfy"
# NSUB of MF1 MT451 is 10 IPART + ITYPE: 5 for spontaneous, 11 for neutron-induced fission
# yields.
sublibrary_number(system::FissioningSystem) = system.spontaneous ? 5 : 11

const libraries = (
    (
        token = "ENDF-B-VIII.0",
        name = "ENDF/B-VIII.0",
        release = "2018-02-02",
        citation = "D. A. Brown et al., Nucl. Data Sheets 148, 1 (2018)",
        doi = "10.1016/j.nds.2018.02.001",
        dataset_doi = "",
    ),
    (
        token = "JEFF-3.3",
        name = "JEFF-3.3",
        release = "2017-11-20",
        citation = "A. J. M. Plompen et al., Eur. Phys. J. A 56, 181 (2020)",
        doi = "10.1140/epja/s10050-020-00141-9",
        dataset_doi = "10.82555/nhfqy-hvz09",
    ),
)

# Where a library distributes a material: the URL of the archive and, for a zip archive, the
# member that holds the material. JEFF-3.3 distributes each sublibrary as one ENDF-6 tape, which
# is then its own member.
function source(library, system::FissioningSystem)
    s = sublibrary(system)
    if library.token == "ENDF-B-VIII.0"
        directory = "ENDF-B-VIII.0_$(s)"
        file = "$(s)-$(lpad(system.Z, 3, '0'))_$(system.symbol)_$(system.A).endf"
        url = "https://www.nndc.bnl.gov/endf-b8.0/zips/$(directory).zip"
        return (url = url, member = "$(directory)/$(file)")
    end
    url = "https://data.oecd-nea.org/records/nhfqy-hvz09/files/JEFF33-$(s).asc?download=1"
    return (url = url, member = archive_name(url))
end

archive_name(url::AbstractString) = String(first(split(basename(url), '?')))

# A source that could not be had: a refused download, or a file without the material.
struct SourceError <: Exception
    message::String
end

Base.showerror(io::IO, exception::SourceError) = print(io, exception.message)

checksum(path::AbstractString) = bytes2hex(open(sha256, path))
utc(time::Dates.DateTime) = string(round(time, Dates.Second)) * "Z"

function retrieve(library, url::AbstractString, directory::AbstractString)
    path = joinpath(directory, archive_name(url))
    try
        Downloads.download(url, path)
    catch exception
        exception isa Downloads.RequestError || rethrow()
        detail = isempty(exception.message) ? "" : " ($(exception.message))"
        throw(
            SourceError(
                "$(url): HTTP status $(exception.response.status)$(detail); fetch it in a \
                 browser and pass --file $(library.token)=<path>",
            ),
        )
    end
    return (
        path = path,
        sha256 = checksum(path),
        retrieved = utc(Dates.now(Dates.UTC)),
        retrieval = "download",
    )
end

function local_source(path::AbstractString)
    return (
        path = path,
        sha256 = checksum(path),
        retrieved = utc(Dates.unix2datetime(mtime(path))),
        retrieval = "local file; retrieved_utc is its modification time",
    )
end

is_zip(path::AbstractString) =
    open(io -> read(io, 4) == UInt8[0x50, 0x4b, 0x03, 0x04], path)

# The ENDF-6 text of a source: the zip member, read through 7-Zip without unpacking the
# archive, or the whole file. A zip without the member gives an empty text.
function endf_text(path::AbstractString, member::AbstractString)
    is_zip(path) || return read(path, String)
    command = `$(p7zip()) e -so -bso0 -bsp0 -bse0 -- $(path) $(member)`
    return read(pipeline(command; stderr = devnull), String)
end

# Columns first to last of an 80-column record, by byte, blank where the line is short.
function columns(line::AbstractString, first::Int, last::Int)
    units = codeunits(line)
    first > length(units) && return ""
    return String(units[first:min(last, length(units))])
end

field(line::AbstractString, k::Int) = columns(line, 11k - 10, 11k)

function control(line::AbstractString)
    MAT = tryparse(Int, strip(columns(line, 67, 70)))
    MF = tryparse(Int, strip(columns(line, 71, 72)))
    MT = tryparse(Int, strip(columns(line, 73, 75)))
    (MAT === nothing || MF === nothing || MT === nothing) && return nothing
    return (MAT, MF, MT)
end

# An ENDF-6 real: the Fortran forms with an exponent letter, and the compact form without one
# that the format uses to fit seven digits into eleven columns (9.223300+4, 1.5-5). Blank is 0.
function endf_real(text::AbstractString)
    s = replace(strip(text), " " => "", 'D' => 'E', 'd' => 'e')
    isempty(s) && return 0.0
    value = tryparse(Float64, s)
    value === nothing || return value
    m = match(r"^([+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+))([+-][0-9]+)$", s)
    m === nothing && throw(ArgumentError("not an ENDF-6 real: \"$(text)\""))
    return parse(Float64, string(m[1], "e", m[2]))
end

function endf_integer(text::AbstractString)
    s = strip(text)
    return isempty(s) ? 0 : parse(Int, s)
end

# The material of a text whose MF1 MT451 names the target in its ground state and the
# sublibrary: its MAT, the NVER and LREL it states, and the records of its MF8 MT454; nothing
# where the text holds no such material.
function locate_material(text::AbstractString, system::FissioningSystem)
    headers = Dict{Int, Vector{String}}()
    for line in eachline(IOBuffer(text))
        record = control(line)
        (record === nothing || record[2:3] != (1, 451)) && continue
        lines = get!(headers, record[1], String[])
        length(lines) < 3 && push!(lines, line)
    end
    matches = [
        MAT for (MAT, lines) in headers if length(lines) == 3 &&
            round(Int, endf_real(field(lines[1], 1))) == target_ZA(system) &&
            endf_integer(field(lines[2], 4)) == 0 &&
            endf_integer(field(lines[3], 5)) == sublibrary_number(system)
    ]
    isempty(matches) && return nothing
    length(matches) == 1 || throw(
        ArgumentError(
            "$(length(matches)) materials of ZA = $(target_ZA(system)) in one file",
        ),
    )
    MAT = only(matches)
    MAT == system.MAT || throw(
        ArgumentError(
            "ZA = $(target_ZA(system)) is filed as MAT $(MAT), not the standard $(system.MAT)",
        ),
    )
    section = [line for line in eachline(IOBuffer(text)) if control(line) == (MAT, 8, 454)]
    isempty(section) && throw(ArgumentError("MAT $(MAT) carries no MF8 MT454"))
    header = headers[MAT][3]
    return (
        MAT = MAT,
        NVER = endf_integer(field(header, 6)),
        LREL = endf_integer(field(header, 3)),
        section = section,
    )
end

# The LIST records of MF8 MT454, one per incident energy: (ZAFP, FPS, YI, DYI) per product, the
# yields as fractions per fission.
function independent_yields(section::Vector{String}, system::FissioningSystem)
    fields = [field(line, k) for line in section for k in 1:6]
    ZA = round(Int, endf_real(fields[1]))
    ZA == target_ZA(system) ||
        throw(ArgumentError("MF8 MT454 is of ZA = $(ZA), not $(target_ZA(system))"))
    energies = endf_integer(fields[3])
    points = NamedTuple{(:energy, :products), Tuple{Float64, Vector{NTuple{4, Float64}}}}[]
    position = 7
    for _ in 1:energies
        position + 5 <= length(fields) ||
            throw(ArgumentError("MF8 MT454 ends before its $(energies) LIST records"))
        energy = endf_real(fields[position])
        NN = endf_integer(fields[position + 4])
        NFP = endf_integer(fields[position + 5])
        NN == 4NFP || throw(ArgumentError("LIST at $(energy) eV: NN = $(NN), NFP = $(NFP)"))
        position += 6
        position + NN - 1 <= length(fields) ||
            throw(ArgumentError("LIST at $(energy) eV is truncated"))
        values = map(endf_real, fields[position:(position + NN - 1)])
        products = [Tuple(values[i:(i + 3)]) for i in 1:4:NN]
        push!(points, (energy = energy, products = products))
        position += 6cld(NN, 6)
    end
    return points
end

# The point a comparison takes: the single point of spontaneous fission, the thermal point of
# induced fission.
function select_point(points, system::FissioningSystem)
    energies = [point.energy for point in points]
    if system.spontaneous
        length(points) == 1 ||
            throw(ArgumentError("spontaneous fission with energy points $(energies) eV"))
        return only(points)
    end
    index = findfirst(E -> isapprox(E, thermal_energy; rtol = 1e-6), energies)
    index === nothing &&
        throw(ArgumentError("no thermal point among the energy points $(energies) eV"))
    return points[index]
end

# Per-cent yields per (Z, A_p), ordered by Z and then A_p: the isomeric states summed, their
# uncertainties in quadrature (MT454 states no correlations), the ternary particles counted
# apart.
function charge_resolved(products)
    yields = Dict{Tuple{Int, Int}, Tuple{Float64, Float64, Int}}()
    ternary_entries = 0
    ternary_yield = 0.0
    for (ZAFP, FPS, YI, DYI) in products
        (isinteger(ZAFP) && isinteger(FPS)) ||
            throw(ArgumentError("ZAFP = $(ZAFP), FPS = $(FPS) is not a nuclide state"))
        Z, A = divrem(Int(ZAFP), 1000)
        1 <= Z <= A || throw(ArgumentError("ZAFP = $(ZAFP) is not a nuclide"))
        if A < ternary_mass_limit
            ternary_entries += 1
            ternary_yield += YI
            continue
        end
        Y, variance, states = get(yields, (Z, A), (0.0, 0.0, 0))
        yields[(Z, A)] = (Y + YI, variance + DYI^2, states + 1)
    end
    rows = sort!([(Z, A, 100Y, 100sqrt(variance)) for ((Z, A), (Y, variance, _)) in yields])
    return (
        rows = rows,
        isomeric = count(entry -> entry[3] > 1, values(yields)),
        ternary_entries = ternary_entries,
        ternary_yield = 100ternary_yield,
    )
end

function mass_resolved(rows)
    yields = Dict{Int, Tuple{Float64, Float64}}()
    for (_, A, Y, ΔY) in rows
        y, variance = get(yields, A, (0.0, 0.0))
        yields[A] = (y + Y, variance + ΔY^2)
    end
    return sort!([(A, Y, sqrt(variance)) for (A, (Y, variance)) in yields])
end

# Mean light and heavy product masses, the groups divided at half the compound nucleus mass; a
# mass at the division counts half to each.
function mean_masses(masses, A₀::Int)
    weight(A) = A < A₀ / 2 ? 1.0 : A == A₀ / 2 ? 0.5 : 0.0
    light = sum(weight(A) * Y for (A, Y, _) in masses)
    heavy = sum((1 - weight(A)) * Y for (A, Y, _) in masses)
    return (
        light = sum(weight(A) * Y * A for (A, Y, _) in masses) / light,
        heavy = sum((1 - weight(A)) * Y * A for (A, Y, _) in masses) / heavy,
    )
end

format(x::Integer) = string(x)
format(x::Real) = string(round(x; sigdigits = significant_digits))

function table_text(header::AbstractString, rows)
    return join([header; [join(map(format, row), ' ') for row in rows]], '\n') * '\n'
end

# Never overwrite a table. An existing file moves to <name>#1.dat, <name>#2.dat, …: the
# numbering scripts/run.jl applies to a run directory, kept before the extension.
function move_aside(file::AbstractString)
    stem, extension = splitext(file)
    index = 1
    while ispath("$(stem)#$(index)$(extension)")
        index += 1
    end
    destination = "$(stem)#$(index)$(extension)"
    mv(file, destination)
    return destination
end

# Writes one table and its provenance entry. An identical table already in place is kept with
# its record; a different one is moved aside, and its entry with it.
function stage_table(directory::AbstractString, name::AbstractString, text, record)
    mkpath(directory)
    file = joinpath(directory, name)
    path = joinpath(directory, "provenance.toml")
    provenance = isfile(path) ? TOML.parsefile(path) : Dict{String, Any}()
    if isfile(file)
        read(file, String) == text && return "unchanged"
        displaced = basename(move_aside(file))
        haskey(provenance, name) && (provenance[displaced] = pop!(provenance, name))
    end
    write(file, text)
    provenance[name] = record
    partial = path * ".partial"
    open(partial, "w") do io
        println(
            io,
            "# Provenance of the evaluated product yields in this directory, per file.",
        )
        TOML.print(io, provenance; sorted = true)
    end
    mv(partial, path; force = true)
    return "written"
end

# The source of a material: a local file given for the library that holds it, or else the
# library's archive, downloaded once per run.
function obtain(library, system::FissioningSystem, files, downloads, retrieved)
    (; url, member) = source(library, system)
    given = get(files, library.token, String[])
    for path in given
        material = locate_material(endf_text(path, member), system)
        material === nothing || return (origin = local_source(path), url, member, material)
    end
    isempty(given) || throw(
        SourceError(
            "no file given for $(library.token) holds MAT $(system.MAT) of $(system.token); \
             its source is $(url)",
        ),
    )
    origin = get!(() -> retrieve(library, url, downloads), retrieved, url)
    material = locate_material(endf_text(origin.path, member), system)
    material === nothing && throw(
        SourceError(
            "$(url) returned no MAT $(system.MAT) with sublibrary \
             $(sublibrary_number(system)); if the server sent a page instead of the file, \
             fetch it in a browser and pass --file $(library.token)=<path>",
        ),
    )
    return (origin = origin, url, member, material)
end

function stage(data::AbstractString, library, system::FissioningSystem, obtained, commit)
    (; origin, url, member, material) = obtained
    points = independent_yields(material.section, system)
    point = select_point(points, system)
    charges = charge_resolved(point.products)
    masses = mass_resolved(charges.rows)
    total = sum(Y for (_, Y, _) in masses)
    common = Dict{String, Any}(
        "library" => library.name,
        "sublibrary" => uppercase(sublibrary(system)),
        "release" => library.release,
        "NVER" => material.NVER,
        "LREL" => material.LREL,
        "url" => url,
        "archive" => basename(origin.path),
        "member" => member,
        "retrieval" => origin.retrieval,
        "sha256" => origin.sha256,
        "retrieved_utc" => origin.retrieved,
        "MAT" => material.MAT,
        "ZA" => target_ZA(system),
        "system" => system.token,
        "section" => "MF8 MT454, independent yields",
        "incident_energy_ev" => point.energy,
        "energy_points_in_file_ev" => [p.energy for p in points],
        "MT454_entries" => length(point.products),
        "isomer_handling" => "the isomeric states (FPS) of each (Z, A_p) are summed, their \
                              uncertainties added in quadrature",
        "isomeric_entries_summed" => charges.isomeric,
        "ternary_handling" => "products with A_p < $(ternary_mass_limit) are ternary \
                               particles and are not written",
        "ternary_entries_excluded" => charges.ternary_entries,
        "ternary_Y_excluded" =>
            round(charges.ternary_yield; sigdigits = significant_digits),
        "units" => "per cent per fission, 200 over both products",
        "sum_Y" => round(total; sigdigits = significant_digits),
        "citation" => library.citation,
        "doi" => library.doi,
        "script" => relpath(@__FILE__, projectdir()),
        "commit" => commit,
    )
    isempty(library.dataset_doi) || (common["dataset_doi"] = library.dataset_doi)

    zero_uncertainty(rows) = count(row -> iszero(last(row)), rows)
    by_charge = merge(
        common,
        Dict{String, Any}(
            "columns" => "Z A_p Y Y_uncertainty",
            "rows" => length(charges.rows),
            "uncertainty_zero_rows" => zero_uncertainty(charges.rows),
        ),
    )
    by_mass = merge(
        common,
        Dict{String, Any}(
            "columns" => "A_p Y Y_uncertainty",
            "rows" => length(masses),
            "uncertainty_zero_rows" => zero_uncertainty(masses),
            "charge_handling" => "Y(A_p) is the sum over Z of Y(Z, A_p), the uncertainties \
                                  added in quadrature",
        ),
    )
    name = "$(library.token).dat"
    directory = joinpath(data, system.token)
    status_charge = stage_table(
        joinpath(directory, "Y_vs_Z_A_p"),
        name,
        table_text("Z A_p Y Y_uncertainty", charges.rows),
        by_charge,
    )
    status_mass = stage_table(
        joinpath(directory, "Y_vs_A_p"),
        name,
        table_text("A_p Y Y_uncertainty", masses),
        by_mass,
    )

    means = mean_masses(masses, compound_mass(system))
    largest = sort(masses; by = row -> row[2], rev = true)[1:min(5, length(masses))]
    @printf(
        "  %s, %s: MAT %d at %g eV, %d MT454 entries, %d isomeric (Z, A_p) summed\n",
        system.token,
        library.token,
        material.MAT,
        point.energy,
        length(point.products),
        charges.isomeric
    )
    @printf(
        "    %d (Z, A_p) entries, %d masses, sum of Y(A_p) = %.6f %%\n",
        length(charges.rows),
        length(masses),
        total
    )
    charges.ternary_entries > 0 && @printf(
        "    %d ternary particles left out, %.6f %%\n",
        charges.ternary_entries,
        charges.ternary_yield
    )
    println(
        "    largest Y(A_p): ",
        join(["$(round(Y; digits = 3)) % at $(A)" for (A, Y, _) in largest], ", "),
    )
    @printf(
        "    <A_L> = %.3f, <A_H> = %.3f, sum %.3f\n",
        means.light,
        means.heavy,
        means.light + means.heavy
    )
    println("    ", status_mass, " ", joinpath(directory, "Y_vs_A_p", name))
    println("    ", status_charge, " ", joinpath(directory, "Y_vs_Z_A_p", name))
    return nothing
end

# The items whose tokens were selected; all of them where none was.
chosen(items, tokens) = [item for item in items if isempty(tokens) || item.token in tokens]

function parse_arguments(arguments::Vector{String})
    data = nothing
    selected_libraries = String[]
    selected_systems = String[]
    files = Dict{String, Vector{String}}()
    downloads = nothing
    library_tokens = [library.token for library in libraries]
    system_tokens = [system.token for system in systems]
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument in ("-h", "--help")
            return nothing
        elseif argument in ("--library", "--system", "--file", "--downloads")
            index < length(arguments) || throw(ArgumentError("$(argument) takes a value"))
            value = arguments[index + 1]
            if argument == "--library"
                value in library_tokens || throw(
                    ArgumentError(
                        "--library $(value): one of $(join(library_tokens, ", "))",
                    ),
                )
                push!(selected_libraries, value)
            elseif argument == "--system"
                value in system_tokens || throw(
                    ArgumentError("--system $(value): one of $(join(system_tokens, ", "))"),
                )
                push!(selected_systems, value)
            elseif argument == "--file"
                token, path = occursin('=', value) ? split(value, '='; limit = 2) : ("", "")
                token in library_tokens ||
                    throw(ArgumentError("--file $(value): expected <library>=<path>"))
                isfile(path) || throw(ArgumentError("--file $(value): no such file"))
                push!(get!(files, String(token), String[]), abspath(path))
            else
                downloads = abspath(value)
            end
            index += 2
        elseif startswith(argument, "-")
            throw(ArgumentError("unknown option $(argument)"))
        else
            data === nothing || throw(ArgumentError("more than one data directory given"))
            isdir(argument) || throw(ArgumentError("no data directory $(argument)"))
            data = abspath(argument)
            index += 1
        end
    end
    return (
        data = something(data, datadir()),
        libraries = chosen(libraries, selected_libraries),
        systems = chosen(systems, selected_systems),
        files = files,
        downloads = downloads,
    )
end

function main(arguments::Vector{String})
    options = try
        parse_arguments(arguments)
    catch exception
        exception isa ArgumentError || rethrow()
        println(stderr, exception.msg)
        println(stderr, usage)
        return 1
    end
    if options === nothing
        println(usage)
        return 0
    end
    downloads = if options.downloads === nothing
        mktempdir()
    else
        mktempdir(mkpath(options.downloads); prefix = "retrieval_", cleanup = false)
    end
    # The commit of the script, with -dirty where the checkout has uncommitted changes.
    commit = something(gitdescribe(projectdir()), "")
    failures = 0
    for library in options.libraries
        retrieved = Dict{String, Any}()
        for system in options.systems
            obtained = try
                obtain(library, system, options.files, downloads, retrieved)
            catch exception
                exception isa SourceError || rethrow()
                failures += 1
                println(stderr, "  $(system.token), $(library.token): $(exception.message)")
                continue
            end
            stage(options.data, library, system, obtained, commit)
        end
    end
    failures > 0 && println(stderr, "  ", failures, " material(s) not staged")
    return failures == 0 ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
