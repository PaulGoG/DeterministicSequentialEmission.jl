# Entry point: julia scripts/run.jl config/<case>.toml [data directory]
#
# Solves every configuration the domain and the kinetic energy grid span, contracts the results
# against the yield weights, and writes the tables and the provenance under
# data/sims/<system>/<run>, the DrWatson layout. Nothing is plotted here — plotting/figures.jl
# reads what this writes, so a figure can be redrawn without solving the model again.
#
# A run directory is never overwritten: one that already exists is moved aside as <run>#1,
# <run>#2, … before the new one is written. Set DRWATSON_STOREPATCH=true to record the
# uncommitted diff of a dirty checkout in the metadata beside the commit.

include(joinpath(@__DIR__, "..", "activate.jl"))

using DeterministicSequentialEmission
using DrWatson: datadir, projectdir, tag!
using Dates: Dates
using InteractiveUtils: versioninfo
using LinearAlgebra: BLAS
using Printf: @printf
using TOML: TOML

function main(arguments::Vector{String})
    if !(1 <= length(arguments) <= 2)
        println("usage: julia scripts/run.jl <configuration.toml> [data directory]")
        return 1
    end
    path = isabspath(arguments[1]) ? arguments[1] : projectdir(arguments[1])
    data = length(arguments) == 2 ? arguments[2] : datadir()

    configuration = load_configuration(path; data_directory = data)
    println(configuration)
    result = run_pipeline(configuration)

    destination =
        datadir("sims", system_label(configuration.system), run_identifier(configuration))
    move_aside(destination)
    write_results(result, destination)
    write_provenance(destination, result, path, data)

    @printf(
        "  %d of %d configurations solved (%.1f %%)\n",
        result.configurations_solved,
        result.configurations_attempted,
        100 * result.configurations_solved / result.configurations_attempted
    )
    @printf(
        "  the segmented curve covers %.0f %% of the heavy masses\n",
        100 * result.segment_coverage
    )
    @printf(
        "  the kinetic energy grid reaches %.1f %% of the yield\n",
        100 * result.yield_reached
    )
    @printf("  ⟨ν⟩ = %.4f per fission\n", 2 * result.total_multiplicity)
    @printf("  ⟨TKE⟩ = %.3f MeV\n", result.mean_kinetic_energy)
    isnan(result.mean_neutron_lab_energy) ||
        @printf("  ⟨E⟩ = %.4f MeV in the laboratory\n", result.mean_neutron_lab_energy)
    println("  written to ", destination)
    return 0
end

# Never overwrite a run. An existing directory moves to <name>#1, <name>#2, …: the numbering
# DrWatson's safesave applies to a file, applied here to a whole run directory.
function move_aside(directory::AbstractString)
    isdir(directory) || return nothing
    index = 1
    while ispath("$(directory)#$(index)")
        index += 1
    end
    mv(directory, "$(directory)#$(index)")
    return nothing
end

# The provenance record: what the library knows about the run, the identifier and the paths,
# the commit through DrWatson's tag!, and the platform, written as metadata.toml; beside it,
# copies of the configuration and of the resolved Manifest. A result then traces to
# configuration, commit, environment and hardware from its own directory.
function write_provenance(
    directory::AbstractString,
    result::RunResult,
    configuration_path::AbstractString,
    data_directory::AbstractString,
)
    metadata = run_metadata(result)
    run = Dict{String, Any}(
        "identifier" => run_identifier(result.configuration),
        "configuration" => relpath(configuration_path, projectdir()),
        "data_directory" => relpath(data_directory, projectdir()),
        "script" => relpath(@__FILE__, projectdir()),
        "written" => string(Dates.now()),
        "package_version" => string(pkgversion(DeterministicSequentialEmission)),
    )
    tag!(run; gitpath = projectdir())
    metadata["run"] = run
    processors = Sys.cpu_info()
    metadata["platform"] = Dict{String, Any}(
        "julia" => string(VERSION),
        "hostname" => gethostname(),
        "cpu" => isempty(processors) ? "unknown" : String(first(processors).model),
        "cpu_threads" => Sys.CPU_THREADS,
        "threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "total_memory_gb" => round(Sys.total_memory() / 2^30; digits = 2),
        "versioninfo" => sprint(versioninfo),
    )
    open(joinpath(directory, "metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted = true)
    end
    cp(configuration_path, joinpath(directory, "configuration.toml"))
    manifest = projectdir("Manifest.toml")
    isfile(manifest) && cp(manifest, joinpath(directory, "Manifest.toml"))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
