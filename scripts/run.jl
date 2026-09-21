# Entry point: julia scripts/run.jl config/<case>.toml [data directory]
#
# Solves every configuration the domain and the kinetic energy grid span, contracts the results
# against the yield weights, and writes the tables under results/<system>/<run>. Nothing is
# plotted here — plotting/figures.jl reads what this writes, so a figure can be redrawn without
# solving the model again.

include(joinpath(@__DIR__, "..", "activate.jl"))

using DeterministicSequentialEmission
using Printf

function main(arguments::Vector{String})
    if !(1 <= length(arguments) <= 2)
        println("usage: julia scripts/run.jl <configuration.toml> [data directory]")
        return 1
    end
    root = dirname(@__DIR__)
    path = isabspath(arguments[1]) ? arguments[1] : joinpath(root, arguments[1])
    data = length(arguments) == 2 ? arguments[2] : joinpath(root, "data")

    configuration = load_configuration(path; data_directory = data)
    println(configuration)
    result = run_pipeline(configuration; write = false)
    destination = write_results(result)

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

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
