include(joinpath(@__DIR__, "activate.jl"))

using DeterministicSequentialEmission
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :numeric)

DocMeta.setdocmeta!(
    DeterministicSequentialEmission,
    :DocTestSetup,
    :(using DeterministicSequentialEmission);
    recursive = true,
)

makedocs(;
    modules = [DeterministicSequentialEmission],
    authors = "Paul-Adrian Gogîță",
    sitename = "DeterministicSequentialEmission.jl",
    repo = Documenter.Remotes.GitHub("PaulGoG", "DeterministicSequentialEmission.jl"),
    format = Documenter.HTML(;
        canonical = "https://PaulGoG.github.io/DeterministicSequentialEmission.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Running a case" => "running.md",
        "Naming" => "naming.md",
        "System and inputs" => "system.md",
        "Fragmentation and charge" => "fragmentation.md",
        "Energetics and partition" => "energetics.md",
        "Emission and spectrum" => "emission.md",
        "Yields and observables" => "observables.md",
        "References" => "references.md",
    ],
    plugins = [bib],
)

deploydocs(;
    repo = "github.com/PaulGoG/DeterministicSequentialEmission.jl",
    devbranch = "main",
)
