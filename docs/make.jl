include(joinpath(@__DIR__, "activate.jl"))

using DeterministicSequentialEmission
using FissionFragmentsDomain
using Documenter
using DocumenterCitations
using DocumenterInterLinks

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :numeric)

# Symbols of the domain package resolve against its deployed documentation.
links = InterLinks(
    "FissionFragmentsDomain" => "https://PaulGoG.github.io/FissionFragmentsDomain.jl/stable/",
)
fallbacks = ExternalFallbacks(; automatic = true)

DocMeta.setdocmeta!(
    DeterministicSequentialEmission,
    :DocTestSetup,
    :(using DeterministicSequentialEmission, FissionFragmentsDomain);
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
        "The fragmentation domain" => "domain.md",
        "Partition" => "energetics.md",
        "Emission and spectrum" => "emission.md",
        "References" => "references.md",
    ],
    plugins = [bib, links, fallbacks],
)

deploydocs(;
    repo = "github.com/PaulGoG/DeterministicSequentialEmission.jl",
    devbranch = "main",
)
