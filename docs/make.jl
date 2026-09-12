using DeterministicSequentialEmission
using Documenter

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
    # Stated explicitly: the repository is local-only until release, so Documenter cannot
    # infer the remote from a configured origin.
    repo = Documenter.Remotes.GitHub("PaulGoG", "DeterministicSequentialEmission.jl"),
    format = Documenter.HTML(;
        canonical = "https://PaulGoG.github.io/DeterministicSequentialEmission.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/PaulGoG/DeterministicSequentialEmission.jl",
    devbranch = "main",
)
