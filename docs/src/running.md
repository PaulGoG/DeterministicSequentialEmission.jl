```@meta
CurrentModule = DeterministicSequentialEmission
```

# Running a case

One TOML file states everything a run needs, and is validated in full before anything is
solved.

```
julia scripts/run.jl config/Cf252_sf.toml
```

`config/Cf252_sf.toml` is the annotated reference. Each key carries its allowed values or its
bounds as a comment, and [`load_configuration`](@ref) enforces exactly what the comment says,
naming the offending key when it cannot.

[`run_pipeline`](@ref) sweeps the fragmentation domain against the kinetic energy grid, solves
the emission chain at every configuration the physics admits, contracts the result against the
yield weights, and returns a [`RunResult`](@ref). A configuration is dropped only where the
physics excludes it — a non-positive `TXE`, a partition the segments do not cover, a fragment
the mass table cannot reach — and the counts of attempted and solved configurations are written
alongside the results, so the coverage of a run can be audited from its own output.

`scripts/run.jl` writes the result on the DrWatson layout. The tables go to
`data/sims/<system>/<run>` through [`write_results`](@ref), where the system is
[`system_label`](@ref) and the run is [`run_identifier`](@ref) — `savename` over
[`run_parameters`](@ref) — so two runs of one system that differ in anything the identifier
covers land in different directories, and a directory name states what produced it. Beside the
tables the script writes `metadata.toml` — [`run_metadata`](@ref) from the library, plus the
identifier, the paths, the commit through DrWatson's `tag!`, and the platform — and copies of
the configuration and the resolved `Manifest.toml`, so a result traces to configuration,
commit, environment and hardware from its own directory. A directory that already holds a run
is moved aside as `<run>#1`, `<run>#2`, … rather than overwritten.

```@autodocs
Modules = [DeterministicSequentialEmission]
Pages = ["configuration.jl", "pipeline.jl"]
```

