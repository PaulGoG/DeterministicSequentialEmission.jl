```@meta
CurrentModule = DeterministicSequentialEmission
```

# Running a case

One TOML file states everything a run needs, and is validated in full before anything is
solved.

```
julia --project scripts/run.jl config/Cf252_0f.toml
```

`config/Cf252_0f.toml` is the annotated reference. Each key carries its allowed values or its
bounds as a comment, and [`load_configuration`](@ref) enforces exactly what the comment says,
naming the offending key when it cannot.

[`run_pipeline`](@ref) sweeps the fragmentation domain against the kinetic energy grid, solves
the emission chain at every configuration the physics admits, contracts the result against the
yield weights, and hands a [`RunResult`](@ref) to [`write_results`](@ref). A configuration is
dropped only where the physics excludes it — a non-positive `TXE`, a partition the segments do
not cover, a fragment the mass table cannot reach — and the counts of attempted and solved
configurations are written alongside the results, so the coverage of a run can be audited from
its own output.

```@autodocs
Modules = [DeterministicSequentialEmission]
Pages = ["configuration.jl", "pipeline.jl"]
```

