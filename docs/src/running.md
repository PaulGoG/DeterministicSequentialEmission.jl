```@meta
CurrentModule = DeterministicSequentialEmission
```

# Running a case

Simulation runs are fully defined by a single TOML configuration file, which is validated before
computation begins.

```
julia scripts/run.jl config/Cf252_sf.toml
```

`config/Cf252_sf.toml` serves as the annotated reference. Configuration keys declare allowed value sets
or physical bounds, and [`load_configuration`](@ref) verifies these constraints on load, reporting specific
validation errors if an invalid parameter is encountered.

[`run_pipeline`](@ref) sweeps the fragmentation domain across the kinetic energy grid, evaluates the
emission cascade at each physical configuration, contracts the observables against the fragment yield
distribution, and returns a [`RunResult`](@ref). Configurations are excluded only when unphysical:
negative `TXE`, excitation partitions exceeding the validity domain, or nuclei outside the mass table.
Diagnostic counts of attempted and solved configurations are recorded alongside the results, enabling
direct inspection of domain coverage.

`scripts/run.jl` organizes output tables following the DrWatson convention. Results are written to
`data/sims/<system>/<run>` via [`write_results`](@ref), where the system is [`system_label`](@ref) and
the run directory is designated by [`run_identifier`](@ref) (constructed via `savename` over
[`run_parameters`](@ref)). Consequently, calculations differing in any functional parameter are placed
in distinct directories. Along with data tables, each run directory preserves `metadata.toml` as well
as exact copies of the configuration file and resolved `Manifest.toml`. `metadata.toml` records
[`run_metadata`](@ref), platform details, input paths, and the exact git commit via DrWatson's `tag!`.
Existing run directories are preserved by renaming previous calculations with sequential numerical suffixes
(`<run>#1`, `<run>#2`, ...) rather than overwriting.

```@autodocs
Modules = [DeterministicSequentialEmission]
Pages = ["configuration.jl", "pipeline.jl"]
```

