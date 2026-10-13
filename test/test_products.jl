# The post-neutron product distributions: the product count on chains built by hand, the
# published count against the archived ²⁵²Cf runs, the invariants of the contraction on the same
# sweeps, and the configuration section.

"A chain of `possible` sequences from `fragment`, the first `solved` of them with a temperature."
function hand_built_chain(fragment::Nuclide, possible::Integer, solved::Integer)
    sequences = EmissionSequence[
        k <= solved ? EmissionSequence(k, 6.0, 13.0, 1.0 - 0.2k, 10.0 - 2k, 1.5) :
        EmissionSequence(k, 6.0, 13.0, nothing, nothing, nothing) for k in 1:possible
    ]
    return EmissionChain(fragment, sequences)
end

@testset "product masses" begin
    fragment = Nuclide(50, 132)

    # A chain that emits nothing leaves the fragment as it is, under either count.
    cold = hand_built_chain(fragment, 0, 0)
    @test product_masses(cold, EmissionProbabilityCount()) == [(132, 1.0)]
    @test product_masses(cold, SequenceCount()) == [(132, 1.0)]

    single = hand_built_chain(fragment, 1, 1)
    @test product_masses(single, EmissionProbabilityCount()) == [(131, 1.0)]
    @test product_masses(single, SequenceCount()) == [(131, 1.0)]

    # Three possible sequences, one of them without a temperature: n counts all three.
    hot = hand_built_chain(fragment, 3, 2)
    @test hot.possible == 3
    @test hot.solved == 2
    spread = product_masses(hot, EmissionProbabilityCount())
    @test first.(spread) == [131, 130, 129]
    @test all(p -> p ≈ 1 / 3, last.(spread))
    @test sum(last, spread) ≈ 1 rtol = RTOL
    @test product_masses(hot, SequenceCount()) == [(129, 1.0)]
    @test product_masses(hot, SequenceCount()) == [(post_neutron_mass(hot), 1.0)]

    # The mean mass loss under the emission probability count is the multiplicity.
    for chain in (cold, single, hot)
        loss = sum(
            p * (fragment.A - A_p) for
            (A_p, p) in product_masses(chain, EmissionProbabilityCount())
        )
        @test loss ≈ neutron_multiplicity(chain) atol = 1e-15
    end
end

guarded(
    CHARGE_DISTRIBUTION_AVAILABLE && YIELD_AVAILABLE,
    "post-neutron products against the archived runs",
) do
    @testset "post-neutron products against the archived runs" begin
        table = mass_table()
        distribution = cf252_charge_distribution()
        domain = cf252_domain()
        model = BackShiftedFermiGas(table)
        cross_section = VariableCrossSection()
        yields = build_fragment_yield(
            domain,
            distribution,
            read_mass_energy_yield(CF252_YIELD_FILE),
        )
        weights = weights_from(yields)
        A₀ = CF252.compound.A

        for (method, file) in (
            (
                ExcitationEnergyRatio(
                    SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
                ),
                "Y_post_vs_A_p_excitation_ratio.dat",
            ),
            (
                TemperatureRatio(
                    SegmentedCurve([
                        (126, 1.0),
                        (130, 1.734),
                        (134, 1.28),
                        (138, 1.28),
                        (150, 0.9),
                        (167, 0.68),
                        (174, 0.59),
                    ]),
                    model,
                ),
                "Y_post_vs_A_p_temperature_ratio.dat",
            ),
        )
            # The sweep of run_pipeline, on the archived grid.
            chains = Tuple{
                Fragmentation,
                Float64,
                EmissionChain,
                EmissionChain,
                ExcitationEnergies,
            }[]
            ν = FragmentQuantity{Float64}()
            for entry in domain.entries, TKE in CF252_TKE
                balance = energetics(table, CF252, entry, TKE)
                balance === nothing && continue
                split = partition(method, entry, value(balance.TXE))
                split === nothing && continue
                light, heavy = (
                    emission_chain(
                        table,
                        model,
                        cross_section,
                        fragment,
                        E;
                        maximum_temperature = 2.0,
                    ) for (fragment, E) in
                    ((entry.light, split.light), (entry.heavy, split.heavy))
                )
                ν[(entry.light.A, entry.light.Z, TKE)] = neutron_multiplicity(light)
                ν[(entry.heavy.A, entry.heavy.Z, TKE)] = neutron_multiplicity(heavy)
                push!(chains, (entry, TKE, light, heavy, split))
            end
            contract(count) = DeterministicSequentialEmission._product_distributions(
                chains,
                weights,
                A₀,
                count,
                yields.total,
            )

            # The published count, A_p = A − n, reproduces the archived Y_p(A_p).
            published = contract(SequenceCount())
            reference = read_reference(file, [:A_p, :Y], 1)
            got = Dict(Int(A_p) => Y for (A_p, Y) in published.mass_yield)
            @test Set(keys(got)) == Set(A_p for (A_p,) in keys(reference))
            # Relative where the yield is resolved, absolute in the far tail.
            relative = maximum(
                abs(get(got, A_p, NaN) - Y) / Y for ((A_p,), Y) in reference if Y > 1e-10
            )
            absolute = maximum(
                (abs(get(got, A_p, NaN) - Y) for ((A_p,), Y) in reference if Y <= 1e-10);
                init = 0.0,
            )
            # A constant factor would mean a normalisation, not a count, differs; say so.
            relative < 1e-8 || @error "Y_p(A_p) departs from $file" relative extrema(
                get(got, A_p, NaN) / Y for ((A_p,), Y) in reference if Y > 1e-3
            )
            @test relative < 1e-8
            @test absolute < 1e-10

            # The invariants, under either count.
            solved = Set(
                (fragment.A, fragment.Z, TKE) for (entry, TKE, _, _, _) in chains for
                fragment in (entry.light, entry.heavy)
            )
            solved_weight = sum(get(weights, key, 0.0) for key in solved)
            max_KE = maximum(TKE * (A₀ - A) / A₀ for (A, _, TKE) in solved)
            for count in (EmissionProbabilityCount(), SequenceCount())
                products = count isa SequenceCount ? published : contract(count)

                # The products of the solved configurations carry the whole yield.
                @test sum(last, products.mass_yield) ≈ yields.total atol = 1e-10

                # Y_p(A_p) is Y_p(A_p, Z) summed over Z, and Y_p(N_p) the same yields regrouped by
                # N_p = A_p − Z.
                by_mass = Dict{Float64, Float64}()
                by_neutron = Dict{Float64, Float64}()
                by_charge = Dict{Int, Float64}()
                for (A_p, Z, Y) in products.yield
                    by_mass[A_p] = get(by_mass, A_p, 0.0) + Y
                    by_neutron[A_p - Z] = get(by_neutron, A_p - Z, 0.0) + Y
                    by_charge[Z] = get(by_charge, Z, 0.0) + Y
                end
                @test all(
                    isapprox(by_mass[A_p], Y; atol = 1e-12) for
                    (A_p, Y) in products.mass_yield
                )
                @test all(
                    isapprox(by_neutron[N_p], Y; atol = 1e-12) for
                    (N_p, Y) in products.neutron_yield
                )
                @test length(by_neutron) == length(products.neutron_yield)

                # Emission leaves the charge alone: Σ_{A_p} Y_p(A_p, Z) is the input Y(Z) over the
                # solved configurations, renormalised as the products are,
                # (total / Σ_solved Y) Σ_{(A, TKE) solved} Y(A, Z, TKE).
                input_by_charge = Dict{Int, Float64}()
                for (A, Z, TKE) in solved
                    input_by_charge[Z] =
                        get(input_by_charge, Z, 0.0) + get(weights, (A, Z, TKE), 0.0)
                end
                @test Set(keys(by_charge)) == Set(keys(input_by_charge))
                @test all(
                    isapprox(
                        by_charge[Z],
                        yields.total * Y / solved_weight;
                        rtol = 1e-10,
                        atol = 1e-12,
                    ) for (Z, Y) in input_by_charge
                )

                # Neutrons carry kinetic energy away, never bring it.
                @test all(KE <= max_KE for (_, KE) in products.kinetic_energy)
                @test all(KE <= max_KE for (_, KE) in products.kinetic_energy_by_charge)
                @test all(
                    TKE <= last(CF252_TKE) for (_, TKE) in products.total_kinetic_energy
                )
            end

            # Under the default count the products lose exactly the multiplicity; under the
            # published one they lose ⟨n⟩, about 2ν − 1.
            @test contract(EmissionProbabilityCount()).mass_loss ≈ total_average(ν, weights) atol =
                1e-12
            @test published.mass_loss > total_average(ν, weights) + 0.5
        end
    end
end

@testset "[products] configuration" begin
    mktempdir() do directory
        write(
            joinpath(directory, "SEGMENTS"),
            "A_H,excitation_ratio,excitation_ratio_uncertainty\n126,0.5,0.0\n150,0.607,0.0\n\
             174,0.865,0.0\n",
        )
        isfile(MASS_EXCESS_FILE) || return nothing
        cp(MASS_EXCESS_FILE, joinpath(directory, "MASSES"))

        # A configuration without the section loads with the defaults.
        default = _configuration(directory)
        @test default.products_enabled
        @test default.product_count isa EmissionProbabilityCount
        @test run_parameters(default)["pc"] == "emission_probability"

        published = _configuration(
            directory,
            "[output]" => "[products]\ncount = \"sequence_count\"\n\n[output]",
        )
        @test published.product_count isa SequenceCount
        @test run_parameters(published)["pc"] == "sequence_count"
        @test occursin("pc=sequence_count", run_identifier(published))
        @test run_identifier(published) != run_identifier(default)

        disabled = _configuration(
            directory,
            "[output]" => "[products]\nenabled = false\n\n[output]",
        )
        @test !disabled.products_enabled
        # The count names the run whether or not the products are written, as the sequence
        # weighting does whether or not the spectrum is.
        @test run_parameters(disabled)["pc"] == "emission_probability"

        @test _rejects(
            directory,
            "products.count",
            "[output]" => "[products]\ncount = \"mean\"\n\n[output]",
        )
        @test _rejects(
            directory,
            "products.enabled",
            "[output]" => "[products]\nenabled = \"yes\"\n\n[output]",
        )
        @test _rejects(
            directory,
            "products.weighting",
            "[output]" => "[products]\nweighting = \"uniform\"\n\n[output]",
        )
    end
end
