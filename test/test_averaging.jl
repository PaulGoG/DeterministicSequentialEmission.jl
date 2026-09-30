guarded(
    CHARGE_DISTRIBUTION_AVAILABLE && YIELD_AVAILABLE,
    "averaged quantities against the archived run",
) do
    @testset "averaged quantities against the archived run" begin
        table = mass_table()
        distribution = cf252_charge_distribution()
        domain = cf252_domain()
        model = BackShiftedFermiGas(table)
        cross_section = VariableCrossSection()
        method = ExcitationEnergyRatio(
            SegmentedCurve([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
        )
        yields = build_fragment_yield(
            domain,
            distribution,
            read_mass_energy_yield(CF252_YIELD_FILE),
        )
        weights = weights_from(yields)

        ν = FragmentQuantity{Float64}()
        T = FragmentQuantity{Float64}()
        excitation = FragmentQuantity{Float64}()
        for entry in domain.entries, TKE in CF252_TKE
            balance = energetics(table, CF252, entry, TKE)
            balance === nothing && continue
            split = partition(method, entry, value(balance.TXE))
            split === nothing && continue
            for (fragment, E) in ((entry.light, split.light), (entry.heavy, split.heavy))
                chain = emission_chain(
                    table,
                    model,
                    cross_section,
                    fragment,
                    E;
                    maximum_temperature = 2.0,
                )
                key = (fragment.A, fragment.Z, TKE)
                ν[key] = neutron_multiplicity(chain)
                excitation[key] = E
                t = sequence_average(chain, :temperature)
                t === nothing || (T[key] = t)
            end
        end

        for (quantity, file, columns, over) in (
            (ν, "nu_vs_A.dat", [:A, :nu], (:charge, :kinetic_energy)),
            (T, "T_vs_A.dat", [:A, :T], (:charge, :kinetic_energy)),
            (excitation, "E_star_vs_A.dat", [:A, :E_star], (:charge, :kinetic_energy)),
            (excitation, "E_star_vs_TKE.dat", [:TKE, :E_star], (:mass, :charge)),
        )
            reference = read_delimited_table(
                joinpath(REFERENCE, file),
                TableSpec(columns; skip = 1),
            )
            arguments = column(reference, columns[1], Float64)
            expected = column(reference, columns[2], Float64)
            got = Dict(average(quantity, weights; over = over))

            @test length(got) == length(reference)
            @test issubset(arguments, keys(got))
            # Reproduced in double precision, not to a tolerance chosen to pass.
            @test maximum(
                abs(get(got, x, NaN) - y) for (x, y) in zip(arguments, expected)
            ) < 1e-11
        end

        # The total averages of the archive's own summary file. The light group runs to A₀/2
        # inclusive and so does the heavy one: both fragments of a symmetric split carry that
        # mass, which is why the two groups sum to slightly more than twice the whole.
        @test total_average(ν, weights; masses = 78:126) ≈ 1.7790091580039273 rtol = 1e-9
        @test total_average(ν, weights; masses = 126:174) ≈ 1.7010552393293843 rtol = 1e-9
        @test total_average(ν, weights) ≈ 1.7398212275804683 rtol = 1e-9
    end
end
