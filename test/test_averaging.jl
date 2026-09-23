@testset "contractions of a rank-3 quantity" begin
    # Synthetic, so the contraction itself is exercised without any input data.
    q = FragmentQuantity{Float64}()
    w = FragmentQuantity{Float64}()
    for A in (100, 140), Z in (40, 41), TKE in (150.0, 170.0)
        q[(A, Z, TKE)] = A + Z / 100 + TKE / 1000
        w[(A, Z, TKE)] = Z == 40 ? 3.0 : 1.0
    end

    # Reducing over two indices leaves the third, ascending.
    by_mass = @inferred average(q, w; over = (:charge, :kinetic_energy))
    @test first.(by_mass) == [100.0, 140.0]
    @test length(by_mass) == 2

    # The weighted mean, worked by hand at A = 100: Z = 40 carries three times the weight.
    expected =
        (
            3 * (100.40 + 0.150) +
            1 * (100.41 + 0.150) +
            3 * (100.40 + 0.170) +
            1 * (100.41 + 0.170)
        ) / 8
    @test last(by_mass[1]) ≈ expected rtol = RTOL

    by_charge = average(q, w; over = (:mass, :kinetic_energy))
    @test first.(by_charge) == [40.0, 41.0]
    by_energy = average(q, w; over = (:mass, :charge))
    @test first.(by_energy) == [150.0, 170.0]

    # Reducing one index leaves two, which this signature does not return.
    @test_throws ArgumentError average(q, w; over = (:mass, :charge, :kinetic_energy))
    @test_throws ArgumentError average(q, w; over = (:mass, :mass))
    @test_throws ArgumentError average(q, w; over = (:temperature,))

    # The total is the same contraction carried to a scalar, and restricting the masses is how
    # the light and heavy groups are separated.
    whole = total_average(q, w)
    @test whole ≈ (sum(w[k] * q[k] for k in keys(q)) / sum(w[k] for k in keys(q))) rtol =
        RTOL
    @test total_average(q, w; masses = 100:100) ≈ last(by_mass[1]) rtol = RTOL
    @test total_average(q, w; masses = 200:300) === nothing
end

@testset "a configuration enters both sums or neither" begin
    # The quantity is tabulated on a coarser grid than the weight, which is the archived
    # configuration's own situation. A cell the model never solved must not dilute the mean by
    # entering the denominator alone.
    q = FragmentQuantity{Float64}()
    w = FragmentQuantity{Float64}()
    q[(140, 54, 150.0)] = 2.0
    q[(140, 54, 170.0)] = 4.0
    for TKE in (150.0, 160.0, 170.0)
        w[(140, 54, TKE)] = 1.0
    end

    # Averaging over the two cells that exist gives 3.0, not 2.0 as a diluted denominator would.
    @test last(only(average(q, w; over = (:charge, :kinetic_energy)))) ≈ 3.0 rtol = RTOL
    @test total_average(q, w) ≈ 3.0 rtol = RTOL

    # And a weight with no quantity contributes nothing at all.
    @test length(average(q, w; over = (:charge, :kinetic_energy))) == 1
end

@testset "uncertainty through a weighted mean" begin
    # Both terms of the propagation: the quantity's own uncertainty and the weights'.
    q = FragmentQuantity{Measurement{Float64}}()
    w = FragmentQuantity{Measurement{Float64}}()
    q[(140, 54, 150.0)] = measurement(2.0, 0.1)
    q[(140, 54, 170.0)] = measurement(4.0, 0.2)
    w[(140, 54, 150.0)] = measurement(1.0, 0.05)
    w[(140, 54, 170.0)] = measurement(3.0, 0.05)

    mean = total_average(q, w)
    @test value(mean) ≈ 3.5 rtol = RTOL

    # Worked from the standard expression, both terms:
    #   δ² = Σ (Wᵢ/ΣW · δqᵢ)² + Σ ((qᵢ − ⟨q⟩)/ΣW · δWᵢ)²
    from_quantity = (1 / 4 * 0.1)^2 + (3 / 4 * 0.2)^2
    from_weights = ((2.0 - 3.5) / 4 * 0.05)^2 + ((4.0 - 3.5) / 4 * 0.05)^2
    @test uncertainty(mean) ≈ sqrt(from_quantity + from_weights) rtol = 1e-10

    # A quantity with no uncertainty leaves only the weight term, and a bare Float64 quantity
    # against Measurement weights still propagates the weight term.
    exact = FragmentQuantity{Float64}()
    exact[(140, 54, 150.0)] = 2.0
    exact[(140, 54, 170.0)] = 4.0
    only_weights = total_average(exact, w)
    @test uncertainty(only_weights) ≈ sqrt(from_weights) rtol = 1e-10

    # An absent uncertainty propagates as absent rather than as zero, which is how the
    # prototype marks a yield that arrived without one.
    unmeasured = FragmentQuantity{Measurement{Float64}}()
    unmeasured[(140, 54, 150.0)] = measurement(1.0, NaN)
    unmeasured[(140, 54, 170.0)] = measurement(3.0, 0.05)
    blurred = total_average(q, unmeasured)
    @test isfinite(value(blurred))
    @test isnan(uncertainty(blurred))
end

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
