@testset "neutron strength function" begin
    # Piecewise linear with plateaux, continuous at the joints the parameterization declares.
    for A in (23, 55, 69, 79, 90, 110, 120, 135, 146, 152, 158, 173, 186, 208)
        left = neutron_strength_function(A)
        right = neutron_strength_function(A + 1)
        @test isfinite(left) && isfinite(right)
        @test left > 0 && right > 0
    end
    @test neutron_strength_function(10) == 2.5e-5
    @test neutron_strength_function(300) == 1e-4

    # The 3s and 4s giant resonances show as maxima near A ≈ 50 and A ≈ 150.
    @test neutron_strength_function(50) > neutron_strength_function(100)
    @test neutron_strength_function(150) > neutron_strength_function(110)
end

@testset "evaporation moments" begin
    constant = ConstantCrossSection()
    variable = VariableCrossSection()

    @test mean_neutron_energy(constant, 1.0) == 2.0
    @test mean_neutron_energy(constant, 0.5) == 1.0

    # As α → 0 the variable cross section reduces to the constant one.
    @test mean_neutron_energy(variable, 1.0e-12, 1.0) ≈ 2.0 atol = 1e-6

    # A non-zero α lowers the mean energy, the low-energy enhancement of σ_c.
    @test mean_neutron_energy(variable, 0.3, 1.0) < 2.0

    @test fermi_gas_energy(13.0, 0.5) ≈ 13.0 * 0.25
end

@testset "residual temperature" begin
    constant = ConstantCrossSection()
    variable = VariableCrossSection()

    # The closed-form root satisfies its own equation.
    a = 13.0
    available = 10.0
    T = residual_temperature(constant, a, available)
    @test a * T^2 + 2T ≈ available

    # The bisection root satisfies the transcendental equation it was found for.
    α = alpha_parameter(133)
    T = residual_temperature(variable, a, available; α = α)
    @test a * T^2 + mean_neutron_energy(variable, α, T) ≈ available atol = 1e-8

    # Below the Fermi-gas floor there is no nuclear temperature to report.
    @test residual_temperature(constant, a, 1e-9) === nothing
    @test residual_temperature(variable, a, 1e-9; α = α) === nothing
    @test residual_temperature(constant, a, -1.0) === nothing

    # The optional ceiling rejects rather than truncates, and is off by default.
    hot = residual_temperature(constant, a, 100.0)
    @test hot > 2.0
    @test residual_temperature(constant, a, 100.0; maximum_temperature = 2.0) === nothing
end

@testset "sequential emission against the archived run" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    cross_section = VariableCrossSection()

    reference = read_table(
        joinpath(REFERENCE, "CF52_sequences.dat"),
        TableSpec([:A_H, :Z_H, :TKE, :side, :k, :S_n, :a, :T, :E_r, :eps]; skip = 1),
    )
    masses = column(reference, :A_H, Int)
    charges = column(reference, :Z_H, Int)
    energies = column(reference, :TKE, Float64)
    sides = column(reference, :side)
    indices = column(reference, :k, Int)
    expected_S = column(reference, :S_n, Float64)
    expected_a = column(reference, :a, Float64)
    expected_T = column(reference, :T, Float64)
    expected_E = column(reference, :E_r, Float64)
    expected_ε = column(reference, :eps, Float64)

    method = ExcitationEnergyRatio(
        JointedSegments([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )

    chains = Dict{Tuple{Int, Int, Float64, String}, EmissionChain}()
    for row in eachindex(masses)
        key = (masses[row], charges[row], energies[row], sides[row])
        haskey(chains, key) && continue
        fragmentation = fragmentation_at(masses[row], charges[row])
        balance = energetics(table, CF252, fragmentation, energies[row])
        split = partition(method, fragmentation, value(balance.TXE))
        fragment = sides[row] == "L" ? fragmentation.light : fragmentation.heavy
        excitation = sides[row] == "L" ? split.light : split.heavy
        # The archived run was produced with the ceiling on the residual temperature active;
        # see the assertion below, which fails without it.
        chains[key] = emission_chain(
            table,
            model,
            cross_section,
            fragment,
            excitation;
            maximum_temperature = 2.0,
        )
    end

    worst_S = 0.0
    worst_a = 0.0
    worst_T = 0.0
    worst_E = 0.0
    worst_ε = 0.0
    compared = 0
    for row in eachindex(masses)
        chain = chains[(masses[row], charges[row], energies[row], sides[row])]
        @test indices[row] <= length(chain)
        indices[row] <= length(chain) || continue
        sequence = chain.sequences[indices[row]]
        @test sequence.index == indices[row]
        worst_S = max(worst_S, abs(sequence.separation_energy - expected_S[row]))

        # A sequence the archive wrote with a zero temperature is one the equation could not
        # solve; the rewrite reports it as absent rather than as a zero that averages in.
        if expected_T[row] > 0
            @test is_solved(sequence)
            is_solved(sequence) || continue
            worst_a = max(worst_a, abs(sequence.level_density_parameter - expected_a[row]))
            worst_T = max(worst_T, abs(sequence.temperature - expected_T[row]))
            worst_E = max(worst_E, abs(sequence.residual_energy - expected_E[row]))
            worst_ε = max(worst_ε, abs(sequence.mean_neutron_energy - expected_ε[row]))
            compared += 1
        else
            @test !is_solved(sequence)
        end
    end

    @test compared > 150
    @test worst_S < 1e-8
    @test worst_a < 1e-8
    # The archived temperatures are rounded to eight decimals; the root itself is found by
    # bisection to a tighter tolerance than that.
    @test worst_T < 1e-7
    @test worst_E < 1e-6
    @test worst_ε < 1e-7
end

@testset "the archived run used the 2 MeV ceiling" begin
    # The prototype source carries the ceiling commented out, but the archived run predates
    # that edit: there are configurations whose first residual temperature exceeds 2 MeV, and
    # the archive records the whole chain as unsolved. Without the ceiling the same chain
    # solves, so the reference itself distinguishes the two settings.
    table = mass_table()
    model = BackShiftedFermiGas(table)
    cross_section = VariableCrossSection()
    method = ExcitationEnergyRatio(
        JointedSegments([(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]),
    )

    # A_H = 134, Z_H = 52 at TKE = 150 MeV: the archive writes zeros for the light fragment.
    fragmentation = fragmentation_at(134, 52)
    balance = energetics(table, CF252, fragmentation, 150.0)
    split = partition(method, fragmentation, value(balance.TXE))

    capped = emission_chain(
        table,
        model,
        cross_section,
        fragmentation.light,
        split.light;
        maximum_temperature = 2.0,
    )
    uncapped = emission_chain(table, model, cross_section, fragmentation.light, split.light)

    @test capped.solved == 0
    @test uncapped.solved > 0
    @test uncapped.sequences[1].temperature > 2.0

    # The ceiling costs the whole chain, not just the sequence that tripped it: the first
    # failure ends the solving stage. That is the information loss it was flagged for.
    @test capped.possible == uncapped.possible
    @test sequence_average(capped, :temperature) === nothing
    @test sequence_average(uncapped, :temperature) !== nothing
end

@testset "chain bookkeeping" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    chain = emission_chain(table, model, VariableCrossSection(), Nuclide(52, 134), 20.0)

    # n counts every prompt neutron; k_max only the sequences with a temperature.
    @test chain.possible >= chain.solved
    @test chain.possible == length(chain.sequences)
    @test chain.solved == count(is_solved, chain.sequences)

    # A_p = A − n is an integer by construction, which is the point of the treatment.
    @test post_neutron_mass(chain) == 134 - chain.possible
    @test post_neutron_mass(chain) isa Int

    # ν is a different quantity from n and must not be confused with it.
    @test neutron_multiplicity(chain) == (chain.possible + 1) / 2

    # Sequence averages run over k_max, not over n.
    average = sequence_average(chain, :temperature)
    @test average ≈
          sum(s.temperature for s in chain.sequences if is_solved(s)) / chain.solved

    @test_throws ArgumentError sequence_average(chain, :nonsense)

    # A fragment with too little excitation to pay a separation energy emits nothing.
    cold = emission_chain(table, model, VariableCrossSection(), Nuclide(52, 134), 0.5)
    @test cold.possible == 0
    @test neutron_multiplicity(cold) == 0.0
    @test post_neutron_mass(cold) == 134
    @test sequence_average(cold, :temperature) === nothing
end

@testset "energy is conserved along a chain" begin
    table = mass_table()
    model = BackShiftedFermiGas(table)
    cross_section = VariableCrossSection()
    chain = emission_chain(table, model, cross_section, Nuclide(52, 134), 25.0)

    # Each solved sequence spends exactly Sₙ + ⟨ε⟩ of the residual energy.
    previous = 25.0
    for sequence in chain.sequences
        is_solved(sequence) || break
        spent = sequence.separation_energy + sequence.mean_neutron_energy
        @test sequence.residual_energy ≈ previous - spent atol = 1e-9
        previous = sequence.residual_energy
    end

    # Temperatures fall monotonically along the chain, as the residual cools.
    solved = [s.temperature for s in chain.sequences if is_solved(s)]
    @test length(solved) >= 2
    @test issorted(solved; rev = true)
end
