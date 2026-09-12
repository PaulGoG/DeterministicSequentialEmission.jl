@testset "charge polarization" begin
    polarization = cf252_polarization()
    @test polarization_at(polarization, 126) ≈ 0.41847
    @test width_at(polarization, 126) ≈ 0.50029

    # A mass outside the tabulated range falls back to the conventional means, which is a
    # documented behaviour rather than a silent zero.
    @test polarization_at(polarization, 300) == -0.5
    @test width_at(polarization, 300) == 0.6

    mean = mean_charge_polarization()
    @test polarization_at(mean, 126) == -0.5
    @test width_at(mean, 126) == 0.6
end

@testset "most probable charge" begin
    # Z_UCD(A) = A Z₀/A₀, exact at the symmetric split.
    @test unchanged_charge_distribution(CF252, 126) ≈ 49.0
    @test most_probable_charge(CF252, 126, 0.0) ≈ 49.0
    @test most_probable_charge(CF252, 126, 0.41847) ≈ 49.41847
end

@testset "isobaric charge distribution" begin
    # A normalized Gaussian: the peak sits at 1/(√(2π) rms) and it is symmetric about Zₚ.
    @test charge_probability(50.0, 50.0, 0.6) ≈ 1 / (sqrt(2π) * 0.6)
    @test charge_probability(49.0, 50.0, 0.6) ≈ charge_probability(51.0, 50.0, 0.6)

    # Summed over a wide enough integer range it approaches unity, but not exactly: sampling a
    # Gaussian of width 0.6 on an integer lattice overshoots by 2exp(−2π²rms²) ≈ 1.6e-3 by
    # Poisson summation. The excess is a property of the discretization, not of the expression.
    total = sum(charge_probability(Z, 50.0, 0.6) for Z in 40:60)
    @test total ≈ 1.0 atol = 3e-3
    @test total > 1.0

    @test_throws ArgumentError charge_probability(50.0, 50.0, 0.0)
    @test_throws ArgumentError charge_probability(50.0, 50.0, -0.6)
end

@testset "charge numbers per mass" begin
    @test charge_numbers(49.41847, 5) == 47:51
    @test charge_numbers(49.6, 1) == 50:50
    @test length(charge_numbers(49.4, 3)) == 3

    # An even count cannot be centred on the rounded most probable charge without biasing
    # the domain to one side.
    @test_throws ArgumentError charge_numbers(49.4, 4)
    @test_throws ArgumentError charge_numbers(49.4, 0)
end
