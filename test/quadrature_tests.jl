using AppleAccelerate: integrate
using AppleAccelerate.LibAccelerate: QUADRATURE_SUCCESS

@testset "Quadrature" begin
    @testset "polynomial / smooth (finite interval)" begin
        r = integrate(x -> x^2, 0, 1)
        @test r.value ≈ 1/3 atol=1e-8
        @test r.status == QUADRATURE_SUCCESS
        @test r.abserr < 1e-7

        @test integrate(sin, 0, π).value ≈ 2.0 atol=1e-8
        @test integrate(x -> 3x^2 + 2x + 1, -1, 2).value ≈ 15.0 atol=1e-8
    end

    @testset "endpoint singularity (QAGS)" begin
        # ∫₀¹ 1/√x dx = 2 — integrable singularity at the left endpoint
        @test integrate(x -> 1/sqrt(x), 0, 1).value ≈ 2.0 atol=1e-6
    end

    @testset "infinite bounds (QAGS)" begin
        @test integrate(x -> exp(-x^2), -Inf, Inf).value ≈ sqrt(π) atol=1e-7
        @test integrate(x -> exp(-x), 0, Inf).value ≈ 1.0 atol=1e-7
    end

    @testset "reversed bounds" begin
        # ∫₁⁰ x² dx = -∫₀¹ x² dx
        @test integrate(x -> x^2, 1, 0).value ≈ -1/3 atol=1e-8
    end

    @testset "integrator selection" begin
        for intg in (:qng, :qag, :qags)
            @test integrate(cos, 0, π/2; integrator = intg).value ≈ 1.0 atol=1e-7
        end
        @test_throws ArgumentError integrate(sin, 0, 1; integrator = :nope)
    end
end
