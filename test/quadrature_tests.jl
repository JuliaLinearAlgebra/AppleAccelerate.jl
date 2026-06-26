using AppleAccelerate: integrate
using AppleAccelerate.LibAccelerate: QUADRATURE_SUCCESS

# Custom exception type for the throwing-integrand test (structs can't be
# defined in local/testset scope).
struct _QuadTestError <: Exception end

@testset "Quadrature" begin
    @testset "polynomial / smooth (finite interval)" begin
        # ∫₀¹ x² dx = 1/3
        r = integrate(x -> x^2, 0, 1)
        @test r.value ≈ 1/3 atol=1e-8
        @test r.status == QUADRATURE_SUCCESS
        @test r.abserr < 1e-7
        @test r.abserr ≥ 0

        # ∫₀^π sin x dx = 2
        @test integrate(sin, 0, π).value ≈ 2.0 atol=1e-8
        # ∫_{-1}^{2} (3x²+2x+1) dx = [x³+x²+x] = (8+4+2)-(-1+1-1) = 14-(-1) = 15
        @test integrate(x -> 3x^2 + 2x + 1, -1, 2).value ≈ 15.0 atol=1e-8
        # ∫₀¹ eˣ dx = e − 1
        @test integrate(exp, 0, 1).value ≈ (exp(1) - 1) atol=1e-8
        # ∫₀¹ 1/(1+x²) dx = π/4
        @test integrate(x -> 1/(1 + x^2), 0, 1).value ≈ π/4 atol=1e-8
    end

    @testset "endpoint singularity (QAGS)" begin
        # ∫₀¹ 1/√x dx = 2 — integrable singularity at the left endpoint
        @test integrate(x -> 1/sqrt(x), 0, 1).value ≈ 2.0 atol=1e-6
    end

    @testset "infinite and semi-infinite bounds (QAGS)" begin
        # ∫_{-∞}^{∞} e^{−x²} dx = √π
        @test integrate(x -> exp(-x^2), -Inf, Inf).value ≈ sqrt(π) atol=1e-7
        # ∫₀^∞ e^{−x} dx = 1   (semi-infinite, upper bound +Inf)
        rsemi = integrate(x -> exp(-x), 0, Inf)
        @test rsemi.value ≈ 1.0 atol=1e-7
        @test rsemi.status == QUADRATURE_SUCCESS
        @test rsemi.abserr < 1e-6
        # ∫_{-∞}^{0} e^{x} dx = 1   (semi-infinite, lower bound −Inf)
        @test integrate(x -> exp(x), -Inf, 0).value ≈ 1.0 atol=1e-7
        # ∫_{-∞}^{∞} 1/(1+x²) dx = π  (Cauchy / arctan)
        @test integrate(x -> 1/(1 + x^2), -Inf, Inf).value ≈ π atol=1e-7
    end

    @testset "reversed bounds" begin
        # ∫₁⁰ x² dx = −∫₀¹ x² dx = −1/3
        @test integrate(x -> x^2, 1, 0).value ≈ -1/3 atol=1e-8
        # Reversed semi-infinite: ∫_{∞}^{0} e^{−x} dx = −1
        @test integrate(x -> exp(-x), Inf, 0).value ≈ -1.0 atol=1e-7
    end

    @testset "integrator selection" begin
        # ∫₀^{π/2} cos x dx = 1, across every exposed algorithm.
        for intg in (:qng, :qag, :qags)
            r = integrate(cos, 0, π/2; integrator = intg)
            @test r.value ≈ 1.0 atol=1e-7
            @test r.status == QUADRATURE_SUCCESS
        end
        @test_throws ArgumentError integrate(sin, 0, 1; integrator = :nope)
    end

    @testset "qag_points validation" begin
        # Every valid Gauss-Kronrod rule integrates ∫₀^{π/2} cos x dx = 1 correctly.
        for pts in (15, 21, 31, 41, 51, 61)
            r = integrate(cos, 0, π/2; integrator = :qag, qag_points = pts)
            @test r.value ≈ 1.0 atol=1e-7
            @test r.status == QUADRATURE_SUCCESS
        end
        # qag_points = 0 (library default) is also valid.
        @test integrate(cos, 0, π/2; integrator = :qag, qag_points = 0).value ≈ 1.0 atol=1e-7
        # An invalid rule must throw (the library would otherwise return 0.0).
        @test_throws ArgumentError integrate(cos, 0, π/2; integrator = :qag, qag_points = 7)
        @test_throws ArgumentError integrate(cos, 0, π/2; integrator = :qag, qag_points = 100)
        @test_throws ArgumentError integrate(cos, 0, π/2; integrator = :qag, qag_points = 16)
    end

    @testset "custom tolerances and max_intervals" begin
        # Exercise the non-default keyword paths; ∫₀¹ x² dx = 1/3.
        r = integrate(x -> x^2, 0, 1; abstol = 1e-10, reltol = 1e-10, max_intervals = 50)
        @test r.value ≈ 1/3 atol=1e-9
        @test r.status == QUADRATURE_SUCCESS
    end

    @testset "type-unstable / non-Float64 integrand returns" begin
        # Integrand returning Int (∫₀¹ 1 dx = 1) — exercises the trampoline `convert`.
        @test integrate(x -> 1, 0, 1).value ≈ 1.0 atol=1e-8
        # Integrand returning a Float32 (∫₀¹ 2x dx = 1).
        @test integrate(x -> Float32(2x), 0, 1).value ≈ 1.0 atol=1e-6
        # Integrand returning a Rational (∫₀¹ x² dx = 1/3) via a type-unstable branch.
        @test integrate(x -> x < 0 ? 0 : x^2, 0, 1).value ≈ 1/3 atol=1e-8
        # Genuinely type-unstable: returns Int or Float64 depending on x.
        @test integrate(x -> x < 0.5 ? 1 : 1.0, 0, 1).value ≈ 1.0 atol=1e-8
    end

    @testset "throwing integrand surfaces a Julia exception" begin
        # The integrand error must come back as a Julia exception (not a crash
        # from unwinding through the C QUADPACK frames).
        @test_throws DomainError integrate(x -> x < 0.5 ? x : throw(DomainError(x)), 0, 1)
        # Custom exception type is preserved.
        @test_throws _QuadTestError integrate(x -> throw(_QuadTestError()), 0, 1)
        # Subsequent integrate calls still work (no corrupted state).
        @test integrate(x -> x^2, 0, 1).value ≈ 1/3 atol=1e-8
    end

    @testset "reentrancy (nested integrate)" begin
        # Outer integrand itself runs an inner integration:
        #   g(t) = ∫₀ᵗ 2x dx = t²   ⇒   ∫₀¹ g(t) dt = ∫₀¹ t² dt = 1/3
        g(t) = integrate(x -> 2x, 0, t).value
        @test integrate(g, 0, 1).value ≈ 1/3 atol=1e-7
        # Doubly nested: h(s) = ∫₀ˢ g(t) dt = ∫₀ˢ t² dt = s³/3
        #   ⇒ ∫₀¹ h(s) ds = ∫₀¹ s³/3 ds = 1/12
        h(s) = integrate(g, 0, s).value
        @test integrate(h, 0, 1).value ≈ 1/12 atol=1e-6
    end
end
