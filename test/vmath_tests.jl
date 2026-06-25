@testset "vmath / LibAccelerate migration" begin

    # The vmath subframework coverage notes are documented on this constant.
    @test isdefined(AppleAccelerate, :VMATH_COVERAGE)
    @test !isempty(string(@doc AppleAccelerate.VMATH_COVERAGE))

    # These tests double as a regression check for the internals migration of the
    # vForce (vv*) and vDSP wrappers in `array.jl` onto the generated
    # `LibAccelerate` raw layer: the idiomatic API must still match Base.
    for T in (Float32, Float64)
        @testset "vForce idiomatic API (migrated path)::$T" begin
            X = abs.(randn(T, 64)) .+ T(0.1)
            Y = abs.(randn(T, 64)) .+ T(0.1)

            @test AppleAccelerate.exp(X)   ≈ exp.(X)
            @test AppleAccelerate.log(X)   ≈ log.(X)
            @test AppleAccelerate.sqrt(X)  ≈ sqrt.(X)
            @test AppleAccelerate.sin(X)   ≈ sin.(X)
            @test AppleAccelerate.cos(X)   ≈ cos.(X)
            @test AppleAccelerate.cbrt(X)  ≈ cbrt.(X)
            @test AppleAccelerate.pow(X, Y) ≈ X .^ Y
            @test AppleAccelerate.pows(X, T(2)) ≈ X .^ T(2)
            @test AppleAccelerate.atan(X, Y) ≈ atan.(X, Y)

            s, c = AppleAccelerate.sincos(X)
            @test s ≈ sin.(X)
            @test c ≈ cos.(X)

            @test AppleAccelerate.cis(X) ≈ cis.(X)

            # mutating variant writes into the provided output buffer
            out = similar(X)
            @test AppleAccelerate.exp!(out, X) === out
            @test out ≈ exp.(X)
        end

        @testset "vDSP idiomatic API (migrated path)::$T" begin
            X = randn(T, 128)
            Y = randn(T, 128)
            @test AppleAccelerate.vadd(X, Y) ≈ X .+ Y
            @test AppleAccelerate.vsub(X, Y) ≈ X .- Y
            @test AppleAccelerate.vmul(X, Y) ≈ X .* Y
            @test AppleAccelerate.sum(X)     ≈ sum(X)
            @test AppleAccelerate.maximum(X) ≈ maximum(X)
            @test AppleAccelerate.vsadd(X, T(3)) ≈ X .+ T(3)
        end
    end
end
