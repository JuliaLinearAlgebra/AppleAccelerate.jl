using AppleAccelerate
using Test

# Edge-case tests: empty inputs and NaN/Inf propagation.

@testset "Empty arrays" begin
    @testset "Real ops ($T)" for T in (Float32, Float64)
        e = T[]

        # Element-wise ops return empty results
        @test AppleAccelerate.vadd(e, e) == T[]
        @test AppleAccelerate.vsub(e, e) == T[]
        @test AppleAccelerate.vmul(e, e) == T[]
        @test AppleAccelerate.vneg(e) == T[]
        @test AppleAccelerate.vabs(e) == T[]
        @test AppleAccelerate.exp(e) == T[]

        # Sum-type reductions match Base on empty input
        @test AppleAccelerate.sum(e) === zero(T)
        @test AppleAccelerate.dot(e, e) === zero(T)

        # vDSP would silently return ±Inf (and index 0) here; we guard and
        # throw ArgumentError to match Base behavior.
        @test_throws ArgumentError AppleAccelerate.maximum(e)
        @test_throws ArgumentError AppleAccelerate.minimum(e)
        @test_throws ArgumentError AppleAccelerate.findmax(e)
        @test_throws ArgumentError AppleAccelerate.findmin(e)
    end

    @testset "Complex ops ($T)" for T in (Float32, Float64)
        ec = Complex{T}[]
        e = T[]

        @test AppleAccelerate.vneg(ec) == Complex{T}[]
        @test AppleAccelerate.vconj(ec) == Complex{T}[]
        @test AppleAccelerate.vmul(ec, ec) == Complex{T}[]
        @test AppleAccelerate.vabs(ec) == T[]
        @test AppleAccelerate.zvadd(ec, ec) == Complex{T}[]
        @test AppleAccelerate.zvsub(ec, ec) == Complex{T}[]
        @test AppleAccelerate.dot(ec, ec) === zero(Complex{T})
        @test AppleAccelerate.dotu(ec, ec) === zero(Complex{T})
        @test AppleAccelerate.polar(ec) == (T[], T[])
        @test AppleAccelerate.ctoz(ec) == (T[], T[])
        @test AppleAccelerate.ztoc(e, e) == Complex{T}[]
    end
end

@testset "NaN/Inf propagation" begin
    @testset "$T" for T in (Float32, Float64)
        # exp: NaN → NaN, Inf → Inf, -Inf → 0 (matches Base)
        X = T[1.0, NaN, Inf, -Inf, 0.0]
        @test isequal(AppleAccelerate.exp(X), Base.exp.(X))

        # sqrt: NaN → NaN, Inf → Inf (matches Base on its valid domain)
        X = T[NaN, Inf, 0.0, 4.0]
        @test isequal(AppleAccelerate.sqrt(X), Base.sqrt.(X))

        # log: NaN → NaN, Inf → Inf (matches Base on its valid domain)
        X = T[NaN, Inf, 1.0]
        @test isequal(AppleAccelerate.log(X), Base.log.(X))

        # sin: NaN → NaN; finite values match Base
        X = T[NaN, 0.5, -1.25]
        @test isnan(AppleAccelerate.sin(X)[1])
        @test AppleAccelerate.sin(X)[2:3] ≈ Base.sin.(X[2:3])
        # Base.sin(Inf) throws DomainError; Accelerate returns NaN instead
        @test isnan(AppleAccelerate.sin(T[Inf])[1])
    end
end
