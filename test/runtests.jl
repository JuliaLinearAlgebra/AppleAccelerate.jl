using AppleAccelerate

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

srand(7)
N = 1_000

@testset "Rounding" begin
    X = 100*randn(N)
    @testset "Testing $f" for f in [:floor,:ceil,:trunc,:round]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X) ≈ fb(X)
    end
end


@testset "Logarithmic" begin
    X = exp(10*randn(N))
    @testset "Testing $f" for f in [:log,:log2,:log10]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X) ≈ fb(X)
    end

    Y = expm1(10*randn(N))
    @testset "Testing $f" for f in [:log1p]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(Y) ≈ fb(Y)
    end
end


@testset "Exponential" begin
    X = 100*randn(N)
    @testset "Testing $f" for f in [:exp,:exp2,:expm1]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X) ≈ fb(X)
    end
end


@testset "Trigonometric" begin
    X = 10*randn(N)
    @testset "Testing $f" for f in [:sin,:sinpi,:cos,:cospi,:tan,:atan] # tanpi not defined in Base
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X) ≈ fb(X)
    end

    Y = 10*randn(N)
    @testset "Testing $f" for f in [:atan2]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X,Y) ≈ fb(X,Y)
    end

    Z = 2*rand(N)-1
    @testset "Testing $f" for f in [:asin,:acos]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(Z) ≈ fb(Z)
    end
end


@testset "Hyperbolic" begin
    X = 10*randn(N)
    @testset "Testing $f" for f in [:sinh,:cosh,:tanh,:asinh]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X) ≈ fb(X)
    end

    Y = exp(10*randn(N))+1
    @testset "Testing $f" for f in [:acosh]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(Y) ≈ fb(Y)
    end

    Z = 2*rand(N)-1
    @testset "Testing $f" for f in [:atanh]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(Z) ≈ fb(Z)
    end
end


@testset "DCT" begin
    r=rand(Float32,2^16)
    d1=dct(r)
    d2=AppleAccelerate.dct(r)
    @test norm(d1[2]/d2[2]*d2[2:end]-d1[2:end])≤1000eps(Float32)
end


@testset "Misc" begin
    X = exp(10*randn(N))
    @testset "Testing $f" for f in [:sqrt]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X) ≈ fb(X)
    end

    Y = 10*randn(N)
    @testset "Testing $f" for f in [:exponent, :abs]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(Y) ≈ fb(Y)
    end

    Z = 10*randn(N)
    @testset "Testing $f" for f in [:copysign]
        @eval fb = $f
        @eval fa = AppleAccelerate.$f
        @test fa(X,Y) ≈ fb(X,Y)
    end

    # @fact AppleAccelerate.rem(X,Y) => [rem(X[i], Y[i]) for i=1:length(X)] # no vectorized rem

end


@testset "Extra" begin
    X = randn(N)
    Y = abs(randn(N))

    @test AppleAccelerate.rec(X) ≈ 1./X
    @test AppleAccelerate.rsqrt(Y) ≈ 1./sqrt(Y)
    @test AppleAccelerate.pow(Y,X) ≈ Y.^X
    @test AppleAccelerate.fdiv(X,Y) ≈ X./Y

#    @test [AppleAccelerate.sincos(X)...] ≈ [sin(X);cos(X)]
    @test AppleAccelerate.cis(X) ≈ cis(X)

end


@testset "Replace Base" begin
    X = randn(N)
    Y = abs(randn(N))

    AppleAccelerate.@replaceBase(sin,atan2,./,.^)
    # @test sin(X) => AppleAccelerate.sin(X)
    # @test atan2(X,Y) => AppleAccelerate.atan2(X,Y)
    # @test X ./ Y => AppleAccelerate.fdiv(X,Y)
    # @test Y .^ X => AppleAccelerate.pow(Y,X)
end
