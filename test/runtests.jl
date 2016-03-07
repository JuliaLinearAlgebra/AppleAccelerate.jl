using AppleAccelerate

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

srand(7)
N = 1_000


for T in (Float32, Float64)
    @testset "Rounding::$T" begin
        X::Array{T} = 100*randn(N)
        @testset "Testing $f:$T" for f in [:floor,:ceil,:trunc,:round]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end
    end
end


for T in (Float32, Float64)
    @testset "Logarithmic::$T" begin
        X::Array{T} = exp(10*randn(N))
        @testset "Testing $f::$T" for f in [:log,:log2,:log10, :log1p]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end
    end
end


for T in (Float32, Float64)
    @testset "Exponential::$T" begin
        @testset "Testing $f::$T" for f in [:exp,:exp2,:expm1]
            X::Array{T} = 10*randn(N)
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end
    end
end


for T in (Float32, Float64)
    X::Array{T} = 10*randn(N)
    @testset "Trigonometric::$T" begin
        @testset "Testing $f::$T" for f in [:sin,:sinpi,:cos,:cospi,:tan,:atan] # tanpi not defined in Base
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end

        Y::Array{T} = 10*randn(N)
        @testset "Testing $f::$T" for f in [:atan2]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X,Y) ≈ fb(X,Y)
        end

        Z::Array{T} = 2*rand(N)-1
        @testset "Testing $f::$T" for f in [:asin,:acos]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Z) ≈ fb(Z)
        end
    end
end


for T in (Float32, Float64)
    @testset "Hyperbolic::$T" begin
        X = 10*randn(N)
        @testset "Testing $f::$T" for f in [:sinh,:cosh,:tanh,:asinh]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end

        Y = exp(10*randn(N))+1
        @testset "Testing $f::$T" for f in [:acosh]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Y) ≈ fb(Y)
        end

        Z = 2*rand(N)-1
        @testset "Testing $f::$T" for f in [:atanh]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Z) ≈ fb(Z)
        end
    end
end


@testset "DCT" begin
    r=rand(Float32,2^16)
    d1=dct(r)
    d2=AppleAccelerate.dct(r)
    @test norm(d1[2]/d2[2]*d2[2:end]-d1[2:end])≤1000eps(Float32)
end


for T in (Float32, Float64)
    @testset "Convolution & Correlation::$T" begin
        X::Array{T} = randn(N)
        Y::Array{T} = randn(N)
        @testset "Testing $f::$T" for f in [:conv, :xcorr]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fb(X, Y) ≈ fa(X, Y)
        end

        @testset "Testing $f::$T" for f in [:xcorr]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fb(X, copy(X)) ≈ fa(X)
        end
    end
end


for T in (Float32, Float64)
    @testset "Misc::$T" begin
        X::Array{T} = exp(10*randn(N))
        @testset "Testing $f::$T" for f in [:sqrt]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end

        Y::Array{T} = 10*randn(N)
        @testset "Testing $f::$T" for f in [:exponent, :abs]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Y) ≈ fb(Y)
        end

        Z::Array{T} = 10*randn(N)
        @testset "Testing $f::$T" for f in [:copysign]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X,Y) ≈ fb(X,Y)
        end

        @test AppleAccelerate.rem(X,Y) == [rem(X[i], Y[i]) for i=1:length(X)]

    end
end


for T in (Float32, Float64)
    @testset "Extra::$T" begin
        X::Array{T} = randn(N)
        Y::Array{T} = abs(randn(N))

        @test AppleAccelerate.rec(X) ≈ 1./X
        @test AppleAccelerate.rsqrt(Y) ≈ 1./sqrt(Y)
        @test AppleAccelerate.pow(Y,X) ≈ Y.^X
        @test AppleAccelerate.fdiv(X,Y) ≈ X./Y

        @test AppleAccelerate.sincos(X)[1] ≈ sin(X)
        @test AppleAccelerate.sincos(X)[2] ≈ cos(X)
        @test AppleAccelerate.cis(X) ≈ cis(X)

    end
end


for T in (Float32, Float64)
    @testset "Replace Base::$T" begin
        X::Array{T} = randn(N)
        Y::Array{T} = abs(randn(N))

        AppleAccelerate.@replaceBase(sin, atan2, ./)
        @test Base.sin(X) == AppleAccelerate.sin(X)
        @test Base.atan2(X, Y) == AppleAccelerate.atan2(X, Y)
        @test X ./ Y  == AppleAccelerate.fdiv(X, Y)

    end
end
