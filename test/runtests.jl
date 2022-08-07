using AppleAccelerate
using DSP, Test, Random, LinearAlgebra, Statistics

if Sys.isapple()

    Random.seed!(7)
    N = 1_000

    for T in (Float32, Float64)
        @testset "Element-wise Operators::$T" begin
            X::Vector{T} = randn(N)
            Y::Vector{T} = randn(N)
            @test (X .+ Y) ≈ AppleAccelerate.vadd(X, Y)
            @test (X .- Y) ≈ AppleAccelerate.vsub(X, Y)
            @test (X .* Y) ≈ AppleAccelerate.vmul(X, Y)
            @test (X ./ Y) ≈ AppleAccelerate.vdiv(X, Y)
        end
    end

    for T in (Float32, Float64)
        @testset "Rounding::$T" begin
            X::Array{T} = 100*randn(N)
            @testset "Testing $f::$T" for f in [:floor,:ceil,:trunc,:round]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb.(X)
                @test fa.(X) ≈ fb.(X)
            end
        end
    end


    for T in (Float32, Float64)
        @testset "Logarithmic::$T" begin
            X::Array{T} = exp.(10*randn(N))
            @testset "Testing $f::$T" for f in [:log,:log2,:log10, :log1p]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb.(X)
                @test fa.(X) ≈ fb.(X)
            end
        end
    end


    for T in (Float32, Float64)
        @testset "Exponential::$T" begin
            @testset "Testing $f::$T" for f in [:exp,:exp2,:expm1]
                X::Array{T} = 10*randn(N)
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb.(X)
                @test fa.(X) ≈ fb.(X)
            end
        end
    end


    for T in (Float32, Float64)
        X::Array{T} = 10*randn(N)
        @testset "Trigonometric::$T" begin
            @testset "Testing $f::$T" for f in [:sin,:sinpi,:cos,:cospi,:tan,:atan] # tanpi not defined in Base
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb.(X)
                @test fa.(X) ≈ fb.(X)
            end

            Y::Array{T} = 10*randn(N)
            @testset "Testing $f::$T" for f in [:atan]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X,Y) ≈ fb.(X,Y)
                @test fa.(X,Y) ≈ fb.(X,Y)
            end

            Z::Array{T} = 2*rand(N) .- 1
            @testset "Testing $f::$T" for f in [:asin,:acos]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(Z) ≈ fb.(Z)
                @test fa.(Z) ≈ fb.(Z)
            end
        end
    end


    for T in (Float32, Float64)
        @testset "Hyperbolic::$T" begin
            X = 10*randn(N)
            @testset "Testing $f::$T" for f in [:sinh,:cosh,:tanh,:asinh]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb.(X)
                @test fa.(X) ≈ fb.(X)
            end

            Y = exp.(10*randn(N)) .+ 1
            @testset "Testing $f::$T" for f in [:acosh]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(Y) ≈ fb.(Y)
                @test fa.(Y) ≈ fb.(Y)
            end

            Z = 2*rand(N) .- 1
            @testset "Testing $f::$T" for f in [:atanh]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(Z) ≈ fb.(Z)
                @test fa.(Z) ≈ fb.(Z)
            end
        end
    end


    @testset "DCT::Float32" begin
        r=rand(Float32,2^16)
        d1=DSP.dct(r)
        plan_accel = AppleAccelerate.plan_dct(length(r), 2)
        d2=AppleAccelerate.dct(r, plan_accel)
        @test norm(d1[2]/d2[2]*d2[2:end]-d1[2:end])≤1000eps(Float32)
    end


    for T in (Float32,  Float64)
        @testset "Convolution & Correlation::$T" begin
            X::Vector{T} = randn(N)
            Y::Vector{T} = randn(N)
            @testset "Testing $f::$T" for f in [:conv, :xcorr]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X, Y) ≈ fb(X, Y)
            end

            @testset "Testing $f::$T" for f in [:xcorr]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb(X, copy(X))
            end
        end
    end

    for T in (Float64, )
        @testset "Biquadratic Flitering::$T" begin
            @testset "Single Section::$T" begin
                X::Vector{T} = randn(10)
                d::Vector{T} = zeros(4)
                c::Vector{T} = [x%0.5 for x in randn(5)]
                fdsp = DSP.Biquad(c[1], c[2], c[3], c[4], c[5])
                fa = AppleAccelerate.biquadcreate(c, 1)
                @test DSP.filt(fdsp, X) ≈ AppleAccelerate.biquad(X, d, length(X), fa)
            end
        end
    end


    for T in (Float32, Float64)
        @testset "Window Functions::$T" begin
            N = 64
            @testset "Testing blackman::$T" begin
                Wa = AppleAccelerate.blackman(N, T)
                Wb = 0.42 .- (0.5cos.(2pi.*(0:(N-1))./N)) .+ (0.08cos.(4pi.*(0:(N-1))./N))
                @test Wa ≈ Wb
            end

            @testset "Testing hamming::$T" begin
                Wa = AppleAccelerate.hamming(N, T)
                Wb = 0.54 .- 0.46cos.(2pi.*(0:(N-1))./N)
                @test Wa ≈ Wb
            end

            @testset "Testing hanning::$T" begin
                Wa = AppleAccelerate.hanning(N, T)
                Wb = 0.5(1 .- cos.(2pi.*(0:(N-1))./N))
                @test Wa ≈ Wb
            end
        end
    end

    for T in (Float32, Float64)
        @testset "Array Properties::$T" begin
            X::Array{T} = randn(N)
            @testset "Testing $f::$T" for f in [:maximum, :minimum, :mean, :sum]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb(X)
            end

            @testset "Testing $f::$T" for f in [:findmax, :findmin]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X)[1] ≈ fb(X)[1]
                @test fa(X)[2] ≈ fb(X)[2]
            end

            @testset "Testing meansqr::$T" begin
                @test AppleAccelerate.meansqr(X) ≈ mean(X .*X)
            end

            @testset "Testing meanmag::$T" begin
                @test AppleAccelerate.meanmag(X) ≈ mean(abs.(X))
            end

        end
    end

    for T in (Float32, Float64)
        @testset "Misc::$T" begin
            X::Array{T} = exp.(10*randn(N))
            @testset "Testing $f::$T" for f in [:sqrt]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X) ≈ fb.(X)
                @test fa.(X) ≈ fb.(X)
            end

            Y::Array{T} = 10*randn(N)
            @testset "Testing $f::$T" for f in [:exponent, :abs]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(Y) ≈ fb.(Y)
                @test fa.(Y) ≈ fb.(Y)
            end

            Z::Array{T} = 10*randn(N)
            @testset "Testing $f::$T" for f in [:copysign]
                @eval fb = $f
                @eval fa = AppleAccelerate.$f
                @test fa(X,Y) ≈ fb.(X,Y)
                @test fa.(X,Y) ≈ fb.(X,Y)
            end

            @test AppleAccelerate.rem(X,Y) == rem.(X, Y)

        end
    end


    for T in (Float32, Float64)
        @testset "Extra::$T" begin
            X::Array{T} = randn(N)
            Y::Array{T} = abs.(randn(N))

            @test AppleAccelerate.rec(X) ≈ 1 ./ X
            @test AppleAccelerate.rsqrt(Y) ≈ 1 ./ sqrt.(Y)
            @test AppleAccelerate.pow(Y,X) ≈ Y.^X
            @test AppleAccelerate.div_float(X,Y) ≈ X./Y

            @test AppleAccelerate.sincos(X)[1] ≈ sin.(X)
            @test AppleAccelerate.sincos(X)[2] ≈ cos.(X)
            @test AppleAccelerate.cis(X) ≈ cis.(X)

        end
    end

    #=
    AppleAccelerate.@replaceBase(sin, atan, /)

    @testset "Replace Base::$T" for T in (Float32, Float64)
    X::Array{T} = randn(N)
    Y::Array{T} = abs.(randn(N))

    @test Base.sin.(X) == AppleAccelerate.sin(X)
    @test Base.atan.(X, Y) == AppleAccelerate.atan(X, Y)
    @test X ./ Y  == AppleAccelerate.div_float(X, Y)
    end
    =#
    
end # Sys.isapple()

