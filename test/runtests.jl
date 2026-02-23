using LinearAlgebra
using AppleAccelerate
using AbstractFFTs
using DSP, FFTW, Test, Random, Statistics

if !Sys.isapple()
    @info("AppleAccelerate.jl will be tested only on macOS. Exiting.")
    exit(0)
end

include("libSparseTests.jl")

Random.seed!(7)
N = 1_000

@testset "AppleAccelerate.jl" begin
for T in (Float32, Float64)
    @testset "Element-wise Operators::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)
        # Vector-vector
        @test (X .+ Y) ≈ AppleAccelerate.vadd(X, Y)
        @test (X .- Y) ≈ AppleAccelerate.vsub(X, Y)
        @test (X .* Y) ≈ AppleAccelerate.vmul(X, Y)
        @test (X ./ Y) ≈ AppleAccelerate.vdiv(X, Y)

        # Vector-vector non-allocating
        AppleAccelerate.vadd!(Z, X, Y)
        @test (X .+ Y) ≈ Z
        AppleAccelerate.vsub!(Z, X, Y)
        @test (X .- Y) ≈ Z
        AppleAccelerate.vmul!(Z, X, Y)
        @test (X .* Y) ≈ Z
        AppleAccelerate.vdiv!(Z, X, Y)
        @test (X ./ Y) ≈ Z

        # Vector-vector broadcasting
        @test (X .+ Y) ≈ AppleAccelerate.vadd.(X, Y)
        @test (X .- Y) ≈ AppleAccelerate.vsub.(X, Y)
        @test (X .* Y) ≈ AppleAccelerate.vmul.(X, Y)
        @test (X ./ Y) ≈ AppleAccelerate.vdiv.(X, Y)

        #Vector-scalar
        c::T         = randn()
        @test (X .+ c) ≈ AppleAccelerate.vsadd.(X, c)
        @test (X .- c) ≈ AppleAccelerate.vssub.(X, c)
        @test (c .- X) ≈ AppleAccelerate.svsub.(X, c)
        @test (X .* c) ≈ AppleAccelerate.vsmul.(X, c)
        @test (X ./ c) ≈ AppleAccelerate.vsdiv.(X, c)

        #Vector-scalar non-allocating
        AppleAccelerate.vsadd!(Y, X, c)
        @test (X .+ c) ≈ Y
        AppleAccelerate.vssub!(Y, X, c)
        @test (X .- c) ≈ Y
        AppleAccelerate.svsub!(Y, X, c)
        @test (c .- X) ≈ Y
        AppleAccelerate.vsmul!(Y, X, c)
        @test (X .* c) ≈ Y
        AppleAccelerate.vsdiv!(Y, X, c)
        @test (X ./ c) ≈ Y

        #Vector-scalar broadcasting
        @test (X .+ c) ≈ AppleAccelerate.vsadd.(X, c)
        @test (X .- c) ≈ AppleAccelerate.vssub.(X, c)
        @test (c .- X) ≈ AppleAccelerate.svsub.(X, c)
        @test (X .* c) ≈ AppleAccelerate.vsmul.(X, c)
        @test (X ./ c) ≈ AppleAccelerate.vsdiv.(X, c)

        @test (X .+ Y .+ Y) ≈ AppleAccelerate.vadd.(X, Y .+ Y)
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

@testset "fft::Float64" begin
    for n in (16, 256, 1024)
        r = randn(ComplexF64, n)
        @test AppleAccelerate.fft(r) ≈ FFTW.fft(r)
        # With explicit plan
        setup = AppleAccelerate.plan_fft(r)
        @test AppleAccelerate.fft(r, setup) ≈ FFTW.fft(r)
    end
end

@testset "fft::Float32" begin
    for n in (16, 256, 1024)
        r = randn(ComplexF32, n)
        @test AppleAccelerate.fft(r) ≈ FFTW.fft(r) rtol=sqrt(eps(Float32))
    end
end

@testset "fft known values" begin
    # DFT of unit impulse is all ones
    @test AppleAccelerate.fft(ComplexF64[1, 0, 0, 0]) ≈ ComplexF64[1, 1, 1, 1]
    # DFT of constant is scaled impulse
    @test AppleAccelerate.fft(ComplexF64[1, 1, 1, 1]) ≈ ComplexF64[4, 0, 0, 0]
    # DFT of single frequency: e^{2πi k/N} for k=0..N-1 has DFT = N δ_{1}
    r = [exp(2π * im * k / 4) for k in 0:3]
    result = AppleAccelerate.fft(ComplexF64.(r))
    @test abs(result[2]) ≈ 4.0 atol=1e-12
    @test abs(result[1]) < 1e-12
    @test abs(result[3]) < 1e-12
    @test abs(result[4]) < 1e-12

    # Parseval's theorem: sum(|x|^2) == sum(|X|^2) / N
    x = randn(ComplexF64, 256)
    X = AppleAccelerate.fft(x)
    @test sum(abs2, x) ≈ sum(abs2, X) / 256

    # Linearity: FFT(a*x + b*y) == a*FFT(x) + b*FFT(y)
    x1 = randn(ComplexF64, 128)
    x2 = randn(ComplexF64, 128)
    a, b = 3.0 + 1.0im, -2.0 + 0.5im
    @test AppleAccelerate.fft(ComplexF64.(a .* x1 .+ b .* x2)) ≈
          a .* AppleAccelerate.fft(x1) .+ b .* AppleAccelerate.fft(x2)

    # Circular shift property: left shift by m multiplies FFT by e^{+2πi m k/N}
    x = randn(ComplexF64, 64)
    m = 7
    x_shifted = circshift(x, -m)
    X = AppleAccelerate.fft(x)
    X_shifted = AppleAccelerate.fft(x_shifted)
    phase = [exp(2π * im * m * k / 64) for k in 0:63]
    @test X_shifted ≈ X .* phase
end

@testset "bfft and ifft roundtrip" begin
    for T in (ComplexF64, ComplexF32)
        F = real(T)
        @testset "$T" begin
            for n in (16, 256, 1024)
                r = randn(T, n)
                fwd = AppleAccelerate.fft(r)
                # bfft is unnormalized inverse
                if F == Float64
                    @test AppleAccelerate.bfft(fwd) ./ n ≈ r
                else
                    @test AppleAccelerate.bfft(fwd) ./ n ≈ r rtol=sqrt(eps(Float32))
                end
                # ifft is normalized inverse
                if F == Float64
                    @test AppleAccelerate.ifft(fwd) ≈ r
                else
                    @test AppleAccelerate.ifft(fwd) ≈ r rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "fft 2D::Float64" begin
    for (nr, nc) in ((4, 4), (8, 16), (16, 8), (32, 32))
        r = randn(ComplexF64, nr, nc)
        @test AppleAccelerate.fft(r) ≈ FFTW.fft(r)
        # With explicit plan
        setup = AppleAccelerate.plan_fft(r)
        @test AppleAccelerate.fft(r, setup) ≈ FFTW.fft(r)
    end
end

@testset "fft 2D::Float32" begin
    for (nr, nc) in ((4, 4), (8, 16), (16, 8), (32, 32))
        r = randn(ComplexF32, nr, nc)
        @test AppleAccelerate.fft(r) ≈ FFTW.fft(r) rtol=sqrt(eps(Float32))
    end
end

@testset "bfft and ifft 2D roundtrip" begin
    for T in (ComplexF64, ComplexF32)
        F = real(T)
        @testset "$T" begin
            for (nr, nc) in ((4, 8), (16, 16))
                r = randn(T, nr, nc)
                fwd = AppleAccelerate.fft(r)
                ntotal = nr * nc
                if F == Float64
                    @test AppleAccelerate.bfft(fwd) ./ ntotal ≈ r
                    @test AppleAccelerate.ifft(fwd) ≈ r
                else
                    @test AppleAccelerate.bfft(fwd) ./ ntotal ≈ r rtol=sqrt(eps(Float32))
                    @test AppleAccelerate.ifft(fwd) ≈ r rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "fft! in-place" begin
    for T in (ComplexF64, ComplexF32)
        F = real(T)
        @testset "1D $T" begin
            for n in (16, 256, 1024)
                x = randn(T, n)
                expected = AppleAccelerate.fft(x)
                y = copy(x)
                AppleAccelerate.fft!(y)
                if F == Float64
                    @test y ≈ expected
                else
                    @test y ≈ expected rtol=sqrt(eps(Float32))
                end
            end
        end
        @testset "2D $T" begin
            for (nr, nc) in ((4, 4), (8, 16))
                x = randn(T, nr, nc)
                expected = AppleAccelerate.fft(x)
                y = copy(x)
                AppleAccelerate.fft!(y)
                if F == Float64
                    @test y ≈ expected
                else
                    @test y ≈ expected rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "bfft! and ifft! in-place" begin
    for T in (ComplexF64, ComplexF32)
        F = real(T)
        @testset "$T" begin
            for n in (16, 256)
                x = randn(T, n)
                fwd = AppleAccelerate.fft(x)
                y = copy(fwd)
                AppleAccelerate.ifft!(y)
                if F == Float64
                    @test y ≈ x
                else
                    @test y ≈ x rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "rfft" begin
    for T in (Float64, Float32)
        @testset "$T" begin
            for n in (16, 256, 1024)
                x = randn(T, n)
                result = AppleAccelerate.rfft(x)
                expected = FFTW.rfft(x)
                if T == Float64
                    @test result ≈ expected
                else
                    @test result ≈ expected rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "brfft and irfft roundtrip" begin
    for T in (Float64, Float32)
        @testset "$T" begin
            for n in (16, 256, 1024)
                x = randn(T, n)
                X = AppleAccelerate.rfft(x)
                if T == Float64
                    @test AppleAccelerate.brfft(X, n) ./ n ≈ x
                    @test AppleAccelerate.irfft(X, n) ≈ x
                else
                    @test AppleAccelerate.brfft(X, n) ./ n ≈ x rtol=sqrt(eps(Float32))
                    @test AppleAccelerate.irfft(X, n) ≈ x rtol=sqrt(eps(Float32))
                end
            end
        end
    end
end

@testset "AbstractFFTs extension" begin
    @testset "Forward FFT matches FFTW" begin
        for T in (ComplexF64, ComplexF32)
            F = real(T)
            @testset "$T" begin
                for n in (16, 64, 256, 1024)
                    x = randn(T, n)
                    result = AbstractFFTs.fft(x)
                    expected = FFTW.fft(x)
                    if F == Float64
                        @test result ≈ expected
                    else
                        @test result ≈ expected rtol=sqrt(eps(Float32))
                    end
                end
            end
        end
    end

    @testset "mul! variant" begin
        x = randn(ComplexF64, 256)
        p = plan_fft(x)
        y = similar(x)
        mul!(y, p, x)
        @test y ≈ FFTW.fft(x)
    end

    @testset "bfft roundtrip" begin
        for T in (ComplexF64, ComplexF32)
            F = real(T)
            @testset "$T" begin
                for n in (16, 256, 1024)
                    x = randn(T, n)
                    X = AbstractFFTs.fft(x)
                    y = AbstractFFTs.bfft(X)
                    if F == Float64
                        @test y ≈ n * x
                    else
                        @test y ≈ n * x rtol=sqrt(eps(Float32))
                    end
                end
            end
        end
    end

    @testset "ifft roundtrip" begin
        for T in (ComplexF64, ComplexF32)
            F = real(T)
            @testset "$T" begin
                for n in (16, 256, 1024)
                    x = randn(T, n)
                    @test AbstractFFTs.ifft(AbstractFFTs.fft(x)) ≈ x rtol=(F == Float64 ? 1e-10 : sqrt(eps(Float32)))
                end
            end
        end
    end

    @testset "inv(plan) and p \\ x" begin
        x = randn(ComplexF64, 256)
        p = plan_fft(x)
        X = p * x
        pinv = inv(p)
        @test pinv * X ≈ x
        @test inv(p) === pinv
        @test p \ X ≈ x
    end

    @testset "Known values" begin
        x_impulse = ComplexF64[1, 0, 0, 0]
        @test AbstractFFTs.fft(x_impulse) ≈ ComplexF64[1, 1, 1, 1]
        x_const = ComplexF64[1, 1, 1, 1]
        @test AbstractFFTs.fft(x_const) ≈ ComplexF64[4, 0, 0, 0]
    end

    @testset "Parseval's theorem" begin
        x = randn(ComplexF64, 256)
        X = AbstractFFTs.fft(x)
        @test sum(abs2, x) ≈ sum(abs2, X) / 256
    end

    @testset "2D Forward FFT matches FFTW" begin
        for T in (ComplexF64, ComplexF32)
            F = real(T)
            @testset "$T" begin
                for (nr, nc) in ((4, 4), (8, 16), (16, 8), (32, 32))
                    x = randn(T, nr, nc)
                    result = AbstractFFTs.fft(x)
                    expected = FFTW.fft(x)
                    if F == Float64
                        @test result ≈ expected
                    else
                        @test result ≈ expected rtol=sqrt(eps(Float32))
                    end
                end
            end
        end
    end

    @testset "2D mul! variant" begin
        x = randn(ComplexF64, 16, 32)
        p = plan_fft(x)
        y = similar(x)
        mul!(y, p, x)
        @test y ≈ FFTW.fft(x)
    end

    @testset "2D ifft roundtrip" begin
        for T in (ComplexF64, ComplexF32)
            F = real(T)
            @testset "$T" begin
                for (nr, nc) in ((4, 8), (16, 16))
                    x = randn(T, nr, nc)
                    @test AbstractFFTs.ifft(AbstractFFTs.fft(x)) ≈ x rtol=(F == Float64 ? 1e-10 : sqrt(eps(Float32)))
                end
            end
        end
    end

    @testset "2D inv(plan) and p \\ x" begin
        x = randn(ComplexF64, 16, 16)
        p = plan_fft(x)
        X = p * x
        pinv = inv(p)
        @test pinv * X ≈ x
        @test p \ X ≈ x
    end

    @testset "Error cases" begin
        # Non-power-of-2
        @test_throws ArgumentError plan_fft(randn(ComplexF64, 100))
        @test_throws ArgumentError fft(randn(ComplexF64, 100))
        # Empty
        @test_throws ArgumentError plan_fft(ComplexF64[])
        # Non-power-of-2 2D
        @test_throws ArgumentError plan_fft(randn(ComplexF64, 3, 4))
        @test_throws ArgumentError plan_fft(randn(ComplexF64, 4, 6))
        # 3D+ arrays
        @test_throws ArgumentError plan_fft(randn(ComplexF64, 4, 4, 4))
        @test_throws ArgumentError fft(randn(ComplexF64, 4, 4, 4))
        @test_throws ArgumentError plan_fft!(randn(ComplexF64, 4, 4, 4))
    end

    @testset "In-place fft!" begin
        x = randn(ComplexF64, 256)
        p = plan_fft!(x)
        y = copy(x)
        p * y
        @test y ≈ FFTW.fft(x)
    end

    @testset "In-place ifft! roundtrip" begin
        x = randn(ComplexF64, 256)
        fwd = copy(x)
        fft!(fwd)
        ifft!(fwd)
        @test fwd ≈ x
    end

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

        @testset "Testing meanmag::$T" begin
            @test AppleAccelerate.meanmag(X) ≈ mean(abs, X)
        end

        @testset "Testing meansqr::$T" begin
            @test AppleAccelerate.meansqr(X) ≈ mean(X .* X)
        end

        @testset "Testing meanssqr::$T" begin
            @test AppleAccelerate.meanssqr(X) ≈ mean(X .* abs.(X))
        end

        @testset "Testing summag::$T" begin
            @test AppleAccelerate.summag(X) ≈ sum(abs, X)
        end

        @testset "Testing sumsqr::$T" begin
            @test AppleAccelerate.sumsqr(X) ≈ sum(abs2, X)
        end

        @testset "Testing sumssqr::$T" begin
            @test AppleAccelerate.sumssqr(X) ≈ sum(X .* abs.(X))
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

for T in (Float32, Float64)
    @testset "vDSP Unary::$T" begin
        X::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)

        @test AppleAccelerate.vneg(X) ≈ -X
        AppleAccelerate.vneg!(Z, X)
        @test Z ≈ -X

        @test AppleAccelerate.vnabs(X) ≈ -abs.(X)
        AppleAccelerate.vnabs!(Z, X)
        @test Z ≈ -abs.(X)

        @test AppleAccelerate.vsq(X) ≈ X .^ 2
        AppleAccelerate.vsq!(Z, X)
        @test Z ≈ X .^ 2

        @test AppleAccelerate.vssq(X) ≈ X .* abs.(X)
        AppleAccelerate.vssq!(Z, X)
        @test Z ≈ X .* abs.(X)

        @test AppleAccelerate.vfrac(X) ≈ X .- trunc.(X)
        AppleAccelerate.vfrac!(Z, X)
        @test Z ≈ X .- trunc.(X)

        @test AppleAccelerate.vabs(X) ≈ abs.(X)
        AppleAccelerate.vabs!(Z, X)
        @test Z ≈ abs.(X)

        Y = copy(X)
        AppleAccelerate.vreverse!(Y)
        @test Y ≈ reverse(X)
        @test AppleAccelerate.vreverse(X) ≈ reverse(X)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Two-Vector::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)

        @test AppleAccelerate.vmax(X, Y) ≈ max.(X, Y)
        AppleAccelerate.vmax!(Z, X, Y)
        @test Z ≈ max.(X, Y)

        @test AppleAccelerate.vmin(X, Y) ≈ min.(X, Y)
        AppleAccelerate.vmin!(Z, X, Y)
        @test Z ≈ min.(X, Y)

        @test AppleAccelerate.vmaxmg(X, Y) ≈ max.(abs.(X), abs.(Y))
        AppleAccelerate.vmaxmg!(Z, X, Y)
        @test Z ≈ max.(abs.(X), abs.(Y))

        @test AppleAccelerate.vminmg(X, Y) ≈ min.(abs.(X), abs.(Y))
        AppleAccelerate.vminmg!(Z, X, Y)
        @test Z ≈ min.(abs.(X), abs.(Y))

        @test AppleAccelerate.vdist(X, Y) ≈ hypot.(X, Y)
        AppleAccelerate.vdist!(Z, X, Y)
        @test Z ≈ hypot.(X, Y)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Scalar-Vector Divide::$T" begin
        X::Vector{T} = randn(N) .+ T(2)  # avoid division by zero
        c::T = T(3.5)
        Z::Vector{T} = similar(X)

        @test AppleAccelerate.svdiv(X, c) ≈ c ./ X
        AppleAccelerate.svdiv!(Z, X, c)
        @test Z ≈ c ./ X
    end
end

for T in (Float32, Float64)
    @testset "vDSP Compound Arithmetic::$T" begin
        A::Vector{T} = randn(N)
        B::Vector{T} = randn(N)
        C::Vector{T} = randn(N)
        D::Vector{T} = randn(N)
        R::Vector{T} = similar(A)

        # vam: (A+B)*C
        @test AppleAccelerate.vam(A, B, C) ≈ (A .+ B) .* C
        AppleAccelerate.vam!(R, A, B, C)
        @test R ≈ (A .+ B) .* C

        # vsbm: (A-B)*C
        @test AppleAccelerate.vsbm(A, B, C) ≈ (A .- B) .* C
        AppleAccelerate.vsbm!(R, A, B, C)
        @test R ≈ (A .- B) .* C

        # vaam: (A+B)*(C+D)
        @test AppleAccelerate.vaam(A, B, C, D) ≈ (A .+ B) .* (C .+ D)
        AppleAccelerate.vaam!(R, A, B, C, D)
        @test R ≈ (A .+ B) .* (C .+ D)

        # vsbsbm: (A-B)*(C-D)
        @test AppleAccelerate.vsbsbm(A, B, C, D) ≈ (A .- B) .* (C .- D)
        AppleAccelerate.vsbsbm!(R, A, B, C, D)
        @test R ≈ (A .- B) .* (C .- D)

        # vasbm: (A+B)*(C-D)
        @test AppleAccelerate.vasbm(A, B, C, D) ≈ (A .+ B) .* (C .- D)
        AppleAccelerate.vasbm!(R, A, B, C, D)
        @test R ≈ (A .+ B) .* (C .- D)

        # vpythg: sqrt((A-C)^2 + (B-D)^2)
        @test AppleAccelerate.vpythg(A, B, C, D) ≈ sqrt.((A .- C) .^ 2 .+ (B .- D) .^ 2)

        # vasm: (A+B)*c
        c::T = T(2.5)
        @test AppleAccelerate.vasm(A, B, c) ≈ (A .+ B) .* c

        # vsbsm: (A-B)*c
        @test AppleAccelerate.vsbsm(A, B, c) ≈ (A .- B) .* c

        # vsma: A*b + C
        b::T = T(1.5)
        @test AppleAccelerate.vsma(A, b, C) ≈ A .* b .+ C

        # vsmsa: A*b + c
        @test AppleAccelerate.vsmsa(A, b, c) ≈ A .* b .+ c

        # vaddsub: returns (A+B, B-A)
        (add_r, sub_r) = AppleAccelerate.vaddsub(A, B)
        @test add_r ≈ A .+ B
        @test sub_r ≈ B .- A
    end
end

for T in (Float32, Float64)
    @testset "vDSP Reductions::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)

        @test AppleAccelerate.dot(X, Y) ≈ LinearAlgebra.dot(X, Y)
        @test AppleAccelerate.distancesq(X, Y) ≈ sum((X .- Y) .^ 2)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Clipping::$T" begin
        X::Vector{T} = randn(N)

        low::T = T(-0.5)
        high::T = T(0.5)
        @test AppleAccelerate.vclip(X, low, high) ≈ clamp.(X, low, high)
        Z = similar(X)
        AppleAccelerate.vclip!(Z, X, low, high)
        @test Z ≈ clamp.(X, low, high)

        # vthr
        thr::T = T(0.0)
        @test AppleAccelerate.vthr(X, thr) ≈ max.(X, thr)

        # vthres
        expected_vthres = [x >= thr ? x : T(0) for x in X]
        @test AppleAccelerate.vthres(X, thr) ≈ expected_vthres
    end
end

@testset "vDSP Type Conversion" begin
    X32 = randn(Float32, N)
    X64 = randn(Float64, N)

    @test AppleAccelerate.vdouble(X32) ≈ Float64.(X32)
    @test AppleAccelerate.vsingle(X64) ≈ Float32.(X64)
end

for T in (Float32, Float64)
    @testset "vDSP Ramp::$T" begin
        start::T = T(1.0)
        step::T = T(0.5)
        n = 100
        ramp = AppleAccelerate.vramp(start, step, n)
        expected = [start + i * step for i in 0:(n-1)]
        @test ramp ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP Integration::$T" begin
        n = 101
        X::Vector{T} = ones(T, n)
        step::T = T(0.1)

        # vtrapz: trapezoidal integration of constant 1 over [0, 10] should give cumulative
        trapz = AppleAccelerate.vtrapz(X, step)
        # Result has same length as X; trapz[1] = 0, trapz[end] ≈ (n-1)*step
        @test trapz[1] ≈ T(0) atol=eps(T)
        @test trapz[end] ≈ T((n - 1) * step) atol=T(0.01)

        # vswsum: sliding window sum
        X2::Vector{T} = ones(T, 20)
        sw = AppleAccelerate.vswsum(X2, 5)
        @test length(sw) == 16
        @test all(x -> isapprox(x, T(5), atol=eps(T)), sw)

        # vswmax: sliding window max
        X3::Vector{T} = T.(collect(1:20))
        swm = AppleAccelerate.vswmax(X3, 5)
        @test length(swm) == 16
        @test swm ≈ T.(collect(5:20))
    end
end

for T in (Float32, Float64)
    @testset "vDSP Interpolation::$T" begin
        A::Vector{T} = randn(N)
        B::Vector{T} = randn(N)
        t::T = T(0.3)

        # vintb: A + t*(B-A)
        @test AppleAccelerate.vintb(A, B, t) ≈ A .+ t .* (B .- A)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Polynomial::$T" begin
        # Evaluate p(x) = 2x^2 + 3x + 1
        # vDSP_vpoly uses coefficients in reverse order: [highest degree first]
        # coeffs = [a_P, a_{P-1}, ..., a_1, a_0] for P-degree poly
        coeffs::Vector{T} = T[2, 3, 1]
        X::Vector{T} = T.(collect(-5.0:0.5:5.0))
        result = AppleAccelerate.vpoly(coeffs, X)
        expected = 2 .* X .^ 2 .+ 3 .* X .+ 1
        @test result ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP Normalization::$T" begin
        X::Vector{T} = randn(N)
        (result, m, s) = AppleAccelerate.vnormalize(X)
        @test m ≈ Statistics.mean(X) atol=T(1e-4)
        if s > 0
            @test result ≈ (X .- m) ./ s atol=T(1e-3)
        end
    end
end

for T in (Float32, Float64)
    @testset "vDSP Decibel Conversion::$T" begin
        X::Vector{T} = abs.(randn(N)) .+ T(0.01)  # positive values
        ref::T = T(1.0)

        # Power: 10*log10(X/ref)
        result_pow = AppleAccelerate.vdbcon(X, ref, true)
        expected_pow = 10 .* log10.(X ./ ref)
        @test result_pow ≈ expected_pow atol=T(1e-3)

        # Amplitude: 20*log10(X/ref)
        result_amp = AppleAccelerate.vdbcon(X, ref, false)
        expected_amp = 20 .* log10.(X ./ ref)
        @test result_amp ≈ expected_amp atol=T(1e-3)
    end
end

end

if AppleAccelerate.get_macos_version() < v"13.4"
    @info("AppleAccelerate.jl needs macOS >= 13.4 for BLAS forwarding. Not testing forwarding capabilities.")
    exit(0)
end

# Set up a debugging fallback function that prints out a stacktrace if the LinearAlgebra
# tests end up calling a function that we don't have forwarded.
accelerate_header = """
import LinearAlgebra
function debug_missing_function()
    println("Missing BLAS/LAPACK function!")
    display(stacktrace())
end
LinearAlgebra.BLAS.lbt_set_default_func(@cfunction(debug_missing_function, Cvoid, ()))
import AppleAccelerate
"""
eval(Meta.parseall(accelerate_header))

config = BLAS.get_config()
@info("Running with $(length(config.loaded_libs)) libraries loaded:")
display(config.loaded_libs)

using Test
@testset "Accelerate Forwarding Sanity Tests" begin
    ver = AppleAccelerate.get_macos_version()
    if ver >= v"13.4"
        @test LinearAlgebra.peakflops() > 0
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :lp64).libname, "Accelerate")
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :ilp64).libname, "Accelerate")

        # Accelerate has `_rook` symbols:
        @test endswith(BLAS.lbt_find_backing_library("dsytrf_rook_", :ilp64).libname, "Accelerate")
    end
end

@testset "CBLAS dot test" begin
    a = ComplexF64[
        1 + 1im,
        2 - 2im,
        3 + 3im
    ]
    @test BLAS.dotc(a, a) ≈ ComplexF64(28)
    @test BLAS.dotu(a, a) ≈ ComplexF64(12im)

    a = Float32[1, 2, 3]
    @test BLAS.dot(a, a) ≈ 14f0
end

@testset "BLAS threading tests" begin
    if AppleAccelerate.get_macos_version() >= v"15"
        @test AppleAccelerate.set_num_threads(1) == 1
        @test AppleAccelerate.get_num_threads() == 1
        @test AppleAccelerate.set_num_threads(4) > 1
        @test AppleAccelerate.get_num_threads() > 1
    else
        @test AppleAccelerate.get_num_threads() == 1
        @test AppleAccelerate.set_num_threads(1) == 1
    end
end

linalg_stdlib_test_path = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")

@testset verbose=false "LinearAlgebra.jl BLAS tests" begin
    joinpath(linalg_stdlib_test_path, "blas.jl") |> include
end

@testset verbose=false "LinearAlgebra.jl LAPACK tests" begin
    joinpath(linalg_stdlib_test_path, "lapack.jl") |> include
end
