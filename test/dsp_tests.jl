@testset "Signal Processing" begin

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

@testset "Biquadratic Filtering::Float64" begin
    @testset "Single Section::Float64" begin
        X::Vector{Float64} = randn(10)
        d::Vector{Float64} = zeros(4)
        c::Vector{Float64} = [x%0.5 for x in randn(5)]
        fdsp = DSP.Biquad(c[1], c[2], c[3], c[4], c[5])
        fa = AppleAccelerate.biquadcreate(c, 1)
        @test DSP.filt(fdsp, X) ≈ AppleAccelerate.biquad(X, d, length(X), fa)
    end
end

@testset "Biquadratic Filtering::Float32" begin
    @testset "Single Section::Float32" begin
        c = [x%0.5 for x in randn(5)]  # Float64 coefficients (required by Apple API)
        X32 = Float32.(randn(10))
        d32 = zeros(Float32, 4)
        fa32 = AppleAccelerate.biquadcreate(c, 1, Float32)
        result32 = AppleAccelerate.biquad(X32, d32, length(X32), fa32)

        # Compare with Float64 result (allowing for precision loss)
        X64 = Float64.(X32)
        d64 = zeros(Float64, 4)
        fa64 = AppleAccelerate.biquadcreate(c, 1, Float64)
        result64 = AppleAccelerate.biquad(X64, d64, length(X64), fa64)
        @test Float32.(result64) ≈ result32 rtol=sqrt(eps(Float32))
    end
end

@testset "Multi-channel Biquad::Float32" begin
    # Simple lowpass filter coefficients (1 section per channel)
    c = [0.1, 0.2, 0.1, -0.5, 0.2,   # channel 1
         0.1, 0.2, 0.1, -0.5, 0.2]    # channel 2
    channels = 2
    sections = 1
    setup = AppleAccelerate.biquadm_create(c, channels, sections, Float32)
    X = [Float32.(randn(64)), Float32.(randn(64))]
    Y = AppleAccelerate.biquadm(X, 64, setup)
    @test length(Y) == 2
    @test length(Y[1]) == 64
    @test length(Y[2]) == 64
    # Results should not be all zeros (filter was applied)
    @test !all(iszero, Y[1])
    @test !all(iszero, Y[2])
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

@testset "Spectral Analysis" begin
    for T in (Float32, Float64)
        CT = Complex{T}
        @testset "zaspec (autospectrum)::$T" begin
            A = CT.(randn(64) .+ im .* randn(64))
            C = zeros(T, 64)
            AppleAccelerate.zaspec!(C, A)
            @test C ≈ abs2.(A)

            # Accumulating: calling twice should double
            C2 = zeros(T, 64)
            AppleAccelerate.zaspec!(C2, A)
            AppleAccelerate.zaspec!(C2, A)
            @test C2 ≈ 2 .* abs2.(A)

            # Allocating version
            @test AppleAccelerate.zaspec(A) ≈ abs2.(A)
        end

        @testset "zcoher (coherence)::$T" begin
            n = 64
            A = abs.(randn(T, n)) .+ T(0.01)  # power spectrum 1, positive
            B = abs.(randn(T, n)) .+ T(0.01)  # power spectrum 2, positive
            C = CT.(randn(n) .+ im .* randn(n))  # cross-spectrum

            D = AppleAccelerate.zcoher(A, B, C)
            expected = abs2.(C) ./ (A .* B)
            @test D ≈ expected
        end

        @testset "ztrans (transfer function)::$T" begin
            n = 64
            A = abs.(randn(T, n)) .+ T(0.01)  # real power spectrum, positive
            B = CT.(randn(n) .+ im .* randn(n))  # complex cross-spectrum

            C = AppleAccelerate.ztrans(A, B)
            expected = B ./ A
            @test C ≈ expected
        end

        @testset "zcspec (cross-spectrum)::$T" begin
            n = 64
            A = CT.(randn(n) .+ im .* randn(n))
            B = CT.(randn(n) .+ im .* randn(n))

            C = AppleAccelerate.zcspec(A, B)
            expected = conj.(A) .* B
            @test C ≈ expected

            # Accumulating: calling twice
            C2 = zeros(CT, n)
            AppleAccelerate.zcspec!(C2, A, B)
            AppleAccelerate.zcspec!(C2, A, B)
            @test C2 ≈ 2 .* (conj.(A) .* B)
        end
    end
end

@testset "plan_fft region argument" begin
    @testset "2D FFT along dim 1 only" begin
        x = randn(ComplexF64, 16, 8)
        p = plan_fft(x, (1,))
        result = p * x
        expected = FFTW.fft(x, (1,))
        @test result ≈ expected
    end

    @testset "2D FFT along dim 2 only" begin
        x = randn(ComplexF64, 16, 8)
        p = plan_fft(x, (2,))
        result = p * x
        expected = FFTW.fft(x, (2,))
        @test result ≈ expected
    end

    @testset "2D FFT along both dims matches full 2D FFT" begin
        x = randn(ComplexF64, 8, 16)
        p12 = plan_fft(x, 1:2)
        p_default = plan_fft(x)
        @test (p12 * x) ≈ (p_default * x)
    end

    @testset "2D ifft roundtrip with region" begin
        x = randn(ComplexF64, 8, 8)
        for region in ((1,), (2,), 1:2)
            fwd = plan_fft(x, region) * x
            result = plan_bfft(x, region) * fwd
            n = prod(size(x, d) for d in region)
            @test result ./ n ≈ x
        end
    end

    @testset "2D plan_fft region Float32" begin
        x = randn(ComplexF32, 8, 8)
        for region in ((1,), (2,))
            result = plan_fft(x, region) * x
            expected = FFTW.fft(x, collect(region))
            @test result ≈ expected rtol=sqrt(eps(Float32))
        end
    end
end

end # @testset "Signal Processing"
