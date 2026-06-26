@testset "Signal Processing" begin

@testset "DCT" begin
    # vDSP's DCT-II is unnormalized relative to DSP's orthonormal DCT-II by a
    # *fixed* (data-independent) factor: out[1] = sqrt(N)*ortho[1] and
    # out[k] = sqrt(N/2)*ortho[k] for k > 1. Compare against that fixed scaling
    # directly (not a data-derived ratio, which would mask scale-factor errors).
    @testset "$T" for T in (Float32, Float64)
        n = 2^12
        r = rand(T, n)
        d_ortho = DSP.dct(Float64.(r))               # orthonormal reference
        scale = fill(sqrt(n / 2), n); scale[1] = sqrt(n)
        ref = Float64.(d_ortho) .* scale
        # AppleAccelerate.dct only takes Float32; cast the input for both paths.
        got = AppleAccelerate.dct(Float32.(r))
        @test Float64.(got) ≈ ref rtol = (T == Float32 ? 1e-3 : 1e-4)
    end
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

    @testset "wrong-size input to a plan throws" begin
        # A plan built for size N must reject a differently-sized input rather
        # than silently returning a wrong-size result.
        p = plan_fft(randn(ComplexF64, 64))
        @test_throws DimensionMismatch p * randn(ComplexF64, 16)
        y = Vector{ComplexF64}(undef, 16)
        @test_throws DimensionMismatch mul!(y, p, randn(ComplexF64, 16))

        # In-place plan: same guard.
        pin = plan_fft!(randn(ComplexF64, 64))
        @test_throws DimensionMismatch pin * randn(ComplexF64, 16)

        # 2D plan rejects wrong matrix size.
        p2 = plan_fft(randn(ComplexF64, 16, 16))
        @test_throws DimensionMismatch p2 * randn(ComplexF64, 8, 8)
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

        # Regression: an oversized `result` buffer must not cause vDSP to read
        # past the zero-padded input (the padding is sized from length(result)).
        @testset "oversized result buffer (no OOB)::$T" begin
            K = randn(T, 5)
            natural = AppleAccelerate.conv(X, K)         # length N + 5 - 1
            extra = 7
            big = fill(T(NaN), length(natural) + extra)  # deliberately too long
            AppleAccelerate.conv!(big, X, K)
            # Valid region must match the natural-length convolution; the result
            # is well-defined for all length(result) elements (the extra tail
            # reads from the in-bounds zero padding, so it must not be NaN).
            @test big[1:length(natural)] ≈ natural
            @test !any(isnan, big)

            # xcorr! with an oversized buffer likewise stays in-bounds.
            bigx = fill(T(NaN), length(natural) + extra)
            AppleAccelerate.xcorr!(bigx, X, Y[1:5])
            @test bigx[1:length(natural)] ≈ AppleAccelerate.xcorr(X, Y[1:5])
            @test !any(isnan, bigx)
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

@testset "Multi-channel Biquad layout (sections × channels)" begin
    # Distinct per-channel gains with 2 cascaded sections verify both the
    # (sections, channels) setup argument order and the section-major
    # coefficient layout. With sections=1 these would be indistinguishable.
    #   channel 0: two passthrough sections        -> y = x
    #   channel 1: gains of 2 then 3 (cascaded)     -> y = 6x
    x = Float32.(collect(1:8))
    # section-major: [sec0_ch0, sec0_ch1, sec1_ch0, sec1_ch1], each [b0,b1,b2,a1,a2]
    c = [1.0,0,0,0,0,  2.0,0,0,0,0,     # section 0: ch0 passthrough, ch1 ×2
         1.0,0,0,0,0,  3.0,0,0,0,0]     # section 1: ch0 passthrough, ch1 ×3
    setup = AppleAccelerate.biquadm_create(c, 2, 2, Float32)
    Y = AppleAccelerate.biquadm([copy(x), copy(x)], 8, setup)
    @test Y[1] ≈ x
    @test Y[2] ≈ 6 .* x
    # numelem larger than the channel length must be rejected, not read OOB
    @test_throws ErrorException AppleAccelerate.biquadm([copy(x), copy(x)], 9, setup)
end

@testset "Multi-channel Biquad layout::Float64" begin
    # Float64 counterpart of the multichannel test above.
    x64 = Float64.(collect(1:8))
    c64 = [1.0,0,0,0,0,  2.0,0,0,0,0,
           1.0,0,0,0,0,  3.0,0,0,0,0]
    setup64 = AppleAccelerate.biquadm_create(c64, 2, 2, Float64)
    Y64 = AppleAccelerate.biquadm([copy(x64), copy(x64)], 8, setup64)
    @test Y64[1] ≈ x64
    @test Y64[2] ≈ 6 .* x64
    # numelem larger than the channel length must be rejected, not read OOB
    @test_throws ErrorException AppleAccelerate.biquadm([copy(x64), copy(x64)], 9, setup64)
end

@testset "deq22 (recursive filter)" begin
    for T in (Float32, Float64)
        @testset "$T" begin
            N = 64
            A = randn(T, N)
            B = T[0.5, -0.3, 0.2, 0.1, -0.05]

            # Test allocating version
            C = AppleAccelerate.deq22(A, B)
            @test length(C) == N

            # Reference: manually compute difference equation
            # Allocating version pads A with 2 leading zeros
            Apad = [T(0); T(0); A]
            Cref = zeros(T, N + 2)  # C[1:2] = 0 (initial state)
            for n in 3:(N+2)
                Cref[n] = Apad[n]*B[1] + Apad[n-1]*B[2] + Apad[n-2]*B[3] - Cref[n-1]*B[4] - Cref[n-2]*B[5]
            end
            @test C ≈ Cref[3:end]

            # Test mutating version
            Apad2 = [T(0); T(0); A]
            C2 = zeros(T, N + 2)
            AppleAccelerate.deq22!(C2, Apad2, B)
            @test C2[3:end] ≈ Cref[3:end]
        end
    end
end

@testset "desamp (FIR decimation)" begin
    for T in (Float32, Float64)
        @testset "$T" begin
            A = randn(T, 100)
            F = randn(T, 5)
            DF = 3

            C = AppleAccelerate.desamp(A, DF, F)
            P = length(F)
            Nout = div(length(A) - P, DF) + 1
            @test length(C) == Nout

            # Reference: compute via dot products
            Cref = Vector{T}(undef, Nout)
            for n in 0:(Nout-1)
                s = T(0)
                for p in 0:(P-1)
                    s += A[n*DF + p + 1] * F[p + 1]
                end
                Cref[n+1] = s
            end
            @test C ≈ Cref

            # Test DF=1 (pure FIR, no decimation)
            C1 = AppleAccelerate.desamp(A, 1, F)
            Nout1 = div(length(A) - P, 1) + 1
            Cref1 = Vector{T}(undef, Nout1)
            for n in 0:(Nout1-1)
                s = T(0)
                for p in 0:(P-1)
                    s += A[n + p + 1] * F[p + 1]
                end
                Cref1[n+1] = s
            end
            @test C1 ≈ Cref1

            # Test mutating version
            Cout = Vector{T}(undef, Nout)
            AppleAccelerate.desamp!(Cout, A, DF, F)
            @test Cout ≈ Cref
        end
    end

    # Regression: P = UInt64(length(F)); length(A) - P used to underflow when
    # the filter is longer than the signal, instead of erroring clearly.
    @testset "desamp errors when filter longer than signal::$T" for T in (Float32, Float64)
        Ashort = randn(T, 3)
        Flong = randn(T, 5)
        @test_throws ErrorException AppleAccelerate.desamp(Ashort, 1, Flong)
        Cout = Vector{T}(undef, 1)
        @test_throws ErrorException AppleAccelerate.desamp!(Cout, Ashort, 1, Flong)
    end
end

@testset "wiener (Wiener-Levinson)" begin
    for T in (Float32, Float64)
        @testset "$T" begin
            L = 8
            # Build a valid positive-definite Toeplitz autocorrelation
            # Use a simple decaying autocorrelation
            autocorr = T[T(1) / (1 + abs(T(k))) for k in 0:(L-1)]

            # Arbitrary cross-correlation
            crosscorr = randn(T, L)

            F, err = AppleAccelerate.wiener(autocorr, crosscorr)
            @test err == 0
            @test length(F) == L

            # Verify Wiener-Hopf equation: R * F ≈ crosscorr
            # R is the Toeplitz matrix formed from autocorr
            R = zeros(T, L, L)
            for i in 1:L, j in 1:L
                R[i,j] = autocorr[abs(i-j)+1]
            end
            @test R * F ≈ crosscorr rtol=(T == Float32 ? sqrt(eps(Float32)) : 1e-10)

            # Test mutating version
            F2 = Vector{T}(undef, L)
            P2 = Vector{T}(undef, L)
            F2r, err2 = AppleAccelerate.wiener!(F2, P2, autocorr, crosscorr)
            @test err2 == 0
            @test F2r ≈ F
        end
    end
end

@testset "DFT (complex-to-complex)" begin
    for T in (Float32, Float64)
        CT = Complex{T}
        @testset "plan_dft and dft::$T" begin
            # Power-of-2 size: compare against FFT
            n = 64
            x = randn(CT, n)
            setup_fwd = AppleAccelerate.plan_dft(n, AppleAccelerate.DFT_FORWARD, T)
            X_dft = AppleAccelerate.dft(x, setup_fwd)
            X_fft = FFTW.fft(x)
            if T == Float64
                @test X_dft ≈ X_fft
            else
                @test X_dft ≈ X_fft rtol=sqrt(eps(Float32))
            end
        end

        @testset "dft roundtrip::$T" begin
            # Non-power-of-2 sizes: f * 2^n where f ∈ {1,3,5,15}, n ≥ 3
            for n in (24, 40, 120)
                x = randn(CT, n)
                X = AppleAccelerate.dft(x)
                x_recovered = AppleAccelerate.idft(X)
                if T == Float64
                    @test x_recovered ≈ x
                else
                    @test x_recovered ≈ x rtol=sqrt(eps(Float32))
                end
            end
        end

        @testset "dft split-complex::$T" begin
            n = 32
            Ir = randn(T, n)
            Ii = randn(T, n)
            setup = AppleAccelerate.plan_dft(n, AppleAccelerate.DFT_FORWARD, T)
            Or, Oi = AppleAccelerate.dft(Ir, Ii, setup)
            # Compare with interleaved version
            X = complex.(Ir, Ii)
            Y = AppleAccelerate.dft(X, setup)
            @test complex.(Or, Oi) ≈ Y
        end

        @testset "idft with setup::$T" begin
            n = 48  # 3 * 2^4
            x = randn(CT, n)
            setup_fwd = AppleAccelerate.plan_dft(n, AppleAccelerate.DFT_FORWARD, T)
            setup_inv = AppleAccelerate.plan_dft(n, AppleAccelerate.DFT_INVERSE, T)
            X = AppleAccelerate.dft(x, setup_fwd)
            x_rec = AppleAccelerate.idft(X, setup_inv)
            if T == Float64
                @test x_rec ≈ x
            else
                @test x_rec ≈ x rtol=sqrt(eps(Float32))
            end
        end
    end

    @testset "invalid DFT length" begin
        @test_throws ErrorException AppleAccelerate.plan_dft(7, AppleAccelerate.DFT_FORWARD, Float32)
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

        # Regression: a short sibling array must throw DimensionMismatch rather
        # than make vDSP read out of bounds.
        @testset "short-array validation::$T" begin
            n = 64
            A   = CT.(randn(n) .+ im .* randn(n))
            Ar  = abs.(randn(T, n)) .+ T(0.01)
            Br  = abs.(randn(T, n)) .+ T(0.01)
            Cc  = CT.(randn(n) .+ im .* randn(n))
            short_real = zeros(T, n - 1)
            short_cplx = zeros(CT, n - 1)

            @test_throws DimensionMismatch AppleAccelerate.zaspec!(short_real, A)
            @test_throws DimensionMismatch AppleAccelerate.zcoher!(short_real, Ar, Br, Cc)
            @test_throws DimensionMismatch AppleAccelerate.ztrans!(short_cplx, Ar, Cc)
            @test_throws DimensionMismatch AppleAccelerate.zcspec!(short_cplx, A, Cc)
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
