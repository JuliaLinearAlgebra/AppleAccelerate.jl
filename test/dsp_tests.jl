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

    @testset "Float64 input throws (no Float64 DCT in Accelerate)" begin
        # vDSP has no vDSP_DCT_CreateSetupD/vDSP_DCT_ExecuteD; the wrapper must
        # reject Float64 with a clear ArgumentError rather than a MethodError.
        x = rand(Float64, 64)
        @test_throws ArgumentError AppleAccelerate.dct(x)
        @test_throws ArgumentError AppleAccelerate.dct(x, 3)
        @test_throws ArgumentError AppleAccelerate.dct(x, AppleAccelerate.plan_dct(64, 2))
        @test_throws ArgumentError AppleAccelerate.idct(x)
    end

    @testset "idct (inverse DCT via DCT-III)" begin
        n = 2^10
        x = rand(Float32, n)
        # Round-trip: idct is normalized so that idct(dct(x)) ≈ x.
        @test AppleAccelerate.idct(AppleAccelerate.dct(x)) ≈ x rtol=sqrt(eps(Float32))
        # Cross-check against DSP.jl's orthonormal inverse DCT with the fixed
        # vDSP-to-orthonormal scaling (see the DCT testset comment above).
        scale = fill(sqrt(n / 2), n); scale[1] = sqrt(n)
        X_ortho = Float32.(DSP.dct(Float64.(x)))
        @test AppleAccelerate.idct(X_ortho .* Float32.(scale)) ≈ x rtol=1e-3
    end
end

@testset "FFT setup cache" begin
    # Repeated no-plan calls must give identical results and reuse one shared setup.
    for T in (ComplexF32, ComplexF64)
        r = randn(T, 256)
        @test AppleAccelerate.fft(r) == AppleAccelerate.fft(r)
        @test AppleAccelerate.ifft(AppleAccelerate.fft(r)) ≈ r
    end
    # The cache hands back the same FFTSetup object for the same (T, log2n, radix).
    s1 = AppleAccelerate._cached_fftsetup(Float32, 256)
    s2 = AppleAccelerate._cached_fftsetup(Float32, 256)
    @test s1 === s2
    @test s1 isa AppleAccelerate.FFTSetup{Float32}
    @test haskey(AppleAccelerate._FFT_SETUP_CACHE, (Float32, 8, 2))
    @test AppleAccelerate._cached_fftsetup(Float64, 256) !==
          AppleAccelerate._cached_fftsetup(Float64, 512)
    # Real FFT convenience path shares the same cache.
    x = randn(Float64, 128)
    @test AppleAccelerate.rfft(x) == AppleAccelerate.rfft(x)
    @test AppleAccelerate.irfft(AppleAccelerate.rfft(x), 128) ≈ x
    # Mixed-radix (DFT) convenience path is cached too, keyed by (T, n, direction).
    y = randn(ComplexF64, 96)
    @test AppleAccelerate.fft(y) == AppleAccelerate.fft(y)
    @test haskey(AppleAccelerate._DFT_SETUP_CACHE, (Float64, 96, AppleAccelerate.DFT_FORWARD))
    @test AppleAccelerate._cached_dftsetup(Float64, 96, AppleAccelerate.DFT_FORWARD) ===
          AppleAccelerate._cached_dftsetup(Float64, 96, AppleAccelerate.DFT_FORWARD)
    # Explicit plan API still returns fresh, uncached setups.
    @test AppleAccelerate.plan_fft(256, Float32) !== AppleAccelerate.plan_fft(256, Float32)
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

@testset "fft non-power-of-2 (mixed-radix)" begin
    # Supported non-power-of-2 lengths: f * 2^k with f ∈ {3, 5, 15}.
    supported = (3 * 2^4, 3 * 2^5, 5 * 2^4, 5 * 2^5, 15 * 2^3, 15 * 2^4)  # 48,96,80,160,120,240

    @testset "ComplexF64 n=$n" for n in supported
        r = randn(ComplexF64, n)
        ref = FFTW.fft(r)
        @test AppleAccelerate.fft(r) ≈ ref
        # inverse round-trip and match to FFTW's inverse
        @test AppleAccelerate.ifft(AppleAccelerate.fft(r)) ≈ r
        @test AppleAccelerate.ifft(ref) ≈ FFTW.ifft(ref)
        # unnormalized backward transform
        @test AppleAccelerate.bfft(r) ≈ FFTW.bfft(r)
    end

    @testset "ComplexF32 n=$n" for n in supported
        r = randn(ComplexF32, n)
        ref = FFTW.fft(r)
        @test AppleAccelerate.fft(r) ≈ ref rtol=sqrt(eps(Float32))
        @test AppleAccelerate.ifft(AppleAccelerate.fft(r)) ≈ r rtol=sqrt(eps(Float32))
        @test AppleAccelerate.bfft(r) ≈ FFTW.bfft(r) rtol=sqrt(eps(Float32))
    end

    # Power-of-2 fast path is unchanged and still routes through fft().
    @test AppleAccelerate.is_supported_fft_length(1024)
    @test AppleAccelerate.is_supported_fft_length(48)
    @test !AppleAccelerate.is_supported_fft_length(7)
    @test !AppleAccelerate.is_supported_fft_length(11)

    # Unsupported (prime / non-f*2^k) lengths must throw a clear ArgumentError.
    for n in (7, 11, 13, 63)  # 63 = 7*9, odd cofactor 63 ∉ {1,3,5,15}
        @test_throws ArgumentError AppleAccelerate.fft(randn(ComplexF64, n))
        @test_throws ArgumentError AppleAccelerate.fft(randn(ComplexF32, n))
        @test_throws ArgumentError AppleAccelerate.ifft(randn(ComplexF64, n))
        @test_throws ArgumentError AppleAccelerate.bfft(randn(ComplexF64, n))
    end
end

@testset "is_supported_fft_length exact" begin
    # Powers of two are always supported. The only supported non-power-of-2
    # lengths are f*2^k with f ∈ {3, 5, 15} and k ≥ 3 (smallest: 24, 40, 120).
    for n in (3, 5, 6, 10, 12, 15, 20, 30, 60)
        @test !AppleAccelerate.is_supported_fft_length(n)
    end
    for n in (1, 2, 4, 8, 16, 24, 40, 48, 120)
        @test AppleAccelerate.is_supported_fft_length(n)
    end

    # Lengths that route to the mixed-radix DFT with k < 3 must be rejected.
    @test_throws ArgumentError AppleAccelerate.fft(randn(ComplexF64, 12))
    @test_throws ArgumentError AppleAccelerate.fft(randn(ComplexF64, 20))

    # The smallest supported non-power-of-2 length succeeds and matches FFTW.
    let r = randn(ComplexF64, 24)
        @test AppleAccelerate.fft(r) ≈ FFTW.fft(r)
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

@testset "rfft 2D" begin
    for T in (Float64, Float32)
        @testset "$T" begin
            # Non-square sizes (16×32 and 32×16) catch dimension swaps.
            for sz in ((8, 8), (16, 32), (32, 16), (64, 64))
                x = randn(T, sz...)
                result = AppleAccelerate.rfft(x)
                expected = FFTW.rfft(x)
                @test size(result) == (sz[1] ÷ 2 + 1, sz[2])
                if T == Float64
                    @test result ≈ expected
                else
                    @test result ≈ expected rtol=sqrt(eps(Float32))
                end
            end
        end
    end

    @testset "reusable plan" begin
        x = randn(Float64, 16, 32)
        setup = AppleAccelerate.plan_rfft(x)
        @test AppleAccelerate.rfft(x, setup) ≈ FFTW.rfft(x)
        @test AppleAccelerate.rfft(x, setup) ≈ FFTW.rfft(x)  # reuse
    end

    # Non-power-of-2 dimensions are rejected.
    @test_throws AssertionError AppleAccelerate.rfft(randn(Float64, 8, 12))
    @test_throws AssertionError AppleAccelerate.rfft(randn(Float64, 12, 8))
end

@testset "brfft and irfft 2D roundtrip" begin
    for T in (Float64, Float32)
        @testset "$T" begin
            for sz in ((8, 8), (16, 32), (32, 16), (64, 64))
                x = randn(T, sz...)
                X = AppleAccelerate.rfft(x)
                n1 = sz[1]
                if T == Float64
                    @test AppleAccelerate.brfft(X, n1) ./ length(x) ≈ x
                    @test AppleAccelerate.irfft(X, n1) ≈ x
                else
                    @test AppleAccelerate.brfft(X, n1) ./ length(x) ≈ x rtol=sqrt(eps(Float32))
                    @test AppleAccelerate.irfft(X, n1) ≈ x rtol=sqrt(eps(Float32))
                end
            end
        end
    end

    @testset "reusable plan" begin
        x = randn(Float64, 32, 16)
        setup = AppleAccelerate.plan_rfft(x)
        X = AppleAccelerate.rfft(x, setup)
        @test AppleAccelerate.irfft(X, 32, setup) ≈ x
    end
end

@testset "rfft non-power-of-2 (mixed-radix)" begin
    # Lengths f * 2^k with f ∈ {3, 5, 15} supported by vDSP_DFT_zrop.
    supported = (48, 96, 160, 240, 480)
    for T in (Float64, Float32)
        @testset "$T n=$n" for n in supported
            x = randn(T, n)
            result = AppleAccelerate.rfft(x)
            expected = FFTW.rfft(x)
            if T == Float64
                @test result ≈ expected
                @test AppleAccelerate.brfft(result, n) ./ n ≈ x
                @test AppleAccelerate.irfft(result, n) ≈ x
            else
                @test result ≈ expected rtol=sqrt(eps(Float32))
                @test AppleAccelerate.brfft(result, n) ./ n ≈ x rtol=sqrt(eps(Float32))
                @test AppleAccelerate.irfft(result, n) ≈ x rtol=sqrt(eps(Float32))
            end
        end
    end

    @testset "unsupported lengths throw" begin
        for n in (100, 63, 7)  # odd cofactor not in {1,3,5,15} or odd length
            @test_throws ArgumentError AppleAccelerate.rfft(randn(Float64, n))
            @test_throws ArgumentError AppleAccelerate.irfft(randn(ComplexF64, n ÷ 2 + 1), n)
        end
    end
end

@testset "batched fft (fftm)" begin
    for T in (ComplexF64, ComplexF32)
        RT = real(T)
        @testset "$T" begin
            for sz in ((16, 8), (32, 5), (8, 64))
                A = randn(T, sz...)
                for dims in (1, 2)
                    ispow2(size(A, dims)) || continue
                    result = AppleAccelerate.fft(A, dims)
                    expected = FFTW.fft(A, dims)
                    if T == ComplexF64
                        @test result ≈ expected
                        @test AppleAccelerate.bfft(result, dims) ./ size(A, dims) ≈ A
                        @test AppleAccelerate.ifft(result, dims) ≈ A
                    else
                        @test result ≈ expected rtol=sqrt(eps(RT))
                        @test AppleAccelerate.bfft(result, dims) ./ size(A, dims) ≈ A rtol=sqrt(eps(RT))
                        @test AppleAccelerate.ifft(result, dims) ≈ A rtol=sqrt(eps(RT))
                    end
                end
            end
        end
    end

    @testset "reusable plan" begin
        A = randn(ComplexF64, 16, 8)
        setup = AppleAccelerate.plan_fft(16, Float64)
        @test AppleAccelerate.fft(A, 1, setup) ≈ FFTW.fft(A, 1)
        @test AppleAccelerate.ifft(FFTW.fft(A, 1), 1, setup) ≈ A
    end

    @test_throws ArgumentError AppleAccelerate.fft(randn(ComplexF64, 16, 8), 3)
end

@testset "does not hijack FFTW planning (#139)" begin
    # Loading AppleAccelerate must not claim AbstractFFTs.plan_fft, which would
    # break generic FFT planning for sizes vDSP cannot handle.

    # The repro from #139: non-power-of-2, 2D.
    x = rand(ComplexF64, 64, 63)
    @test (FFTW.plan_fft(x) * x) ≈ FFTW.fft(x)

    # 3D: previously threw.
    z = rand(ComplexF64, 4, 4, 4)
    @test (FFTW.plan_fft(z) * z) ≈ FFTW.fft(z)

    # Power-of-2 must also route to FFTW now, not vDSP.
    @test FFTW.plan_fft(rand(ComplexF64, 16)) isa FFTW.cFFTWPlan

    # Fails loudly if an AbstractFFTs extension is ever reintroduced.
    @test !any(m -> occursin("AppleAccelerate", string(m.module)),
               methods(FFTW.plan_fft))

    # The native API is unaffected.
    v = rand(ComplexF64, 16)
    @test AppleAccelerate.fft(v) ≈ FFTW.fft(v)
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

    @testset "Partial numelem returns exactly numelem outputs::Float64" begin
        X::Vector{Float64} = randn(10)
        d::Vector{Float64} = zeros(4)
        c::Vector{Float64} = [x%0.5 for x in randn(5)]
        fa = AppleAccelerate.biquadcreate(c, 1)
        numelem = 6
        result = AppleAccelerate.biquad(X, d, numelem, fa)
        @test length(result) == numelem
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

@testset "Validation / error paths" begin
    # conv! result buffer too small (dsp.jl:66)
    @test_throws ErrorException AppleAccelerate.conv!(zeros(3), randn(10), randn(4))
    # xcorr! result buffer too small (dsp.jl:115)
    @test_throws ErrorException AppleAccelerate.xcorr!(zeros(3), randn(10), randn(4))

    # biquadcreate with too few coefficients (dsp.jl:224)
    @test_throws ErrorException AppleAccelerate.biquadcreate(Float64[1.0, 0, 0], 1)
    @test_throws ErrorException AppleAccelerate.biquadcreate(Float64[1.0, 0, 0], 1, Float32)

    # biquad: incomplete delays (dsp.jl:244) and numelem too large (dsp.jl:248)
    bq = AppleAccelerate.biquadcreate(Float64[1.0, 0, 0, 0, 0], 1)
    X = randn(10)
    @test_throws ErrorException AppleAccelerate.biquad(X, zeros(2), length(X), bq)  # delays too short (need 4)
    @test_throws ErrorException AppleAccelerate.biquad(X, zeros(4), length(X) + 1, bq)  # numelem > length

    bq32 = AppleAccelerate.biquadcreate(Float64[1.0, 0, 0, 0, 0], 1, Float32)
    X32 = randn(Float32, 10)
    @test_throws ErrorException AppleAccelerate.biquad(X32, zeros(Float32, 2), length(X32), bq32)
    @test_throws ErrorException AppleAccelerate.biquad(X32, zeros(Float32, 4), length(X32) + 1, bq32)

    # biquadm_create with too few coefficients (dsp.jl:331)
    @test_throws ErrorException AppleAccelerate.biquadm_create(Float64[1.0, 0, 0, 0, 0], 2, 1, Float32)
    @test_throws ErrorException AppleAccelerate.biquadm_create(Float64[1.0, 0, 0, 0, 0], 2, 1, Float64)

    # plan_dct invalid type (dsp.jl:777) and invalid length (dsp.jl:779)
    @test_throws ErrorException AppleAccelerate.plan_dct(4096, 1)   # unsupported DCT type
    @test_throws ErrorException AppleAccelerate.plan_dct(4096, 5)   # unsupported DCT type
    @test_throws ErrorException AppleAccelerate.plan_dct(7, 2)      # invalid length (f*2^n, n>=4)
    @test_throws ErrorException AppleAccelerate.plan_dct(4097, 2)   # not f*2^n

    # zaspec! output too short already covered; add deq22/wiener validation paths
    @test_throws ErrorException AppleAccelerate.deq22!(zeros(5), randn(5), randn(4))   # B != 5
    @test_throws ErrorException AppleAccelerate.deq22!(zeros(5), randn(2), Float64[1,2,3,4,5])  # A too short
    @test_throws ErrorException AppleAccelerate.deq22!(zeros(6), randn(5), Float64[1,2,3,4,5])  # C != A length
end

@testset "DCT types 2/3/4 vs DSP.dct/idct" begin
    # vDSP DCT-II is unnormalized vs DSP's orthonormal DCT. Cross-validate the
    # supported types (2, 3, 4) against the orthonormal reference with the fixed
    # data-independent scale factors, and exercise setup reuse.
    n = 2^12
    r = Float32.(rand(n))
    rf = Float64.(r)

    # Type II: out[1]=sqrt(N)*ortho[1], out[k]=sqrt(N/2)*ortho[k]
    setup2 = AppleAccelerate.plan_dct(n, 2)
    got2a = AppleAccelerate.dct(r, setup2)
    got2b = AppleAccelerate.dct(r, setup2)   # reuse same setup
    @test got2a == got2b
    ortho2 = DSP.dct(rf)
    scale2 = fill(sqrt(n / 2), n); scale2[1] = sqrt(n)
    @test Float64.(got2a) ≈ ortho2 .* scale2 rtol = 1e-3

    # Type III is the inverse of Type II. vDSP DCT-III applied to the orthonormal
    # coefficients (unscaled here) relates to DSP.idct. Validate round-trip:
    # idct-like recovery via scaled type-III of the forward type-II output.
    setup3 = AppleAccelerate.plan_dct(n, 3)
    # vDSP: DCT-III(DCT-II(x)) = (N/2) * x. Verify the scaling relationship.
    back = AppleAccelerate.dct(got2a, setup3)
    @test Float64.(back) ./ (n / 2) ≈ rf rtol = 1e-3

    # Type IV is its own inverse up to scale: DCT-IV(DCT-IV(x)) = (N/2) * x.
    setup4 = AppleAccelerate.plan_dct(n, 4)
    y4 = AppleAccelerate.dct(r, setup4)
    back4 = AppleAccelerate.dct(y4, setup4)
    @test Float64.(back4) ./ (n / 2) ≈ rf rtol = 1e-3
end

# ---------------------------------------------------------------------
# Additional FFT/DFT variants (temp-buffer, in-place-real, batched, small
# radix, fixed-size, interleaved DFT) and biquad live-state controls.
# ---------------------------------------------------------------------

@testset "temp-buffer FFT variants (zopt/zipt)" begin
    @testset "$T" for T in (Float64, Float32)
        rt = sqrt(eps(T)); CT = Complex{T}
        # 1D
        x = randn(CT, 128)
        s = AppleAccelerate.plan_fft(128, T)
        ws = AppleAccelerate.FFTWorkspace(x)
        @test AppleAccelerate.fft(x, s, ws) ≈ FFTW.fft(x) rtol=rt
        @test AppleAccelerate.ifft(AppleAccelerate.fft(x, s, ws), s, ws) ≈ x rtol=rt
        xc = copy(x); AppleAccelerate.fft!(xc, s, ws); @test xc ≈ FFTW.fft(x) rtol=rt
        xc2 = copy(x); AppleAccelerate.bfft!(xc2, s, ws); @test xc2 ≈ FFTW.bfft(x) rtol=rt
        xc3 = copy(x); AppleAccelerate.ifft!(AppleAccelerate.fft!(xc3, s, ws), s, ws)
        @test xc3 ≈ x rtol=rt
        # 2D (non-square)
        A = randn(CT, 8, 32)
        s2 = AppleAccelerate.plan_fft(32, T); ws2 = AppleAccelerate.FFTWorkspace(A)
        @test AppleAccelerate.fft(A, s2, ws2) ≈ FFTW.fft(A) rtol=rt
        Ac = copy(A); AppleAccelerate.fft!(Ac, s2, ws2); @test Ac ≈ FFTW.fft(A) rtol=rt
        @test AppleAccelerate.ifft(AppleAccelerate.fft(A, s2, ws2), s2, ws2) ≈ A rtol=rt
        # too-small workspace is rejected
        @test_throws DimensionMismatch AppleAccelerate.fft(x, s, AppleAccelerate.FFTWorkspace{T}(4))
    end
end

@testset "temp-buffer & in-place real FFT (zropt/zrip/zript)" begin
    @testset "$T" for T in (Float64, Float32)
        rt = sqrt(eps(T))
        # 1D real, out-of-place + temp
        x = randn(T, 128); s = AppleAccelerate.plan_rfft(x)
        ws = AppleAccelerate.FFTWorkspace(x)
        @test AppleAccelerate.rfft(x, s, ws) ≈ FFTW.rfft(x) rtol=rt
        X = AppleAccelerate.rfft(x, s, ws)
        @test AppleAccelerate.brfft(X, 128, s, ws) ≈ FFTW.brfft(X, 128) rtol=rt
        @test AppleAccelerate.irfft(X, 128, s, ws) ≈ x rtol=rt
        # 1D real, in-place (zrip / zript). rfft! consumes its input.
        @test AppleAccelerate.rfft!(copy(x), s) ≈ FFTW.rfft(x) rtol=rt
        @test AppleAccelerate.rfft!(copy(x), s, ws) ≈ FFTW.rfft(x) rtol=rt
        @test AppleAccelerate.rfft!(copy(x)) ≈ FFTW.rfft(x) rtol=rt
        # 2D real (non-square)
        A = randn(T, 8, 16); s2 = AppleAccelerate.plan_rfft(A)
        ws2 = AppleAccelerate.FFTWorkspace(A)
        @test AppleAccelerate.rfft(A, s2, ws2) ≈ FFTW.rfft(A) rtol=rt
        RA = AppleAccelerate.rfft(A, s2, ws2)
        @test AppleAccelerate.irfft(RA, 8, s2, ws2) ≈ A rtol=rt
        @test AppleAccelerate.brfft(RA, 8, s2, ws2) ≈ FFTW.brfft(RA, 8) rtol=rt
        @test AppleAccelerate.rfft!(copy(A), s2, ws2) ≈ FFTW.rfft(A) rtol=rt
    end
end

@testset "batched FFT extra variants (fftm zip/zipt/zopt/zr*)" begin
    @testset "$T" for T in (Float64, Float32)
        rt = sqrt(eps(T)); CT = Complex{T}
        for sz in ((16, 5), (8, 64)), dims in (1, 2)
            n = size(zeros(sz...), dims); ispow2(n) || continue
            A = randn(CT, sz...)
            s = AppleAccelerate.plan_fft(n, T); ws = AppleAccelerate.FFTWorkspace(A)
            # in-place (zip), temp in-place (zipt), temp out-of-place (zopt)
            Ac = copy(A); AppleAccelerate.fft!(Ac, dims, s); @test Ac ≈ FFTW.fft(A, dims) rtol=rt
            Ac = copy(A); AppleAccelerate.fft!(Ac, dims, s, ws); @test Ac ≈ FFTW.fft(A, dims) rtol=rt
            Ac = copy(A); AppleAccelerate.fft!(Ac, dims); @test Ac ≈ FFTW.fft(A, dims) rtol=rt
            @test AppleAccelerate.fft(A, dims, s, ws) ≈ FFTW.fft(A, dims) rtol=rt
            @test AppleAccelerate.ifft(AppleAccelerate.fft(A, dims, s, ws), dims, s, ws) ≈ A rtol=rt
            Ac = copy(A); AppleAccelerate.ifft!(AppleAccelerate.fft!(Ac, dims), dims); @test Ac ≈ A rtol=rt
        end
        # batched real (zrop/zropt/zrip/zript) along both dims
        for sz in ((16, 5), (8, 7)), dims in (1, 2)
            n = size(zeros(sz...), dims); ispow2(n) || continue
            A = randn(T, sz...)
            s = AppleAccelerate.plan_rfft(A); ws = AppleAccelerate.FFTWorkspace{T}(length(A))
            R = AppleAccelerate.rfft(A, dims)
            @test R ≈ FFTW.rfft(A, dims) rtol=rt
            @test AppleAccelerate.rfft(A, dims, s, ws) ≈ FFTW.rfft(A, dims) rtol=rt
            @test AppleAccelerate.rfft!(copy(A), dims, s) ≈ FFTW.rfft(A, dims) rtol=rt
            @test AppleAccelerate.rfft!(copy(A), dims, s, ws) ≈ FFTW.rfft(A, dims) rtol=rt
            @test AppleAccelerate.irfft(R, n, dims) ≈ A rtol=rt
            @test AppleAccelerate.irfft(R, n, dims, s, ws) ≈ A rtol=rt
            @test AppleAccelerate.brfft(R, n, dims) ≈ FFTW.brfft(R, n, dims) rtol=rt
        end
        @test_throws ArgumentError AppleAccelerate.fft!(randn(CT, 16, 8), 3, AppleAccelerate.plan_fft(16, T))
    end
end

@testset "small-radix FFT (fft3/fft5)" begin
    @testset "$T" for T in (Float64, Float32)
        rt = sqrt(eps(T)); CT = Complex{T}
        for n in (3, 6, 12, 24, 48)
            x = randn(CT, n)
            @test AppleAccelerate.fftradix3(x) ≈ FFTW.fft(x) rtol=rt
            @test AppleAccelerate.bfftradix3(AppleAccelerate.fftradix3(x)) ≈ n .* x rtol=rt
        end
        for n in (5, 10, 20, 40)
            x = randn(CT, n)
            @test AppleAccelerate.fftradix5(x) ≈ FFTW.fft(x) rtol=rt
            @test AppleAccelerate.bfftradix5(AppleAccelerate.fftradix5(x)) ≈ n .* x rtol=rt
        end
        # radix-3 handles length 12 and radix-5 handles length 20, which the
        # mixed-radix DFT path cannot (verified above); guards reject bad lengths.
        @test_throws ArgumentError AppleAccelerate.fftradix3(randn(CT, 7))
        @test_throws ArgumentError AppleAccelerate.fftradix5(randn(CT, 7))
    end
end

@testset "fixed-size FFT16/FFT32 (Float32)" begin
    rt = sqrt(eps(Float32))
    for (N, f, bf) in ((16, AppleAccelerate.fft16, AppleAccelerate.bfft16),
                       (32, AppleAccelerate.fft32, AppleAccelerate.bfft32))
        x = randn(ComplexF32, N)
        @test f(x) ≈ FFTW.fft(x) rtol=rt
        @test bf(f(x)) ≈ N .* x rtol=rt
    end
    @test_throws DimensionMismatch AppleAccelerate.fft16(randn(ComplexF32, 8))
    @test_throws DimensionMismatch AppleAccelerate.fft32(randn(ComplexF32, 16))
end

@testset "interleaved DFT (vDSP_DFT_Interleaved)" begin
    @testset "$T" for T in (Float64, Float32)
        rt = sqrt(eps(T)); CT = Complex{T}
        for n in (16, 24, 32, 48)
            x = randn(CT, n)
            s = AppleAccelerate.plan_dft_interleaved(n, AppleAccelerate.DFT_FORWARD, T)
            @test AppleAccelerate.dft_interleaved(x, s) ≈ FFTW.fft(x) rtol=rt
            @test AppleAccelerate.dft_interleaved(x) ≈ FFTW.fft(x) rtol=rt
            @test AppleAccelerate.idft_interleaved(AppleAccelerate.dft_interleaved(x)) ≈ x rtol=rt
        end
        @test_throws ArgumentError AppleAccelerate.plan_dft_interleaved(5, AppleAccelerate.DFT_FORWARD, T)
    end
end

@testset "biquad single-channel coefficient setter" begin
    # Only the single-precision setup supports live coefficient updates.
    bq = AppleAccelerate.biquadcreate([1.0, 0, 0, 0, 0], 1, Float32)  # passthrough
    x = Float32.(1:8)
    @test AppleAccelerate.biquad(x, zeros(Float32, 4), 8, bq) ≈ x
    AppleAccelerate.biquad_setcoefficients!(bq, Float32[2, 0, 0, 0, 0])   # Single setter
    @test AppleAccelerate.biquad(x, zeros(Float32, 4), 8, bq) ≈ 2 .* x
    AppleAccelerate.biquad_setcoefficients!(bq, [3.0, 0, 0, 0, 0])        # Double setter
    @test AppleAccelerate.biquad(x, zeros(Float32, 4), 8, bq) ≈ 3 .* x
    # Float64 setups have no coefficient setter in vDSP.
    @test_throws ArgumentError AppleAccelerate.biquad_setcoefficients!(
        AppleAccelerate.biquadcreate([1.0, 0, 0, 0, 0], 1, Float64), [1.0, 0, 0, 0, 0])
    @test_throws DimensionMismatch AppleAccelerate.biquad_setcoefficients!(bq, Float32[1, 0, 0], 0, 1)
end

@testset "biquadm live-state controls" begin
    @testset "$T" for T in (Float32, Float64)
        x = T.(1:8)
        # SetCoefficients: retune one channel live (2 channels, 1 section).
        setup = AppleAccelerate.biquadm_create([1.0,0,0,0,0, 1.0,0,0,0,0], 2, 1, T)
        Y = AppleAccelerate.biquadm([copy(x), copy(x)], 8, setup)
        @test Y[1] ≈ x && Y[2] ≈ x
        AppleAccelerate.biquadm_setcoefficients!(setup, T[5,0,0,0,0], 0, 1, 1, 1)  # channel 1 -> 5x
        Y2 = AppleAccelerate.biquadm([copy(x), copy(x)], 8, setup)
        @test Y2[1] ≈ x
        @test Y2[2] ≈ 5 .* x

        # ResetState: a filter with feedback carries state across calls; reset clears it.
        cr = [1.0, 0, 0, -0.5, 0]     # b0=1, a1=-0.5 (has memory)
        sr = AppleAccelerate.biquadm_create(cr, 1, 1, T)
        imp = T[1, 0, 0, 0, 0, 0, 0, 0]
        ya = AppleAccelerate.biquadm([copy(imp)], 8, sr)
        yb = AppleAccelerate.biquadm([copy(imp)], 8, sr)   # continues from leftover state
        AppleAccelerate.biquadm_resetstate!(sr)
        yc = AppleAccelerate.biquadm([copy(imp)], 8, sr)   # fresh again
        @test ya[1] ≈ yc[1]
        @test !(yb[1] ≈ ya[1])

        # CopyState: transplant state so the next block matches.
        s1 = AppleAccelerate.biquadm_create(cr, 1, 1, T)
        s2 = AppleAccelerate.biquadm_create(cr, 1, 1, T)
        AppleAccelerate.biquadm([copy(imp)], 8, s1)        # s1 accumulates state
        AppleAccelerate.biquadm_copystate!(s2, s1)
        nxt = T[1, 0, 0, 0, 0, 0, 0, 0]
        @test AppleAccelerate.biquadm([copy(nxt)], 8, s1)[1] ≈
              AppleAccelerate.biquadm([copy(nxt)], 8, s2)[1]

        # SetActiveFilters: disabling a section bypasses it (2 sections, each ×2).
        sa = AppleAccelerate.biquadm_create([2.0,0,0,0,0, 2.0,0,0,0,0], 1, 2, T)
        @test AppleAccelerate.biquadm([copy(x)], 8, sa)[1] ≈ 4 .* x
        AppleAccelerate.biquadm_setactivefilters!(sa, Bool[true, false])
        AppleAccelerate.biquadm_resetstate!(sa)
        @test AppleAccelerate.biquadm([copy(x)], 8, sa)[1] ≈ 2 .* x

        # SetTargets: drives coefficients toward a target; steady state reaches it.
        st = AppleAccelerate.biquadm_create([1.0,0,0,0,0], 1, 1, T)
        AppleAccelerate.biquadm_settargets!(st, T[4,0,0,0,0], 0.5, 0.0, 0, 0, 1, 1)
        yt = AppleAccelerate.biquadm([ones(T, 200)], 200, st)
        @test yt[1][end] ≈ 4 rtol=sqrt(eps(T))

        # Guards.
        @test_throws DimensionMismatch AppleAccelerate.biquadm_setcoefficients!(setup, T[1,0,0], 0, 0, 1, 1)
        @test_throws ArgumentError AppleAccelerate.biquadm_setcoefficients!(setup, T[1,0,0,0,0], 0, 5, 1, 1)
        @test_throws DimensionMismatch AppleAccelerate.biquadm_setactivefilters!(sa, Bool[true])  # too short
        @test_throws ArgumentError AppleAccelerate.biquadm_copystate!(
            AppleAccelerate.biquadm_create(cr, 1, 1, T),
            AppleAccelerate.biquadm_create([1.0,0,0,0,0, 1.0,0,0,0,0], 2, 1, T))
    end
end

end # @testset "Signal Processing"
