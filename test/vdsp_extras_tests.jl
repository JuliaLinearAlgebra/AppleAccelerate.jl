## vdsp_extras_tests.jl ##
# Tests for src/vdsp_extras.jl: the type-dispatched int ↔ float front end, the
# cross-precision biquadm setters, and the allocation-free fixed-size FFTs.

const _EXTRA_INTS = (Int8, Int16, Int32, UInt8, UInt16, UInt32)

@testset "float_to_int / int_to_float" begin
    @testset "$T → $I" for T in (Float32, Float64), I in _EXTRA_INTS
        # In-range values only (out-of-range results are unspecified by vDSP), chosen
        # to discriminate the rounding mode: exact ties, both signs, near the limits.
        lo = I <: Signed ? T(-100) : T(0)
        x = vcat(T[0.5, 1.5, 2.5, 3.5, 2.7, 2.2, 99.49, 100.5],
                 I <: Signed ? T[-0.5, -1.5, -2.5, -2.7, -2.2, -99.5] : T[],
                 lo .+ (T(100) - lo) .* rand(T, 200))

        @test AppleAccelerate.float_to_int(I, x) == unsafe_trunc.(I, x)
        @test AppleAccelerate.float_to_int(I, x; rounding = :trunc) == unsafe_trunc.(I, x)
        @test AppleAccelerate.float_to_int(I, x; rounding = :nearest) == round.(I, x)   # ties → even
        @test AppleAccelerate.float_to_int(I, x) isa Vector{I}

        # The front end picks exactly the C-named kernel.
        u = I <: Signed ? "" : "u"; w = 8 * sizeof(I)
        @test AppleAccelerate.float_to_int(I, x) ==
              getfield(AppleAccelerate, Symbol("vfix", u, w))(x)
        @test AppleAccelerate.float_to_int(I, x; rounding = :nearest) ==
              getfield(AppleAccelerate, Symbol("vfixr", u, w))(x)

        # Extremes of the integer range that the float type represents exactly.
        ext = T[typemin(I), I === Int32 || I === UInt32 ? T(2)^30 : T(typemax(I))]
        @test AppleAccelerate.float_to_int(I, ext; rounding = :nearest) == I.(ext)

        # Mutating form, strided input and output (every other element, and reversed).
        buf = zeros(I, 2 * length(x))
        C = view(buf, 1:2:length(buf))
        @test AppleAccelerate.float_to_int!(C, x; rounding = :nearest) === C
        @test C == round.(I, x)
        @test all(iszero, view(buf, 2:2:length(buf)))      # untouched between strides
        xs = view(x, length(x):-2:1)
        @test AppleAccelerate.float_to_int(I, xs) == unsafe_trunc.(I, xs)

        # Round trip.
        n = I.(unsafe_trunc.(I, x))
        f = AppleAccelerate.int_to_float(T, n)
        @test f isa Vector{T}
        @test f == T.(n)
        @test AppleAccelerate.float_to_int(I, f) == n

        # int → float over the full integer range, strided.
        full = rand(I, 300)
        push!(full, typemin(I), typemax(I))
        @test AppleAccelerate.int_to_float(T, full) == T.(full)
        fv = view(full, 1:3:length(full))
        out = zeros(T, length(fv))
        @test AppleAccelerate.int_to_float!(out, fv) === out
        @test out == T.(fv)

        # Guards.
        @test_throws DimensionMismatch AppleAccelerate.float_to_int!(zeros(I, 2), x)
        @test_throws DimensionMismatch AppleAccelerate.int_to_float!(zeros(T, 2), full)
        @test_throws ArgumentError AppleAccelerate.float_to_int(I, x; rounding = :up)
        @test isempty(AppleAccelerate.float_to_int(I, T[]))
        @test isempty(AppleAccelerate.int_to_float(T, I[]))
    end

    # Unsupported element types are a MethodError, not a silent reinterpretation.
    @test_throws MethodError AppleAccelerate.float_to_int(Int64, [1.0, 2.0])
    @test_throws MethodError AppleAccelerate.int_to_float(Float64, Int64[1, 2])
end

@testset "biquadm cross-precision setters" begin
    # (setup precision T, coefficient precision S): S ≠ T hits the
    # SetCoefficientsDouble / SetCoefficientsSingleD kernels.
    @testset "setup $T ← coeffs $S" for (T, S) in ((Float32, Float64), (Float64, Float32))
        x = T.(1:8)
        pass = [1.0,0,0,0,0, 1.0,0,0,0,0]                     # 2 channels, 1 section
        setup = AppleAccelerate.biquadm_create(pass, 2, 1, T)
        @test AppleAccelerate.biquadm_setcoefficients!(setup, S[5,0,0,0,0], 0, 1, 1, 1) === setup
        Y = AppleAccelerate.biquadm([copy(x), copy(x)], 8, setup)
        @test Y[1] ≈ x                                        # channel 0 untouched
        @test Y[2] ≈ 5 .* x                                   # channel 1 retuned

        # A filter with memory, against a setup built directly from the same coefficients.
        c = [0.25, 0.5, 0.25, -0.3, 0.1]
        ref = AppleAccelerate.biquadm_create(c, 1, 1, T)
        got = AppleAccelerate.biquadm_create([1.0,0,0,0,0], 1, 1, T)
        AppleAccelerate.biquadm_setcoefficients!(got, S.(c), 0, 0, 1, 1)
        sig = T.(randn(64))
        @test AppleAccelerate.biquadm([copy(sig)], 64, got)[1] ≈
              AppleAccelerate.biquadm([copy(sig)], 64, ref)[1]  rtol = sqrt(eps(Float32))

        # SetTargets: interpolates toward the target; steady state reaches it.
        st = AppleAccelerate.biquadm_create([1.0,0,0,0,0], 1, 1, T)
        @test AppleAccelerate.biquadm_settargets!(st, S[4,0,0,0,0], 0.5, 0.0, 0, 0, 1, 1) === st
        yt = AppleAccelerate.biquadm([ones(T, 200)], 200, st)
        @test yt[1][end] ≈ 4 rtol = sqrt(eps(Float32))

        # Guards (same contract as the matching-precision methods).
        @test_throws DimensionMismatch AppleAccelerate.biquadm_setcoefficients!(setup, S[1,0,0], 0, 0, 1, 1)
        @test_throws ArgumentError AppleAccelerate.biquadm_setcoefficients!(setup, S[1,0,0,0,0], 0, 5, 1, 1)
        @test_throws ArgumentError AppleAccelerate.biquadm_setcoefficients!(setup, S[1,0,0,0,0], 1, 0, 1, 1)
        @test_throws ArgumentError AppleAccelerate.biquadm_setcoefficients!(setup, S[1,0,0,0,0], -1, 0, 1, 1)
        @test_throws DimensionMismatch AppleAccelerate.biquadm_settargets!(setup, S[1,0,0], 0.5, 0.0, 0, 0, 1, 1)
        @test_throws ArgumentError AppleAccelerate.biquadm_settargets!(setup, S[1,0,0,0,0], 0.5, 0.0, 0, 2, 1, 1)
    end
end

@testset "fixed-size FFT, in place (copv)" begin
    @testset "N = $N" for (N, f!, b!, f) in
            ((16, AppleAccelerate.fft16!, AppleAccelerate.bfft16!, AppleAccelerate.fft16),
             (32, AppleAccelerate.fft32!, AppleAccelerate.bfft32!, AppleAccelerate.fft32))
        x = randn(ComplexF32, N)
        ref = FFTW.fft(ComplexF64.(x))
        tol = 10 * sqrt(eps(Float32))

        # Out of place: input preserved, result matches FFTW and the zopv-based wrapper.
        x0 = copy(x); out = zeros(ComplexF32, N)
        @test f!(out, x) === out
        @test x == x0
        @test out ≈ ref rtol = tol
        @test out ≈ f(x)

        # In place.
        y = copy(x)
        @test f!(y) === y
        @test y ≈ ref rtol = tol

        # Unnormalized inverse.
        @test b!(zeros(ComplexF32, N), x) ≈ FFTW.bfft(ComplexF64.(x)) rtol = tol
        @test b!(f!(copy(x))) ≈ N .* x rtol = tol

        # Aligned unit-stride views into a larger buffer work copy-free.
        big = randn(ComplexF32, 4N)
        v = view(big, 2N+1:3N)
        before = copy(big)
        vref = FFTW.fft(ComplexF64.(v))
        f!(v)
        @test v ≈ vref rtol = tol
        @test big[1:2N] == before[1:2N] && big[3N+1:end] == before[3N+1:end]

        # No allocation on the mutating paths.
        f!(out, x); f!(y)
        @test (@allocated f!(out, x)) == 0
        @test (@allocated f!(y)) == 0

        # Guards.
        @test_throws DimensionMismatch f!(zeros(ComplexF32, N - 1))
        @test_throws DimensionMismatch f!(zeros(ComplexF32, N), zeros(ComplexF32, 2N))
        @test_throws ArgumentError f!(view(zeros(ComplexF32, 2N), 1:2:2N))          # non-unit stride
        @test_throws ArgumentError f!(view(big, 2:N+1))                              # misaligned
        @test_throws ArgumentError f!(view(big, 1:N), view(big, 3:N+2))              # partial overlap
    end
end

@testset "VDSP_COVERAGE note" begin
    @test AppleAccelerate.VDSP_COVERAGE === nothing
    doc = string(@doc AppleAccelerate.VDSP_COVERAGE)
    @test occursin("vDSP_fft2d_zrip", doc)
    @test occursin("vDSP_DFT_CreateSetup", doc)
end
