# The array-type contract (see docs/src/array.md, "Array types and strides"):
#
#   * routines whose C prototype takes a stride accept any `StridedVector` and pass the
#     real stride through — a `view` gives the same answer as the `Vector` it selects;
#   * routines that need contiguous memory accept strided types too, but a
#     non-contiguous argument is rejected with an `ArgumentError` instead of being
#     read/written as if it were dense.
#
# Every check here cross-validates a view against the equivalent dense `Vector`/`Matrix`.

@testset "Strided array contract" begin

AA = AppleAccelerate

# Three ways to present the same logical vector `x` without it being a `Vector`.
_contig(x)  = view([x[1:2]; x; x[1:2]], 3:length(x)+2)               # unit stride, offset
_stepped(x) = (b = repeat(x, inner=2); b[2:2:end] .= 0; view(b, 1:2:length(b)))  # stride 2
_reversed(x) = view(reverse(x), length(x):-1:1)                        # stride -1
_VIEWS = (("contiguous view", _contig), ("stride-2 view", _stepped),
                ("reversed view", _reversed))

# A stride-2 output buffer pre-filled with a sentinel, so we can prove that a strided
# `result` is written only at the selected positions.
function _sentinel_out(::Type{S}, n) where S
    buf = fill(S(-7), 2n)
    return buf, view(buf, 1:2:2n)
end

@testset "complexarray.jl stride passthrough::$T" for T in (Float32, Float64)
    CT = Complex{T}
    n = 33
    X, Y, Z = randn(CT, n), randn(CT, n), randn(CT, n)
    R = rand(T, n) .+ T(0.5)
    c = CT(0.3, -1.2)

    cases = (
        (AA.vneg,   (X,)),      (AA.vconj,  (X,)),      (AA.vcopy,  (X,)),
        (AA.vmul,   (X, Y)),    (AA.vdiv,   (X, Y)),    (AA.vabs,   (X,)),
        (AA.vphase, (X,)),      (AA.vmags,  (X,)),      (AA.vmagsa, (X, R)),
        (AA.dot,    (X, Y)),    (AA.dotu,   (X, Y)),
        (AA.zvadd,  (X, Y)),    (AA.zvsub,  (X, Y)),    (AA.zvcmul, (X, Y)),
        (AA.zrvmul, (X, R)),    (AA.zrvdiv, (X, R)),    (AA.zrvadd, (X, R)),
        (AA.zrvsub, (X, R)),    (AA.zvcma,  (X, Y, Z)), (AA.zvma,   (X, Y, Z)),
        (AA.zidotpr, (X, Y)),   (AA.zrdotpr, (X, R)),
        (AA.zvmmaa, (X, Y, Z, X, Y)),
    )
    @testset "$name" for (name, mk) in _VIEWS
        for (f, args) in cases
            @test f(map(mk, args)...) ≈ f(args...)
        end
        @test AA.vsmul(mk(X), c) ≈ AA.vsmul(X, c)
        @test AA.zvsma(mk(X), c, mk(Y)) ≈ AA.zvsma(X, c, Y)
        @test AA.zconv(mk(X), mk(Y[1:5])) ≈ AA.zconv(X, Y[1:5])
    end

    @testset "strided outputs only touch selected elements" begin
        buf, out = _sentinel_out(CT, n)
        AA.vmul!(out, _stepped(X), _reversed(Y))
        @test out ≈ X .* Y
        @test all(==(CT(-7)), buf[2:2:end])

        buf, out = _sentinel_out(CT, n)
        AA.zvadd!(out, X, _contig(Y))
        @test out ≈ X .+ Y
        @test all(==(CT(-7)), buf[2:2:end])

        rbuf, rout = _sentinel_out(T, n)
        AA.vabs!(rout, _stepped(X))
        @test rout ≈ abs.(X)
        @test all(==(T(-7)), rbuf[2:2:end])

        buf, out = _sentinel_out(CT, n)
        AA.zvfill!(out, c)
        @test all(==(c), out)
        @test all(==(CT(-7)), buf[2:2:end])
    end

    @testset "contiguous-only routines" begin
        # accepted as (contiguous) views …
        @test all(AA.polar(_contig(X)) .≈ AA.polar(X))
        @test all(AA.ctoz(_contig(X)) .≈ AA.ctoz(X))
        re, im = AA.ctoz(X)
        @test AA.ztoc(_contig(re), _contig(im)) ≈ X
        # … and rejected, not misread, when non-contiguous.
        @test_throws ArgumentError AA.polar(_stepped(X))
        @test_throws ArgumentError AA.ctoz(_stepped(X))
        @test_throws ArgumentError AA.ztoc(_stepped(re), _stepped(im))
    end

    @testset "complex matrix ops" begin
        A, B = randn(CT, 4, 6), randn(CT, 6, 5)
        big = randn(CT, 8, 6)
        @test AA.zmmul(view(A, :, :), view(B, :, 1:5)) ≈ A * B
        Av = view(big, 1:2:8, :)                       # non-contiguous rows
        @test_throws ArgumentError AA.zmmul(Av, B)
        @test_throws ArgumentError AA.zmmul!(view(randn(CT, 8, 5), 1:2:8, :), A, B)
    end
end

@testset "dsp.jl filters::$T" for T in (Float32, Float64)
    n = 40
    X, K = randn(T, n), randn(T, 5)
    @testset "$name" for (name, mk) in _VIEWS
        @test AA.conv(mk(X), mk(K)) ≈ AA.conv(X, K)
        @test AA.xcorr(mk(X), mk(K)) ≈ AA.xcorr(X, K)

        coeffs = T[0.2, 0.1, 0.05, -0.3, 0.1]
        bq = AA.biquadcreate(Float64.(coeffs), 1, T)
        @test AA.biquad(mk(X), zeros(T, 4), n, bq) ≈ AA.biquad(X, zeros(T, 4), n, bq)

        @test AA.deq22(mk(X), coeffs) ≈ AA.deq22(X, coeffs)
    end

    @testset "strided conv! output only touches selected elements" begin
        buf, out = _sentinel_out(T, n + 4)
        AA.conv!(out, X, K)
        @test out ≈ AA.conv(X, K)
        @test all(==(T(-7)), buf[2:2:end])
    end

    @testset "contiguous-only arguments" begin
        bq = AA.biquadcreate(Float64[0.2, 0.1, 0.05, -0.3, 0.1], 1, T)
        @test_throws ArgumentError AA.biquad(X, _stepped(zeros(T, 4)), n, bq)   # delays
        @test_throws ArgumentError AA.deq22(X, _stepped(T[0.2, 0.1, 0.05, -0.3, 0.1]))
        @test_throws ArgumentError AA.biquadcreate(_stepped(Float64[0.2, 0.1, 0.05, -0.3, 0.1]), 1, T)

        F = randn(T, 4)
        @test AA.desamp(_contig(X), 2, _contig(F)) ≈ AA.desamp(X, 2, F)
        @test_throws ArgumentError AA.desamp(_stepped(X), 2, F)

        w = zeros(T, 16)
        @test AA.hanning!(view(w, 1:16), 16) ≈ AA.hanning(16, T)
        @test_throws ArgumentError AA.hanning!(view(zeros(T, 32), 1:2:32), 16)
        @test_throws ArgumentError AA.hamming!(view(zeros(T, 32), 1:2:32), 16)
        @test_throws ArgumentError AA.blackman!(view(zeros(T, 32), 1:2:32), 16)
    end
end

@testset "dsp.jl transforms::$T" for T in (Float32, Float64)
    CT = Complex{T}
    n = 64
    x, r = randn(CT, n), randn(T, n)
    @testset "$name" for (name, mk) in _VIEWS
        @test AA.fft(mk(x)) ≈ AA.fft(x)
        @test AA.bfft(mk(x)) ≈ AA.bfft(x)
        @test AA.ifft(mk(x)) ≈ AA.ifft(x)
        @test AA.rfft(mk(r)) ≈ AA.rfft(r)
        @test AA.irfft(mk(AA.rfft(r)), n) ≈ AA.irfft(AA.rfft(r), n)
        @test AA.fft(mk(x[1:24])) ≈ AA.fft(x[1:24])    # mixed-radix path
        @test AA.dft(mk(x)) ≈ AA.dft(x)
        @test AA.dft_interleaved(mk(x)) ≈ AA.dft_interleaved(x)
        @test AA.fftradix3(mk(x[1:48])) ≈ AA.fftradix3(x[1:48])

        setup, ws = AA.plan_fft(n, T), AA.FFTWorkspace{T}(n)
        @test AA.fft(mk(x), setup, ws) ≈ AA.fft(x)
        @test AA.rfft(mk(r), setup, ws) ≈ AA.rfft(r)

        y = copy(x); v = mk(y)
        @test AA.fft!(v) === v
        @test v ≈ AA.fft(x)
    end

    @testset "in-place real FFT needs contiguous memory" begin
        y = copy(r)
        @test AA.rfft!(_contig(y)) ≈ AA.rfft!(copy(r))
        @test_throws ArgumentError AA.rfft!(_stepped(y))
    end

    @testset "2-D and batched" begin
        M = randn(CT, 16, 12)
        blk = view(M, :, 1:8)                          # contiguous block
        rows = view(randn(CT, 32, 8) .= 0, 1:2:32, :)  # non-contiguous
        rows .= blk
        @test AA.fft(blk) ≈ AA.fft(Matrix(blk))
        @test AA.fft(rows) ≈ AA.fft(Matrix(blk))
        @test AA.fft(rows, 1) ≈ AA.fft(Matrix(blk), 1)
        Rm = real.(Matrix(blk))
        @test AA.rfft(view(Rm, :, :)) ≈ AA.rfft(Rm)
        @test AA.rfft(view(Rm, :, :), 1) ≈ AA.rfft(Rm, 1)
    end

    if T === Float32
        d = randn(Float32, 64)
        @test AA.dct(_contig(d)) ≈ AA.dct(d)
        @test_throws ArgumentError AA.dct(_stepped(d))
        z = randn(ComplexF32, 16)
        @test AA.fft16(_stepped(z)) ≈ AA.fft16(z)
    end
end

@testset "array.jl matrix ops::$T" for T in (Float32, Float64)
    A, B = randn(T, 4, 6), randn(T, 6, 5)
    @test AA.mmul(view(A, :, :), view(randn(T, 6, 7) .= 0, :, 1:5) .= B) ≈ A * B
    @test AA.mtrans(view(A, :, 1:3)) ≈ permutedims(A[:, 1:3])
    @test AA.mmov(view(A, :, 2:4)) == A[:, 2:4]
    F3 = randn(T, 3, 3)
    img = randn(T, 8, 8)
    @test AA.f3x3(view(img, :, :), F3) ≈ AA.f3x3(img, F3)

    Anc = view(randn(T, 8, 6), 1:2:8, :)               # non-contiguous rows
    @test_throws ArgumentError AA.mmul(Anc, B)
    @test_throws ArgumentError AA.mtrans(Anc)
    @test_throws ArgumentError AA.mmov(Anc)
    @test_throws ArgumentError AA.f3x3(view(randn(T, 16, 8), 1:2:16, :), F3)
    @test_throws DimensionMismatch AA.mmov!(zeros(T, 2, 2), A)
end

end # @testset "Strided array contract"
