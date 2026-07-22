@testset "Complex Array Operations (Batch 5)" begin

for T in (Float32, Float64)
    @testset "Complex-Complex Binary::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # zvadd
        @test AppleAccelerate.zvadd(A, B) ≈ A .+ B
        R = similar(A)
        AppleAccelerate.zvadd!(R, A, B)
        @test R ≈ A .+ B

        # zvsub
        @test AppleAccelerate.zvsub(A, B) ≈ A .- B
        AppleAccelerate.zvsub!(R, A, B)
        @test R ≈ A .- B

        # zvcmul: conj(A) * B
        @test AppleAccelerate.zvcmul(A, B) ≈ conj.(A) .* B
        AppleAccelerate.zvcmul!(R, A, B)
        @test R ≈ conj.(A) .* B
    end
end

for T in (Float32, Float64)
    @testset "Complex-Real Binary::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{T} = randn(T, N)

        # zrvmul: complex * real
        @test AppleAccelerate.zrvmul(A, B) ≈ A .* B
        R = similar(A)
        AppleAccelerate.zrvmul!(R, A, B)
        @test R ≈ A .* B

        # zrvadd: complex + real (adds to real part)
        expected_add = complex.(real.(A) .+ B, imag.(A))
        @test AppleAccelerate.zrvadd(A, B) ≈ expected_add
        AppleAccelerate.zrvadd!(R, A, B)
        @test R ≈ expected_add

        # zrvsub: complex - real (subtracts from real part)
        expected_sub = complex.(real.(A) .- B, imag.(A))
        @test AppleAccelerate.zrvsub(A, B) ≈ expected_sub
        AppleAccelerate.zrvsub!(R, A, B)
        @test R ≈ expected_sub

        # zrvdiv: complex / real (elementwise)
        Bnz::Vector{T} = B .+ T(2)  # ensure nonzero divisor
        expected_div = A ./ Bnz
        @test AppleAccelerate.zrvdiv(A, Bnz) ≈ expected_div
        AppleAccelerate.zrvdiv!(R, A, Bnz)
        @test R ≈ expected_div
    end
end

for T in (Float32, Float64)
    @testset "Complex Compound::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        C::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # zvcma: conj(A)*B + C
        @test AppleAccelerate.zvcma(A, B, C) ≈ conj.(A) .* B .+ C
        R = similar(A)
        AppleAccelerate.zvcma!(R, A, B, C)
        @test R ≈ conj.(A) .* B .+ C

        # zvma: A*B + C
        @test AppleAccelerate.zvma(A, B, C) ≈ A .* B .+ C
        AppleAccelerate.zvma!(R, A, B, C)
        @test R ≈ A .* B .+ C

        # zvsma: A*b + C (b is complex scalar)
        b::Complex{T} = complex(T(2.5), T(-1.3))
        @test AppleAccelerate.zvsma(A, b, C) ≈ A .* b .+ C
        AppleAccelerate.zvsma!(R, A, b, C)
        @test R ≈ A .* b .+ C
    end
end

for T in (Float32, Float64)
    @testset "Complex Dot Products::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # zidotpr: conjugate dot — sum(conj(A) .* B)
        @test AppleAccelerate.zidotpr(A, B) ≈ sum(conj.(A) .* B)

        # zrdotpr: complex-real dot
        Br::Vector{T} = randn(T, N)
        @test AppleAccelerate.zrdotpr(A, Br) ≈ sum(A .* Br)

        # Regression: complex `dot` must be conjugated to match LinearAlgebra.dot
        # (sum(conj(X) .* Y)), and `dotu` exposes the un-conjugated product.
        @test AppleAccelerate.dot(A, B) ≈ LinearAlgebra.dot(A, B) ≈ sum(conj.(A) .* B)
        @test AppleAccelerate.dotu(A, B) ≈ sum(A .* B)
        # The two differ for genuinely complex inputs.
        @test !(AppleAccelerate.dot(A, B) ≈ AppleAccelerate.dotu(A, B))
    end
end

for T in (Float32, Float64)
    @testset "Complex Fill::$T" begin
        c = complex(T(3.0), T(-1.5))
        R = Vector{Complex{T}}(undef, 100)
        AppleAccelerate.zvfill!(R, c)
        @test all(==(c), R)
    end
end

for T in (Float32, Float64)
    @testset "Complex Convolution::$T" begin
        X = complex.(randn(T, 50), randn(T, 50))
        K = complex.(randn(T, 5), randn(T, 5))
        result = AppleAccelerate.zconv(X, K)
        @test length(result) == length(X) + length(K) - 1
        @test eltype(result) == Complex{T}
        # Numerical correctness: vDSP_zconv applies the kernel without reversal
        # (correlation-style), matching DSP.conv with a reversed kernel.
        @test result ≈ DSP.conv(X, reverse(K)) atol=T(1e-3)

        # Regression: an oversized `result` buffer must not make vDSP read past
        # the zero-padded input (padding is sized from length(result)).
        extra = 7
        big = fill(Complex{T}(NaN), length(result) + extra)
        AppleAccelerate.zconv!(big, X, K)
        @test big[1:length(result)] ≈ result atol=T(1e-3)
        @test !any(isnan, big)

        # Too-short result must still be rejected.
        @test_throws ErrorException AppleAccelerate.zconv!(
            similar(big, length(X) + length(K) - 2), X, K)
    end
end

for T in (Float32, Float64)
    @testset "Complex DimensionMismatch::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Rbad = similar(A, N - 1)
        @test_throws DimensionMismatch AppleAccelerate.zvadd!(Rbad, A, B)
        @test_throws DimensionMismatch AppleAccelerate.zvsub!(Rbad, A, B)
    end
end

for T in (Float32, Float64)
    @testset "Complex Matrix Multiply::$T" begin
        m, p, n = 4, 3, 5
        A = complex.(randn(T, m, p), randn(T, m, p))
        B = complex.(randn(T, p, n), randn(T, p, n))

        C = AppleAccelerate.zmmul(A, B)
        @test C ≈ A * B atol=T(1e-3)

        C2 = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmmul!(C2, A, B)
        @test C2 ≈ A * B atol=T(1e-3)
    end

    @testset "Complex Matrix Multiply-Add/Subtract::$T" begin
        m, p, n = 4, 3, 5
        A = complex.(randn(T, m, p), randn(T, m, p))
        B = complex.(randn(T, p, n), randn(T, p, n))
        C = complex.(randn(T, m, n), randn(T, m, n))

        # zmma: D = A*B + C
        D = AppleAccelerate.zmma(A, B, C)
        @test D ≈ A * B + C atol=T(1e-3)
        Dm = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmma!(Dm, A, B, C)
        @test Dm ≈ A * B + C atol=T(1e-3)

        # zmms: D = A*B - C
        E = AppleAccelerate.zmms(A, B, C)
        @test E ≈ A * B - C atol=T(1e-3)
        Em = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmms!(Em, A, B, C)
        @test Em ≈ A * B - C atol=T(1e-3)

        # dimension mismatches
        @test_throws DimensionMismatch AppleAccelerate.zmma(A, complex.(randn(T, p+1, n), randn(T, p+1, n)), C)
        @test_throws DimensionMismatch AppleAccelerate.zmma!(Matrix{Complex{T}}(undef, m, n), A, B, complex.(randn(T, m, n+1), randn(T, m, n+1)))
        @test_throws DimensionMismatch AppleAccelerate.zmms!(Matrix{Complex{T}}(undef, m+1, n), A, B, C)
    end
end

# ============================================================
# Unary complex ops: vneg, vconj, vcopy
# ============================================================
for T in (Float32, Float64)
    @testset "Complex Unary::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # vneg: -X
        @test AppleAccelerate.vneg(X) ≈ -X
        Rc = similar(X)
        @test AppleAccelerate.vneg!(Rc, X) === Rc
        @test Rc ≈ .-X

        # vconj: conj(X)
        @test AppleAccelerate.vconj(X) ≈ conj.(X)
        @test AppleAccelerate.vconj!(Rc, X) === Rc
        @test Rc ≈ conj.(X)

        # vcopy: exact copy (bit-identical, not just approx)
        Xc = AppleAccelerate.vcopy(X)
        @test Xc == X
        @test Xc !== X

        # DimensionMismatch paths
        bad = similar(X, N - 1)
        @test_throws DimensionMismatch AppleAccelerate.vneg!(bad, X)
        @test_throws DimensionMismatch AppleAccelerate.vconj!(bad, X)
    end
end

# ============================================================
# Binary complex ops: vmul, vdiv, vsmul
# ============================================================
for T in (Float32, Float64)
    @testset "Complex Binary::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Y::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        # ensure nonzero divisor magnitudes
        Ynz::Vector{Complex{T}} = Y .+ Complex{T}(2, 2)

        # vmul: elementwise X .* Y
        @test AppleAccelerate.vmul(X, Y) ≈ X .* Y
        Rc = similar(X)
        @test AppleAccelerate.vmul!(Rc, X, Y) === Rc
        @test Rc ≈ X .* Y

        # vdiv: elementwise X ./ Y
        @test AppleAccelerate.vdiv(X, Ynz) ≈ X ./ Ynz
        @test AppleAccelerate.vdiv!(Rc, X, Ynz) === Rc
        @test Rc ≈ X ./ Ynz

        # vsmul: X .* scalar
        c = complex(T(1.7), T(-0.9))
        @test AppleAccelerate.vsmul(X, c) ≈ X .* c
        @test AppleAccelerate.vsmul!(Rc, X, c) === Rc
        @test Rc ≈ X .* c

        # DimensionMismatch paths
        bad = similar(X, N - 1)
        @test_throws DimensionMismatch AppleAccelerate.vmul!(bad, X, Y)
        @test_throws DimensionMismatch AppleAccelerate.vdiv!(bad, X, Y)
        @test_throws DimensionMismatch AppleAccelerate.vsmul!(bad, X, c)
    end
end

# ============================================================
# Complex → Real ops: vabs, vphase, vmags, vmagsa
# ============================================================
for T in (Float32, Float64)
    @testset "Complex To Real::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # vabs: abs (modulus)
        @test AppleAccelerate.vabs(X) ≈ abs.(X)
        Rr = Vector{T}(undef, N)
        @test AppleAccelerate.vabs!(Rr, X) === Rr
        @test Rr ≈ abs.(X)

        # vphase: angle
        @test AppleAccelerate.vphase(X) ≈ angle.(X)
        @test AppleAccelerate.vphase!(Rr, X) === Rr
        @test Rr ≈ angle.(X)

        # vmags: abs2
        @test AppleAccelerate.vmags(X) ≈ abs2.(X)
        @test AppleAccelerate.vmags!(Rr, X) === Rr
        @test Rr ≈ abs2.(X)

        # vmagsa: abs2(X) + B
        B::Vector{T} = randn(T, N)
        @test AppleAccelerate.vmagsa(X, B) ≈ abs2.(X) .+ B
        @test AppleAccelerate.vmagsa!(Rr, X, B) === Rr
        @test Rr ≈ abs2.(X) .+ B

        # DimensionMismatch paths
        badr = Vector{T}(undef, N - 1)
        @test_throws DimensionMismatch AppleAccelerate.vabs!(badr, X)
        @test_throws DimensionMismatch AppleAccelerate.vphase!(badr, X)
        @test_throws DimensionMismatch AppleAccelerate.vmags!(badr, X)
        @test_throws DimensionMismatch AppleAccelerate.vmagsa!(badr, X, randn(T, N))
        @test_throws DimensionMismatch AppleAccelerate.vmagsa!(Rr, X, randn(T, N - 1))
    end
end

# ============================================================
# dot / dotu DimensionMismatch and degenerate-equality cases
# ============================================================
for T in (Float32, Float64)
    @testset "Complex dot edge cases::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        @test_throws DimensionMismatch AppleAccelerate.dot(A, B[1:end-1])
        @test_throws DimensionMismatch AppleAccelerate.dotu(A, B[1:end-1])
        @test_throws DimensionMismatch AppleAccelerate.zidotpr(A, B[1:end-1])
        @test_throws DimensionMismatch AppleAccelerate.zrdotpr(A, randn(T, N - 1))

        # For real-only vectors, dot and dotu coincide (conjugation is a no-op).
        Ar::Vector{Complex{T}} = complex.(randn(T, N))
        Br::Vector{Complex{T}} = complex.(randn(T, N))
        @test AppleAccelerate.dot(Ar, Br) ≈ AppleAccelerate.dotu(Ar, Br)
    end
end

# ============================================================
# polar / rect coordinate conversion and round-trip
# ============================================================
for T in (Float32, Float64)
    @testset "Polar/Rect::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        mags, angs = AppleAccelerate.polar(X)
        @test mags ≈ abs.(X)
        @test angs ≈ angle.(X)

        # rect: (mag, angle) → complex, vs mag .* cis(angle)
        R = AppleAccelerate.rect(mags, angs)
        @test R ≈ mags .* cis.(angs)

        # round-trip polar -> rect recovers X
        @test AppleAccelerate.rect(mags, angs) ≈ X

        @test_throws DimensionMismatch AppleAccelerate.rect(mags, angs[1:end-1])
    end
end

# ============================================================
# ctoz / ztoc split-complex format conversion + round-trip
# (exercises the GC-safe split-view constructors)
# ============================================================
for T in (Float32, Float64)
    @testset "ctoz/ztoc::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        re, im = AppleAccelerate.ctoz(X)
        @test re == real.(X)
        @test im == imag.(X)

        # ztoc reassembles interleaved complex from split-complex
        Xr = AppleAccelerate.ztoc(re, im)
        @test Xr == X

        # full round-trip
        re2, im2 = AppleAccelerate.ctoz(AppleAccelerate.ztoc(re, im))
        @test re2 == re
        @test im2 == im

        # mismatched re/im lengths must throw (count is taken from re, so a
        # shorter im would be read out of bounds)
        @test_throws DimensionMismatch AppleAccelerate.ztoc(re, im[1:end-1])
    end
end

# ============================================================
# zvfill!/zmmul DimensionMismatch, zvcmul/zrvmul/etc mismatch paths
# ============================================================
for T in (Float32, Float64)
    @testset "Batch5 mismatch paths::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        C::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Br::Vector{T} = randn(T, N)
        bad = similar(A, N - 1)
        badr = similar(Br, N - 1)

        @test_throws DimensionMismatch AppleAccelerate.zvcmul!(bad, A, B)
        @test_throws DimensionMismatch AppleAccelerate.zrvmul!(bad, A, Br)
        @test_throws DimensionMismatch AppleAccelerate.zrvdiv!(bad, A, Br)
        @test_throws DimensionMismatch AppleAccelerate.zrvadd!(bad, A, Br)
        @test_throws DimensionMismatch AppleAccelerate.zrvsub!(bad, A, Br)
        @test_throws DimensionMismatch AppleAccelerate.zrvmul!(similar(A), A, badr)
        @test_throws DimensionMismatch AppleAccelerate.zvcma!(bad, A, B, C)
        @test_throws DimensionMismatch AppleAccelerate.zvma!(bad, A, B, C)
        @test_throws DimensionMismatch AppleAccelerate.zvsma!(bad, A, complex(T(1), T(1)), C)

        # zmmul dimension mismatches
        M = complex.(randn(T, 3, 4), randn(T, 3, 4))
        Nn = complex.(randn(T, 5, 2), randn(T, 5, 2))  # 4 ≠ 5
        @test_throws DimensionMismatch AppleAccelerate.zmmul(M, Nn)
        Mok = complex.(randn(T, 3, 4), randn(T, 3, 4))
        Bok = complex.(randn(T, 4, 2), randn(T, 4, 2))
        Cbad = Matrix{Complex{T}}(undef, 3, 3)  # wrong shape
        @test_throws DimensionMismatch AppleAccelerate.zmmul!(Cbad, Mok, Bok)
    end
end

# ============================================================
# Complex matrix multiply-add / multiply-subtract / reverse-subtract
# (zmma, zmms, zmsm), vector multiply-multiply-add-add (zvmmaa), and
# complex-real decimating resample (zrdesamp)
# ============================================================
for T in (Float32, Float64)
    @testset "Complex Matrix Multiply-Add/Subtract::$T" begin
        # non-square, non-symmetric so a layout/order bug fails loudly
        m, p, n = 4, 3, 5
        A = complex.(randn(T, m, p), randn(T, m, p))
        B = complex.(randn(T, p, n), randn(T, p, n))
        C = complex.(randn(T, m, n), randn(T, m, n))
        rtol = T == Float32 ? sqrt(eps(Float32)) : sqrt(eps(Float64))

        # zmma: D = A*B + C
        D = AppleAccelerate.zmma(A, B, C)
        @test D ≈ A * B .+ C rtol=rtol
        D2 = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmma!(D2, A, B, C)
        @test D2 ≈ A * B .+ C rtol=rtol

        # zmms: D = A*B - C
        D3 = AppleAccelerate.zmms(A, B, C)
        @test D3 ≈ A * B .- C rtol=rtol
        D4 = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmms!(D4, A, B, C)
        @test D4 ≈ A * B .- C rtol=rtol

        # zmsm: D = C - A*B (vDSP's reverse subtract, confirmed empirically —
        # NOT (A-B)*C)
        D5 = AppleAccelerate.zmsm(A, B, C)
        @test D5 ≈ C .- A * B rtol=rtol
        D6 = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmsm!(D6, A, B, C)
        @test D6 ≈ C .- A * B rtol=rtol

        # DimensionMismatch guards
        Bbad = complex.(randn(T, p + 1, n), randn(T, p + 1, n))
        @test_throws DimensionMismatch AppleAccelerate.zmma(A, Bbad, C)
        @test_throws DimensionMismatch AppleAccelerate.zmms(A, Bbad, C)
        @test_throws DimensionMismatch AppleAccelerate.zmsm(A, Bbad, C)
        Cbad = complex.(randn(T, m, n + 1), randn(T, m, n + 1))
        @test_throws DimensionMismatch AppleAccelerate.zmma(A, B, Cbad)
        @test_throws DimensionMismatch AppleAccelerate.zmms(A, B, Cbad)
        @test_throws DimensionMismatch AppleAccelerate.zmsm(A, B, Cbad)
        Dbad = Matrix{Complex{T}}(undef, m, n + 1)
        @test_throws DimensionMismatch AppleAccelerate.zmma!(Dbad, A, B, C)
        @test_throws DimensionMismatch AppleAccelerate.zmms!(Dbad, A, B, C)
        @test_throws DimensionMismatch AppleAccelerate.zmsm!(Dbad, A, B, C)
    end

    @testset "Complex Vector Multiply-Multiply-Add-Add (zvmmaa)::$T" begin
        Av::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Bv::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Cv::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Dv::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Ev::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        rtol = T == Float32 ? sqrt(eps(Float32)) : sqrt(eps(Float64))

        ref = Av .* Bv .+ Cv .* Dv .+ Ev
        F = AppleAccelerate.zvmmaa(Av, Bv, Cv, Dv, Ev)
        @test F ≈ ref rtol=rtol

        F2 = similar(Av)
        @test AppleAccelerate.zvmmaa!(F2, Av, Bv, Cv, Dv, Ev) === F2
        @test F2 ≈ ref rtol=rtol

        bad = similar(Av, N - 1)
        @test_throws DimensionMismatch AppleAccelerate.zvmmaa!(bad, Av, Bv, Cv, Dv, Ev)
        @test_throws DimensionMismatch AppleAccelerate.zvmmaa!(F2, bad, Bv, Cv, Dv, Ev)
    end

    @testset "Complex-Real Decimating Resample (zrdesamp)::$T" begin
        Av::Vector{Complex{T}} = complex.(randn(T, 41), randn(T, 41))
        Fv::Vector{T} = randn(T, 7)
        DF = 3
        rtol = T == Float32 ? sqrt(eps(Float32)) : sqrt(eps(Float64))

        P = length(Fv)
        Nout = div(length(Av) - P, DF) + 1

        Cv = AppleAccelerate.zrdesamp(Av, DF, Fv)
        @test length(Cv) == Nout

        # manual decimated-correlation reference (0-indexed n, p)
        Cref = Vector{Complex{T}}(undef, Nout)
        for n in 0:(Nout - 1)
            s = Complex{T}(0)
            for p in 0:(P - 1)
                s += Av[n * DF + p + 1] * Fv[p + 1]
            end
            Cref[n + 1] = s
        end
        @test Cv ≈ Cref rtol=rtol

        Cout = Vector{Complex{T}}(undef, Nout)
        AppleAccelerate.zrdesamp!(Cout, Av, DF, Fv)
        @test Cout ≈ Cref rtol=rtol

        # DF=1 (pure FIR correlation, no decimation)
        C1 = AppleAccelerate.zrdesamp(Av, 1, Fv)
        Nout1 = div(length(Av) - P, 1) + 1
        Cref1 = Vector{Complex{T}}(undef, Nout1)
        for n in 0:(Nout1 - 1)
            s = Complex{T}(0)
            for p in 0:(P - 1)
                s += Av[n * 1 + p + 1] * Fv[p + 1]
            end
            Cref1[n + 1] = s
        end
        @test C1 ≈ Cref1 rtol=rtol

        # errors when filter longer than signal
        Ashort = complex.(randn(T, 3), randn(T, 3))
        Flong = randn(T, 5)
        @test_throws ErrorException AppleAccelerate.zrdesamp(Ashort, 1, Flong)
        Coutbad = Vector{Complex{T}}(undef, 10)
        @test_throws ErrorException AppleAccelerate.zrdesamp!(Coutbad, Ashort, 1, Flong)
    end
end

end # @testset
