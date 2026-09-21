# Plan interface (src/fftplan.jl): `plan * x`, `mul!`, `inv`, `\`, introspection.
# Every path is cross-validated against FFTW for Float32 and Float64. FFTW is always
# called through its one-shot API on the original input, never through a MEASURE plan
# (which would overwrite the buffer it is handed).

const AA = AppleAccelerate

# Run `mul!` once to compile, then measure. A function barrier keeps the globals out.
function _plan_allocs(y, p, x)
    mul!(y, p, x)
    return @allocated mul!(y, p, x)
end

_rtol(::Type{T}) where {T} = T === Float32 ? 2f-3 : 1e-10

@testset "FFTPlan" begin

@testset "1-D complex ($T)" for T in (Float32, Float64)
    # 1, 2, 4: below the interleaved DFT's minimum → stride-2 vDSP_fft_zop.
    # 8 … 4096 and 24/40/120/960: vDSP_DFT_Interleaved. 8192: vDSP_fft_zop.
    # 3072, 15360: mixed radix beyond the interleaved limit → vDSP_DFT_Execute.
    for n in (1, 2, 4, 8, 64, 1024, 4096, 8192, 24, 40, 120, 960, 3072, 15360)
        x = randn(Complex{T}, n)
        xc = copy(x)
        pf = AA.fftplan(x); pb = AA.bfftplan(x); pi = AA.ifftplan(x)
        X = FFTW.fft(x)
        @test pf * x ≈ X rtol=_rtol(T)
        @test x == xc                                  # input untouched
        @test pb * x ≈ FFTW.bfft(x) rtol=_rtol(T)
        @test pi * x ≈ FFTW.ifft(x) rtol=_rtol(T)
        @test pf \ (pf * x) ≈ x rtol=_rtol(T)
        @test inv(pb) * (pb * x) ≈ x rtol=_rtol(T)
        @test inv(pi) * x ≈ X rtol=_rtol(T)            # inv(ifft) is the forward FFT
        @test inv(inv(pf)) * x ≈ X rtol=_rtol(T)

        y = similar(x)
        @test mul!(y, pf, x) === y
        @test y ≈ X rtol=_rtol(T)
        z = copy(x)
        @test mul!(z, pf, z) === z                      # in place
        @test z ≈ X rtol=_rtol(T)
        z = copy(x); mul!(z, pi, z)
        @test z ≈ FFTW.ifft(x) rtol=_rtol(T)
    end

    # Shape constructors build the same plan as the array ones.
    x = randn(Complex{T}, 256)
    @test AA.fftplan(Complex{T}, 256) * x == AA.fftplan(x) * x
    @test AA.ifftplan(Complex{T}, (256,)) * x == AA.ifftplan(x) * x
    @test AA.bfftplan(Complex{T}, 256, 1) * x == AA.bfftplan(x) * x
end

@testset "2-D complex ($T)" for T in (Float32, Float64)
    for sz in ((8, 32), (64, 4), (16, 16), (1, 8), (2, 2))   # non-square catches layout bugs
        x = randn(Complex{T}, sz)
        pf = AA.fftplan(x); pb = AA.bfftplan(x); pi = AA.ifftplan(x)
        X = FFTW.fft(x)
        @test pf * x ≈ X rtol=_rtol(T)
        @test pb * x ≈ FFTW.bfft(x) rtol=_rtol(T)
        @test pi * x ≈ FFTW.ifft(x) rtol=_rtol(T)
        @test pf \ (pf * x) ≈ x rtol=_rtol(T)
        @test inv(pb) * (pb * x) ≈ x rtol=_rtol(T)
        z = copy(x); mul!(z, pf, z)
        @test z ≈ X rtol=_rtol(T)
    end
    x = randn(Complex{T}, 8, 32)
    @test AA.fftplan(Complex{T}, (8, 32)) * x == AA.fftplan(x) * x
end

@testset "batched 1-D complex along dims ($T)" for T in (Float32, Float64)
    for sz in ((8, 32), (64, 4), (16, 5), (3, 16)), dims in (1, 2)
        ispow2(sz[dims]) || continue                   # only the transformed dim must be 2^k
        x = randn(Complex{T}, sz)
        pf = AA.fftplan(x, dims); pb = AA.bfftplan(x, dims); pi = AA.ifftplan(x, dims)
        X = FFTW.fft(x, dims)
        @test pf * x ≈ X rtol=_rtol(T)
        @test pb * x ≈ FFTW.bfft(x, dims) rtol=_rtol(T)
        @test pi * x ≈ FFTW.ifft(x, dims) rtol=_rtol(T)
        @test pf \ (pf * x) ≈ x rtol=_rtol(T)
        @test inv(pb) * (pb * x) ≈ x rtol=_rtol(T)
        @test pf * x ≈ AA.fft(x, dims) rtol=_rtol(T)   # agrees with the one-shot API
        z = copy(x); mul!(z, pf, z)
        @test z ≈ X rtol=_rtol(T)
    end
end

@testset "1-D real ($T)" for T in (Float32, Float64)
    # 48, 80, 96, 160 are mixed-radix lengths vDSP's real DFT accepts in both precisions.
    for n in (2, 4, 8, 64, 1024, 8192, 48, 80, 96, 160)
        x = randn(T, n)
        xc = copy(x)
        pf = AA.rfftplan(x)
        X = FFTW.rfft(x)
        @test AA.output_size(pf) == size(X)
        Y = pf * x
        @test Y ≈ X rtol=_rtol(T)
        @test x == xc
        Xc = copy(X)
        pb = AA.brfftplan(X, n); pi = AA.irfftplan(X, n)
        @test pb * X ≈ FFTW.brfft(X, n) rtol=_rtol(T)
        @test X == Xc                                  # the spectrum is not clobbered
        @test pi * X ≈ x rtol=_rtol(T)
        @test pf \ Y ≈ x rtol=_rtol(T)
        @test inv(pb) * (pb * X) ≈ X rtol=_rtol(T)
        @test inv(pi) * x ≈ X rtol=_rtol(T)
    end
    @test AA.rfftplan(T, 64) * ones(T, 64) ≈ FFTW.rfft(ones(T, 64))
end

@testset "2-D real ($T)" for T in (Float32, Float64)
    for sz in ((8, 16), (32, 4), (2, 2), (16, 16))
        x = randn(T, sz)
        pf = AA.rfftplan(x)
        X = FFTW.rfft(x)
        @test AA.output_size(pf) == size(X)
        @test pf * x ≈ X rtol=_rtol(T)
        @test AA.brfftplan(X, sz[1]) * X ≈ FFTW.brfft(X, sz[1]) rtol=_rtol(T)
        @test AA.irfftplan(X, sz[1]) * X ≈ x rtol=_rtol(T)
        @test pf \ (pf * x) ≈ x rtol=_rtol(T)
    end
end

@testset "DCT" begin
    # vDSP's DCTs are FFTW's REDFT10 / REDFT01 / REDFT11 without the factor of 2.
    kinds = Dict(2 => FFTW.REDFT10, 3 => FFTW.REDFT01, 4 => FFTW.REDFT11)
    for n in (16, 48, 256, 240), t in (2, 3, 4)
        x = randn(Float32, n)
        p = AA.dctplan(x, t)
        @test p * x ≈ AA.dct(x, t)
        @test p * x ≈ Float32.(FFTW.r2r(Float64.(x), kinds[t]) ./ 2) rtol=1f-3
        @test p \ (p * x) ≈ x rtol=1f-3
        @test inv(inv(p)) * x ≈ p * x rtol=1f-3
    end
    @test AA.dctplan(Float32, 64) * ones(Float32, 64) ≈ AA.dct(ones(Float32, 64))
    @test_throws ArgumentError AA.dctplan(randn(64))            # no Float64 DCT in vDSP
    @test_throws ArgumentError AA.dctplan(Float64, 64)
    @test_throws ArgumentError AA.dctplan(randn(Float32, 64), 1)
    @test_throws ArgumentError AA.dctplan(randn(Float32, 24))   # k < 4
end

@testset "mul! does not allocate" begin
    for T in (Float32, Float64)
        for n in (64, 120, 8192)                       # interleaved (2^k, mixed), stride-2 zop
            x = randn(Complex{T}, n); y = similar(x)
            @test _plan_allocs(y, AA.fftplan(x), x) == 0
            @test _plan_allocs(y, AA.ifftplan(x), x) == 0
            @test _plan_allocs(x, AA.fftplan(x), x) == 0          # in place
        end
        x = randn(Complex{T}, 16, 64); y = similar(x)
        @test _plan_allocs(y, AA.fftplan(x), x) == 0
        @test _plan_allocs(y, AA.fftplan(x, 1), x) == 0
        @test _plan_allocs(y, AA.ifftplan(x, 2), x) == 0
        r = randn(T, 1024); R = Vector{Complex{T}}(undef, 513)
        @test _plan_allocs(R, AA.rfftplan(r), r) == 0
        @test _plan_allocs(r, AA.irfftplan(R, 1024), R) == 0
    end
    x = randn(Float32, 256); y = similar(x)
    @test _plan_allocs(y, AA.dctplan(x), x) == 0
end

@testset "introspection and show" begin
    x = randn(ComplexF64, 8, 32)
    p = AA.fftplan(x)
    @test p isa AA.FFTPlan
    @test size(p) == (8, 32) && size(p, 1) == 8 && size(p, 3) == 1
    @test ndims(p) == 2 && length(p) == 256
    @test eltype(p) === ComplexF64
    @test AA.output_size(p) == (8, 32)
    @test occursin("forward 2-D FFT", sprint(show, p))
    @test occursin("dims=2", sprint(show, AA.bfftplan(x, 2)))
    @test occursin("scaled by", sprint(show, AA.ifftplan(x)))

    r = randn(Float32, 64)
    pr = AA.rfftplan(r)
    @test size(pr) == (64,) && AA.output_size(pr) == (33,)
    @test eltype(pr) === Float32
    @test eltype(inv(pr)) === ComplexF32
    @test size(inv(pr)) == (33,) && AA.output_size(inv(pr)) == (64,)
    @test occursin("forward real FFT", sprint(show, pr))
    @test occursin("vDSP_DFT_Interleaved", sprint(show, AA.fftplan(ComplexF32, 64)))
    @test occursin("DCT-III", sprint(show, inv(AA.dctplan(Float32, 64))))
end

@testset "argument checking" begin
    @test_throws ArgumentError AA.fftplan(randn(ComplexF64, 7))          # unsupported length
    @test_throws ArgumentError AA.fftplan(randn(ComplexF64, 8, 12))      # 2-D must be 2^k
    @test_throws ArgumentError AA.fftplan(randn(ComplexF64, 12, 8), 1)
    @test_throws ArgumentError AA.fftplan(randn(ComplexF64, 8, 8), 3)
    @test_throws ArgumentError AA.fftplan(randn(ComplexF64, 8), 2)
    @test_throws ArgumentError AA.fftplan(ComplexF64, (4, 4, 4))
    @test_throws ArgumentError AA.rfftplan(randn(7))
    @test_throws ArgumentError AA.rfftplan(randn(1))
    @test_throws ArgumentError AA.rfftplan(randn(24))                    # complex-only length
    @test_throws ArgumentError AA.rfftplan(randn(8, 12))
    @test_throws DimensionMismatch AA.brfftplan(randn(ComplexF64, 33), 128)
    msg = try AA.fftplan(randn(ComplexF64, 7)); "" catch err; sprint(showerror, err) end
    @test occursin("is_supported_fft_length", msg)

    x = randn(ComplexF64, 64); p = AA.fftplan(x)
    @test_throws DimensionMismatch p * randn(ComplexF64, 32)
    @test_throws DimensionMismatch mul!(zeros(ComplexF64, 32), p, x)
    @test_throws ArgumentError p * randn(ComplexF32, 64)                 # wrong precision
    @test_throws ArgumentError mul!(zeros(ComplexF32, 64), p, x)
    @test_throws ArgumentError p * view(randn(ComplexF64, 128), 1:2:128) # not contiguous

    # Contiguous views are fine (roots the parent across the ccall).
    M = randn(ComplexF64, 64, 3)
    @test p * view(M, :, 2) ≈ FFTW.fft(M[:, 2])
    out = zeros(ComplexF64, 64, 3)
    mul!(view(out, :, 3), p, view(M, :, 2))
    @test out[:, 3] ≈ FFTW.fft(M[:, 2])
    @test all(iszero, out[:, 1:2])
end

@testset "one plan, many tasks" begin
    for x in (randn(ComplexF64, 1024), randn(ComplexF64, 8192), randn(ComplexF64, 32, 64))
        p = AA.fftplan(x)
        ref = FFTW.fft(x)
        xs = [x .* k for k in 1:16]
        ys = Vector{typeof(x)}(undef, 16)
        @sync for k in 1:16
            Threads.@spawn (ys[k] = p * xs[k])
        end
        @test all(ys[k] ≈ ref .* k for k in 1:16)
    end
    # Concurrent plan creation goes through the shared, locked setup caches.
    ps = Vector{Any}(undef, 16)
    @sync for k in 1:16
        Threads.@spawn (ps[k] = AA.fftplan(ComplexF32, 2048))
    end
    @test all(p.setup === ps[1].setup for p in ps)
end

@testset "fallback keyword on the one-shot API" begin
    x = randn(ComplexF64, 7)
    @test_throws ArgumentError AA.fft(x)
    @test AA.fft(x; fallback = FFTW.fft) ≈ FFTW.fft(x)
    @test AA.bfft(x; fallback = FFTW.bfft) ≈ FFTW.bfft(x)
    @test AA.ifft(x; fallback = FFTW.ifft) ≈ FFTW.ifft(x)
    r = randn(7)
    @test_throws ArgumentError AA.rfft(r)
    @test AA.rfft(r; fallback = FFTW.rfft) ≈ FFTW.rfft(r)
    r24 = randn(24)                  # fine for the complex DFT, rejected by the real one
    @test_throws ArgumentError AA.rfft(r24)
    @test AA.rfft(r24; fallback = FFTW.rfft) ≈ FFTW.rfft(r24)

    # The fallback is never consulted when vDSP can do the transform.
    boom(_) = error("fallback must not be called")
    for n in (64, 120)
        y = randn(ComplexF64, n)
        @test AA.fft(y; fallback = boom) ≈ FFTW.fft(y)
        @test AA.ifft(y; fallback = boom) ≈ FFTW.ifft(y)
    end
    @test AA.rfft(randn(64); fallback = boom) isa Vector{ComplexF64}
    @test AA.rfft(randn(48); fallback = boom) isa Vector{ComplexF64}
end

end # @testset "FFTPlan"
