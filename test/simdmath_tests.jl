const SM = AppleAccelerate.SIMDMath

# Inputs inside each function's real domain, so the reference never returns NaN and
# the comparison actually tests something.
function _simd_domain(S::Symbol, ::Type{T}, n::Int) where {T}
    r = rand(T, n)
    S in (:acos, :asin, :atanh)                ? (r .* T(1.8) .- T(0.9)) :
    S in (:acosh,)                             ? (r .* T(4) .+ T(1.1)) :
    S in (:log, :log2, :log10, :cbrt, :tgamma) ? (r .* T(8) .+ T(0.3)) :
    S in (:log1p,)                             ? (r .* T(4)) :
    S in (:exp, :exp2, :exp10)                 ? (r .* T(4) .- T(2)) :
    S in (:sinpi, :cospi, :tanpi)              ? (r .* T(3) .- T(1.5)) :
                                                 (r .* T(2) .- T(1))
end

_ulps(got::T, ref::T) where {T} = got == ref ? 0.0 : Float64(abs(big(got) - big(ref)) / eps(ref))

# Deliberately loose: these routines trade accuracy for speed, and the worst observed
# on macOS 26 / M-series was 3 ULP (tanpi, Float64). The bound exists to catch a
# broken *mapping* in the UNARY / BINARY tables, not to pin an accuracy contract.
const ULP_TOL = 8.0

# The whole point of this module is that the scalar call is replaced by a call to the
# SIMD routine. If that silently stops happening the functions still return correct
# answers, just slowly -- so assert on the emitted code, not only on the numbers.
#
# This has to run in a child process. Vectorisation only happens when `@inbounds` is
# honoured, and `Pkg.test()` defaults to `--check-bounds=yes`, which disables it and
# blocks vectorisation outright. Checking in-process would mean this test silently
# never runs under CI -- precisely where a regression would slip through -- so spawn
# a child with the default bounds-checking setting instead.
const _VECTORISE_CHECK = raw"""
using AppleAccelerate, InteractiveUtils
const SM = AppleAccelerate.SIMDMath
failures = String[]
for (tbl, nargs) in ((SM.UNARY, 1), (SM.BINARY, 2))
    for (jlname, _c64, _c32, simd) in tbl
        f = getfield(SM, jlname)
        for (T, suffix) in ((Float64, "d2"), (Float32, "f4"))
            sym = "_simd_$(simd)_$(suffix)"
            g = nargs == 1 ?
                (o, x) -> (@simd for i in eachindex(x, o); @inbounds o[i] = f(x[i]); end; o) :
                (o, x, y) -> (@simd for i in eachindex(x, o); @inbounds o[i] = f(x[i], y[i]); end; o)
            argtypes = ntuple(_ -> Vector{T}, nargs + 1)
            io = IOBuffer()
            code_native(io, g, argtypes; debuginfo = :none)
            occursin(Regex("\\b_" * sym * "\\b"), String(take!(io))) || push!(failures, "$jlname/$T -> $sym")
        end
    end
end
print(join(failures, ","))
"""

@testset "SIMDMath" begin

    @testset "vectorises to the SIMD routine" begin
        script = tempname() * ".jl"
        write(script, _VECTORISE_CHECK)
        try
            out = read(`$(Base.julia_cmd()) --check-bounds=auto --startup-file=no
                        --project=$(Base.active_project()) $script`, String)
            failed = filter(!isempty, split(out, ","))
            # Name the offenders rather than just asserting a count, so a regression
            # points straight at the broken mapping.
            @test failed == SubString{String}[]
        finally
            rm(script; force = true)
        end
    end

    # Cross-check against Base. This is what catches a wrong entry in the UNARY /
    # BINARY tables -- a mismapped name or (for atan2) a swapped argument order
    # produces a completely different number, not a few-ULP difference.
    @testset "matches Base (1-arg)" begin
        for (jlname, basef) in ((:acos, acos), (:acosh, acosh), (:asin, asin), (:asinh, asinh),
                                (:atan, atan), (:atanh, atanh), (:cbrt, cbrt), (:cos, cos),
                                (:cosh, cosh), (:cospi, cospi), (:exp, exp), (:exp10, exp10),
                                (:exp2, exp2), (:expm1, expm1), (:log, log), (:log10, log10),
                                (:log1p, log1p), (:log2, log2), (:sin, sin), (:sinh, sinh),
                                (:sinpi, sinpi), (:tan, tan), (:tanh, tanh), (:tanpi, tanpi))
            f = getfield(SM, jlname)
            for T in (Float32, Float64)
                X = _simd_domain(jlname, T, N)
                worst = maximum(eachindex(X)) do i
                    ref = basef(X[i])
                    (isfinite(ref) && ref != 0) ? _ulps(f(X[i]), ref) : 0.0
                end
                @test worst <= ULP_TOL
            end
        end
    end

    @testset "matches Base (2-arg)" begin
        for T in (Float32, Float64)
            X = _simd_domain(:atan, T, N)
            Y = abs.(_simd_domain(:atan, T, N)) .+ T(0.5)
            # atan(y, x) must keep Base's argument order, not C's atan2 order by accident
            @test maximum(i -> _ulps(SM.atan(X[i], Y[i]), atan(X[i], Y[i])), eachindex(X)) <= ULP_TOL
            @test maximum(i -> _ulps(SM.hypot(X[i], Y[i]), hypot(X[i], Y[i])), eachindex(X)) <= ULP_TOL
            @test maximum(i -> _ulps(SM.rem(X[i], Y[i]), rem(X[i], Y[i])), eachindex(X)) <= ULP_TOL

            B = abs.(rand(T, N)) .+ T(0.5)
            E = rand(T, N) .* T(3)
            @test maximum(i -> _ulps(SM.pow(B[i], E[i]), B[i]^E[i]), eachindex(B)) <= ULP_TOL
        end
    end

    # No Base equivalent, so pin against the identities instead of a reference impl.
    @testset "erf / erfc / tgamma / remainder / nextafter" begin
        for T in (Float32, Float64)
            X = _simd_domain(:erf, T, 256)
            @test all(i -> SM.erf(X[i]) + SM.erfc(X[i]) ≈ one(T), eachindex(X))
            @test SM.erf(zero(T)) == zero(T)
            # tgamma(n+1) == n!
            for n in 1:8
                @test SM.tgamma(T(n + 1)) ≈ T(factorial(n)) rtol = sqrt(eps(T))
            end
            @test SM.nextafter(one(T), T(2)) === nextfloat(one(T))
            @test SM.nextafter(one(T), zero(T)) === prevfloat(one(T))
            @test SM.remainder(T(5), T(3)) ≈ T(-1)
        end
    end

    # Outside a vectorised loop these must still be correct -- that path calls the
    # scalar libm symbol, which is a different code path from the SIMD routine.
    @testset "scalar path (unvectorised) is correct" begin
        for T in (Float32, Float64)
            @test SM.log(T(1)) == zero(T)
            @test SM.exp(zero(T)) == one(T)
            @test SM.sin(zero(T)) == zero(T)
            @test SM.cos(zero(T)) == one(T)
            @test SM.log(SM.exp(T(2))) ≈ T(2)
            @test SM.sinpi(T(1)) ≈ zero(T) atol = eps(T)
            @test SM.cospi(zero(T)) == one(T)
            @test SM.pow(T(2), T(10)) ≈ T(1024)
            @test SM.hypot(T(3), T(4)) ≈ T(5)
        end
    end

    # A vectorised loop and its scalar remainder take different code paths; a length
    # that is not a multiple of the vector width exercises both in one call.
    @testset "vector body and scalar remainder agree" begin
        for T in (Float32, Float64), n in (1, 2, 3, 5, 7, 8, 15, 31, 33)
            X = abs.(rand(T, n)) .+ T(0.5)
            O = similar(X)
            @simd for i in eachindex(X, O)
                @inbounds O[i] = SM.log(X[i])
            end
            @test all(i -> _ulps(O[i], log(X[i])) <= ULP_TOL, eachindex(X))
        end
    end

    @testset "type stability" begin
        for T in (Float32, Float64)
            @test @inferred(SM.log(one(T))) isa T
            @test @inferred(SM.pow(T(2), T(3))) isa T
            @test @inferred(SM.atan(one(T), one(T))) isa T
        end
        # Float32 in must not silently widen to Float64
        @test SM.log(1.5f0) isa Float32
        @test SM.log(1.5) isa Float64
    end
end
