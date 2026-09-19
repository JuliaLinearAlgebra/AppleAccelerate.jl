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
# This has to run in a child process, because two things the *test harness* does each
# independently stop the loop from vectorising at all:
#
#   * `Pkg.test()` defaults to `--check-bounds=yes`, which disables `@inbounds`;
#   * `julia-actions/julia-runtest` defaults to coverage on, and the counters
#     coverage inserts into the loop body block vectorisation outright.
#
# Neither affects correctness, so the numeric tests above stay green and only this
# check would fail -- in CI only, which is exactly where a real regression needs
# catching. So spawn a child with both explicitly turned off. `Base.julia_cmd()`
# propagates the parent's `--check-bounds` and `--code-coverage`, so both have to be
# overridden here; the last occurrence of a flag wins.
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
# `@simdmath` must reach the same SIMD routines, through its dispatchers, as the
# scoped-import spelling -- and then produce the same bits, since it is the same call.
macro_loop(o, x, y) = (SM.@simdmath @simd for i in eachindex(x, y, o)
    @inbounds o[i] = log(x[i])^y[i] + atan(y[i], x[i]) + x[i]^2 + rem(i, 4)
end; o)
inner_loop(o, x, y) = (@simd for i in eachindex(x, y, o)
    SM.@simdmath @inbounds o[i] = log(x[i])^y[i] + atan(y[i], x[i]) + x[i]^2 + rem(i, 4)
end; o)
plain_loop(o, x, y) = (@simd for i in eachindex(x, y, o)
    @inbounds o[i] = SM.pow(SM.log(x[i]), y[i]) + SM.atan(y[i], x[i]) + x[i]^2 + rem(i, 4)
end; o)
for (T, suffix) in ((Float64, "d2"), (Float32, "f4"))
    for (label, g) in (("@simdmath outside @simd", macro_loop), ("@simdmath inside @simd", inner_loop))
        io = IOBuffer()
        code_native(io, g, ntuple(_ -> Vector{T}, 3); debuginfo = :none)
        asm = String(take!(io))
        for base in ("log", "pow", "atan2")
            occursin(Regex("\\b__simd_$(base)_$(suffix)\\b"), asm) || push!(failures, "$label/$T -> _simd_$(base)_$(suffix)")
        end
        x = rand(T, 1003) .+ T(1.5); y = rand(T, 1003)
        g(similar(x), x, y) == plain_loop(similar(x), x, y) || push!(failures, "$label/$T differs from scoped calls")
    end
end
# Strided access vectorises only when the stride is a compile-time constant. A stride
# known only at run time (`for i in 1:stride:length(x)` with `stride::Int`) leaves the
# call scalar on every LLVM tried so far; that is a vectoriser limitation, not a broken
# mapping, so it is deliberately not asserted either way here.
const_stride(x) = (u = zero(eltype(x)); @simd for i in 1:2:length(x); @inbounds u += SM.log(x[i]); end; u)
for (T, suffix) in ((Float64, "d2"), (Float32, "f4"))
    io = IOBuffer()
    code_native(io, const_stride, (Vector{T},); debuginfo = :none)
    occursin(Regex("\\b__simd_log_$(suffix)\\b"), String(take!(io))) || push!(failures, "constant stride/$T -> _simd_log_$(suffix)")
end
print(join(failures, ","))
"""

@testset "SIMDMath" begin

    @testset "vectorises to the SIMD routine" begin
        script = tempname() * ".jl"
        write(script, _VECTORISE_CHECK)
        try
            out = read(`$(Base.julia_cmd()) --check-bounds=auto --code-coverage=none
                        --startup-file=no --project=$(Base.active_project()) $script`, String)
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

    # No Base equivalent, so pin against known values instead of a reference impl.
    @testset "remainder / nextafter" begin
        for T in (Float32, Float64)
            @test SM.nextafter(one(T), T(2)) === nextfloat(one(T))
            @test SM.nextafter(one(T), zero(T)) === prevfloat(one(T))
            @test SM.remainder(T(5), T(3)) ≈ T(-1)
        end
    end

    # Some functions are intentionally not provided; see the note in src/simdmath.jl.
    # Assert on the tables rather than isdefined(SM, f): every module implicitly does
    # `using Base`, so isdefined(SM, :sqrt) is true regardless of what we define.
    @testset "unaccelerated functions stay out" begin
        provided = Set(Symbol[first(t) for t in (SM.UNARY..., SM.BINARY...)])
        # measured at ~0.98x of a scalar libm loop -- vector signature, scalar guts
        @test :erf ∉ provided
        @test :erfc ∉ provided
        @test :tgamma ∉ provided
        # writes the global signgam, so it cannot be declared memory(none)
        @test :lgamma ∉ provided
        # LLVM lowers these to native instructions; a libcall would be slower
        for f in (:sqrt, :floor, :ceil, :round, :trunc, :fma, :abs, :fabs)
            @test f ∉ provided
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

    # Every table entry names two C symbols; a typo in either only shows up as a
    # link error at first call, or (for the SIMD one) as a silent loss of speed.
    @testset "table symbols resolve" begin
        libm = AppleAccelerate.Libdl.dlopen("/usr/lib/system/libsystem_m.dylib")
        for tbl in (SM.UNARY, SM.BINARY), (_jl, c64, c32, simd) in tbl
            for sym in (c64, c32, "_simd_$(simd)_d2", "_simd_$(simd)_f4")
                @test AppleAccelerate.Libdl.dlsym_e(libm, sym) != C_NULL
            end
        end
    end

    @testset "@simdmath" begin
        rewritten(ex) = SM._rewrite(ex)
        disp(f) = GlobalRef(SM, SM._dispatcher(f))

        @testset "what is rewritten" begin
            @test rewritten(:(log(x))) == Expr(:call, disp(:log), :x)
            @test rewritten(:(atan(y, x))) == Expr(:call, disp(:atan), :y, :x)
            @test rewritten(:(atan(x))) == Expr(:call, disp(:atan), :x)
            @test rewritten(:(x^y)) == Expr(:call, disp(:pow), :x, :y)
            # nested, and through the macros it has to compose with
            @test rewritten(:(exp(log(x)))) == Expr(:call, disp(:exp), Expr(:call, disp(:log), :x))
            loop = quote
                @simd for i in r
                    @inbounds u += log(X[i])
                end
            end
            @test occursin("_simdmath_log", string(rewritten(loop)))
        end

        @testset "what is left alone" begin
            for ex in (:(Base.log(x)), :(log.(x)), :(x .^ y), :(x^2), :(foo(x)), :(sqrt(x)),
                       :(log(x; base = 2)), :(log(xs...)), :(atan(a, b, c)), :(hypot(x)),
                       :(erf(x)), :(:(log(x))), :(quote log(x) end))
                @test rewritten(ex) == ex
            end
            # the signature of a method defined inside the block is a binding, not a call
            @test rewritten(:(log(x) = exp(x))).args[1] == :(log(x))
            @test rewritten(:(function log(x); exp(x); end)).args[1] == :(log(x))
            @test occursin("_simdmath_exp", string(rewritten(:(log(x) = exp(x)))))
        end

        @testset "values" begin
            for T in (Float32, Float64)
                x, y = T(1.7), T(0.3)
                @test (SM.@simdmath log(x)) === SM.log(x)
                @test (SM.@simdmath x^y) === SM.pow(x, y)
                @test (SM.@simdmath atan(y, x)) === SM.atan(y, x)
                @test (SM.@simdmath x^2) === x^2
                @test @inferred(SM._simdmath_log(x)) isa T
                @test @inferred(SM._simdmath_pow(x, y)) isa T
            end
            # Types SIMDMath has no method for keep their Base meaning instead of
            # throwing, which is what makes the macro safe around index arithmetic.
            n = 3
            @test (SM.@simdmath rem(7, 4)) === 3
            @test (SM.@simdmath 2.0^n) === 8.0
            @test (SM.@simdmath 2^n) === 8
            @test (SM.@simdmath log(1)) === 0.0
            @test (SM.@simdmath exp(1.0im)) == exp(1.0im)
            @test (SM.@simdmath hypot(3.0f0, 4.0)) === 5.0
            @test (SM.@simdmath log(big"2.0")) == log(big"2.0")
            # no Base counterpart to fall back to
            @test_throws MethodError SM.@simdmath nextafter(1, 2)
        end

        # The macro must not leak or capture names: it works from a module that has
        # never heard of SIMDMath, and leaves that module's own `log` alone.
        @testset "hygiene" begin
            m = Module()
            Core.eval(m, :(using AppleAccelerate))
            Core.eval(m, :(f(x) = AppleAccelerate.SIMDMath.@simdmath log(x) + exp(x)))
            @test Base.invokelatest(m.f, 2.0) == SM.log(2.0) + SM.exp(2.0)
            @test Core.eval(m, :(log)) === Base.log
        end

        @testset "loop matches explicit SIMDMath calls" begin
            for T in (Float32, Float64), n in (1, 3, 8, 33)
                X = rand(T, n) .+ T(1.5); Y = rand(T, n)
                A = similar(X); B = similar(X)
                SM.@simdmath @simd for i in eachindex(X, Y, A)
                    @inbounds A[i] = log(X[i])^Y[i] + hypot(X[i], Y[i])
                end
                @simd for i in eachindex(X, Y, B)
                    @inbounds B[i] = SM.pow(SM.log(X[i]), Y[i]) + SM.hypot(X[i], Y[i])
                end
                @test A == B
            end
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
