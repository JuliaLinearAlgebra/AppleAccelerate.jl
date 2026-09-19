## simdmath.jl — scalar math functions that auto-vectorise inside `@simd` loops ##
#
# `array.jl` exposes vForce's `vv*` routines, which are *array*-shaped: they need a
# materialised input array and somewhere to put the output. That shape does not fit
# every loop — strided access, values computed on the fly, control flow in the loop
# body, or any iterable that is not a dense `Array`. In those cases the only option
# today is a scalar `Base` loop.
#
# LLVM's loop vectoriser can replace a *scalar* call inside a vectorising loop with
# a call to an equivalent SIMD routine, provided it is told which routine to use.
# `clang` spells that `-fveclib=Accelerate`; the underlying mechanism is the LLVM
# "vector function ABI" (VFABI), which is driven entirely by IR and so is reachable
# from Julia via `llvmcall` — no compiler flag, no patched Julia. This file uses it
# to expose Apple's SIMD math routines as ordinary scalar Julia functions that
# vectorise when (and only when) the surrounding loop does.
#
# ## Scope: this does NOT replace the array API
#
# Measured on an M-series machine, vForce's `vv*` routines beat these SIMD routines
# by roughly 2x per element at *every* array size tested (N = 16 … 500_000) — there
# is no small-N crossover, and vForce still wins even when it needs extra temporary
# arrays and extra passes over memory. `_simd_log_f4` is simply not as well tuned as
# `vvlogf`, which can unroll and pipeline across a long array in a way a single
# 4-lane call cannot.
#
# So: **if the data is a dense array and you can call `AppleAccelerate.log!`, do
# that instead — it is faster.** These functions are for the cases where that is not
# an option, where the honest comparison is against a scalar `Base` loop, and where
# they run about 2–4x faster (Float32) or 1.2–2x faster (Float64).
#
# ## How the generated IR works
#
# For each function three things must be true, or the mapping is silently ignored
# and you get a scalar libm call back:
#
#   1. The scalar call must be to a *named symbol* (`@logf`), not an opaque Julia
#      function, so the vectoriser has something to match on.
#   2. The vector routine (`@_simd_log_f4`) must be *declared in the module*.
#      `VFDatabase` looks it up by name and drops the mapping if it is absent.
#   3. That declaration must be anchored in `@llvm.compiler.used`. Without it the
#      declaration is unreferenced at the point the vectoriser runs, gets stripped
#      by earlier passes, rule 2 fails, and everything silently falls back to
#      scalar. This is the single easiest step to miss.
#
# The mapping itself is the `"vector-function-abi-variant"` attribute on the call,
# whose value is a VFABI-mangled string: `_ZGV_LLVM_N4v_logf(_simd_log_f4)` reads
# as "the `_LLVM_` ISA, Not masked, 4 lanes, one vector operand, mapping scalar
# `logf` onto vector `_simd_log_f4`".
#
# ## Why `_simd_*` and not Accelerate's `vfp.h`
#
# Accelerate's own `vfp.h` SIMD entry points (`vlogf`, `vexpf`, …) are Float32-only
# and cover 22 functions. libsystem_m's `<simd/math.h>` (`_simd_log_f4`,
# `_simd_log_d2`, …) is public SDK API, covers the same ground plus Float64, and is
# what recent LLVM prefers on Darwin. Since Float64 is Julia's default, using one
# source for both widths keeps accuracy consistent across `Float32`/`Float64`
# rather than mixing two vendors' implementations.
#
# ## Accuracy
#
# These routines are *less accurate than Base*. See the `SIMDMath` docstring; the
# tolerances are pinned in `test/simdmath_tests.jl`.
#
# ## Why the package needs Julia >= 1.11
#
# The IR below is LLVM 16+. Julia 1.11 ships LLVM 16, Julia 1.10 ships LLVM 15, and
# two things here do not exist in 15:
#
#   * opaque pointers (`ptr`). On LLVM 15 `@llvm.compiler.used` has to be spelled
#     `[1 x i8*]` with an explicit `bitcast` of the typed function pointer.
#   * the `memory(none)` attribute, which was `readnone` before LLVM 16.
#
# Emitting both dialects is possible but doubles the surface for one superseded
# LLVM, so `julia = "1.11"` is the floor instead. On 1.10 this file fails to parse
# with `error: expected type` at the `ptr` in `[1 x ptr]`.

"""
    AppleAccelerate.SIMDMath

Scalar math functions that LLVM replaces with Apple's SIMD math routines when they
appear inside a vectorised loop.

Unlike the array functions in `AppleAccelerate` proper (which wrap vForce and need
a whole array to work on), these are *scalar* functions called from inside a `@simd`
loop:

```jldoctest simdmath
julia> using AppleAccelerate.SIMDMath: log

julia> function weighted_logsum(X, W)  # no temporary for log.(X), unlike the array API
           u = zero(eltype(X))
           @simd for i in eachindex(X, W)
               @inbounds u += W[i] * log(X[i])
           end
           u
       end;

julia> X = collect(1.0:1000.0); W = fill(0.5, 1000);

julia> weighted_logsum(X, W) ≈ sum(W .* Base.log.(X))
true
```

# When to use this — and when not to

!!! important
    **This is not a faster replacement for the array API.** vForce's `vv*` routines
    (`AppleAccelerate.log!` and friends) are about 2x faster per element than these
    SIMD routines at every array size measured, and win even when they need extra
    temporaries and extra passes over memory. If your data is a dense array, use
    those instead.

`SIMDMath` is for the cases where vForce does not apply:

  * strided or otherwise non-contiguous access;
  * values computed on the fly rather than read from an array;
  * loops with control flow, or over iterables that are not dense `Array`s;
  * hot paths where you cannot allocate and cannot preallocate a temporary.

There, the realistic alternative is a scalar `Base` loop, against which these run
roughly 2–4x faster for `Float32` and 1.2–2x for `Float64`.

# Accuracy

!!! warning
    These routines trade accuracy for speed and are **less accurate than the
    equivalent `Base` functions**, which are correctly rounded (≤ 0.5 ULP).
    Measured worst case is roughly 1–3 ULP depending on the function; see
    `test/simdmath_tests.jl` for the pinned per-function bounds.

    Accuracy is also **not consistent across the same call site**. If the loop
    vectorises you get the SIMD routine; if it does not (an unvectorisable loop,
    or the scalar remainder of one that did) you get the correctly-rounded libm
    call instead. Do not use these where reproducibility across loop shapes
    matters.

# Vectorising

A call only becomes a SIMD call if the surrounding loop actually vectorises, which
in practice means `@simd` (and usually `@inbounds`). Outside a vectorised loop
these are ordinary — correct, but no faster than `Base`, and sometimes slower.
`Float32` runs 4 lanes at a time, `Float64` runs 2.

!!! note
    `--check-bounds=yes` (which `Pkg.test()` passes by default) and `--code-coverage`
    each stop the loop vectorising at all, silently reducing these to scalar libm
    calls. Results stay correct; the speedup disappears.

To confirm you are getting the vectorised form, look for the symbol:

```julia
using InteractiveUtils
@code_native weighted_logsum(X, W)   # expect a call to _simd_log_d2
```

# Available functions

One argument: `acos`, `acosh`, `asin`, `asinh`, `atan`, `atanh`, `cbrt`, `cos`,
`cosh`, `cospi`, `exp`, `exp10`, `exp2`, `expm1`, `log`, `log10`, `log1p`, `log2`,
`sin`, `sinh`, `sinpi`, `tan`, `tanh`, `tanpi`.

Two arguments: `atan(y, x)`, `hypot`, `nextafter`, `pow`, `rem`, `remainder`.

All are defined for `Float32` and `Float64` only. `nextafter`, `pow` and `remainder`
have no `Base` counterpart with matching semantics and are named after their C
equivalents; use `pow(x, y)` rather than `x^y` (or let [`@simdmath`](@ref) rewrite
`x^y` for you).

These are all the one- and two-argument routines `<simd/math.h>` provides that are
worth calling. Names such as `copysign`, `min`/`max` (`fmin`/`fmax`) and `fdim` have
no `_simd_*` entry point at all -- the header implements them inline, and LLVM
already vectorises the `Base` versions natively.

Deliberately absent:

  * `sqrt`, `floor`, `ceil`, `round`, `trunc`, `fma`, `abs` — LLVM already lowers
    these to native instructions, so a library call would be *slower*.
  * `erf`, `erfc`, `tgamma` — `<simd/math.h>` declares these, but they measure at
    0.97–0.99x of a scalar libm loop: they are scalar loops behind a vector
    signature and buy nothing.
  * `lgamma` — writes the global `signgam`, so it cannot be declared `memory(none)`.
"""
module SIMDMath

# Functions LLVM lowers natively (sqrt/floor/ceil/round/trunc/fma/fabs) are
# excluded on purpose: a libcall would be slower than the instruction.

# (Julia name, C scalar symbol for Float64, C scalar symbol for Float32,
#  `<simd/math.h>` base name -- `_simd_<base>_d2` / `_simd_<base>_f4`)
#
# The C scalar name is what the *unvectorised* path calls, so it has to be the real
# libm symbol. It is usually `<base>` for Float64 and `<base>f` for Float32, but
# Darwin spells a few with a leading underscore pair, hence the explicit columns.
const UNARY = (
    # (jlname,   c64,        c32,         simd base)
    (:acos,     "acos",     "acosf",     "acos"),
    (:acosh,    "acosh",    "acoshf",    "acosh"),
    (:asin,     "asin",     "asinf",     "asin"),
    (:asinh,    "asinh",    "asinhf",    "asinh"),
    (:atan,     "atan",     "atanf",     "atan"),
    (:atanh,    "atanh",    "atanhf",    "atanh"),
    (:cbrt,     "cbrt",     "cbrtf",     "cbrt"),
    (:cos,      "cos",      "cosf",      "cos"),
    (:cosh,     "cosh",     "coshf",     "cosh"),
    (:cospi,    "__cospi",  "__cospif",  "cospi"),
    (:exp,      "exp",      "expf",      "exp"),
    (:exp10,    "__exp10",  "__exp10f",  "exp10"),
    (:exp2,     "exp2",     "exp2f",     "exp2"),
    (:expm1,    "expm1",    "expm1f",    "expm1"),
    (:log,      "log",      "logf",      "log"),
    (:log10,    "log10",    "log10f",    "log10"),
    (:log1p,    "log1p",    "log1pf",    "log1p"),
    (:log2,     "log2",     "log2f",     "log2"),
    (:sin,      "sin",      "sinf",      "sin"),
    (:sinh,     "sinh",     "sinhf",     "sinh"),
    (:sinpi,    "__sinpi",  "__sinpif",  "sinpi"),
    (:tan,      "tan",      "tanf",      "tan"),
    (:tanh,     "tanh",     "tanhf",     "tanh"),
    (:tanpi,    "__tanpi",  "__tanpif",  "tanpi"),
)
# `lgamma` is deliberately omitted: it writes the global `signgam`, which the
# `memory(none)` declaration below would let LLVM assume never happens.
#
# `erf`, `erfc` and `tgamma` are omitted despite `<simd/math.h>` declaring
# `_simd_erf_f4` and friends: measured against a plain scalar libm loop they come
# out at 0.97-0.99x, i.e. those routines are scalar loops behind a vector signature
# and buy nothing. Exposing them would advertise an acceleration that is not there.

const BINARY = (
    # (jlname,     c64,          c32,           simd base)
    (:atan,       "atan2",      "atan2f",      "atan2"),
    (:hypot,      "hypot",      "hypotf",      "hypot"),
    (:nextafter,  "nextafter",  "nextafterf",  "nextafter"),
    (:pow,        "pow",        "powf",        "pow"),
    (:rem,        "fmod",       "fmodf",       "fmod"),
    (:remainder,  "remainder",  "remainderf",  "remainder"),
)

"""
    _ir(cname, vecname, W, lt, nargs) -> String

Build the `llvmcall` module for one function: a scalar call to `cname` carrying a
VFABI mapping onto the `W`-lane `vecname`, with `vecname` anchored in
`@llvm.compiler.used` so it survives to the loop vectoriser.
"""
function _ir(cname::String, vecname::String, W::Int, lt::String, nargs::Int)
    formals = join(("$lt %$(i-1)" for i in 1:nargs), ", ")
    actuals = formals
    scalar_tys = join((lt for _ in 1:nargs), ", ")
    vector_tys = join(("<$W x $lt>" for _ in 1:nargs), ", ")
    vtok = "v"^nargs
    """
    @llvm.compiler.used = appending global [1 x ptr] [ptr @"$vecname"], section "llvm.metadata"
    declare <$W x $lt> @"$vecname"($vector_tys)
    declare $lt @"$cname"($scalar_tys) #0
    define $lt @entry($formals) #2 {
    top:
      %r = call $lt @"$cname"($actuals) #1
      ret $lt %r
    }
    attributes #0 = { nounwind willreturn memory(none) }
    attributes #1 = { "vector-function-abi-variant"="_ZGV_LLVM_N$W$(vtok)_$cname($vecname)" }
    attributes #2 = { alwaysinline nounwind willreturn memory(none) }
    """
end

# The scalar declaration is `memory(none)`, i.e. errno is assumed untouched. Julia
# never reads math errno, and this is the same assumption clang makes under
# -fno-math-errno; without it the calls could not be hoisted or vectorised at all.

for (tbl, nargs) in ((UNARY, 1), (BINARY, 2))
    for (jlname, c64, c32, simd) in tbl
        for (T, lt, cname, W, suffix) in ((Float64, "double", c64, 2, "d2"),
                                          (Float32, "float",  c32, 4, "f4"))
            ir = _ir(cname, "_simd_$(simd)_$(suffix)", W, lt, nargs)
            args = [Symbol(:x, i) for i in 1:nargs]
            sig = [:($(args[i])::$T) for i in 1:nargs]
            argtys = [T for _ in 1:nargs]
            @eval @inline function $jlname($(sig...))
                Base.llvmcall(($ir, "entry"), $T, Tuple{$(argtys...)}, $(args...))
            end
        end
    end
end

# ---------------------------------------------------------------------------------
# `@simdmath`
# ---------------------------------------------------------------------------------
#
# `using AppleAccelerate.SIMDMath: log` shadows `Base.log` for the *whole* enclosing
# module, which is a lot of blast radius for one hot loop -- and because the functions
# above only have `Float32`/`Float64` methods, every other `log(::Int)` or
# `rem(i, n)` in that module then throws a `MethodError`. The macro scopes the
# substitution to one expression instead, and routes through the `_simdmath_*`
# dispatchers below so that argument types SIMDMath does not cover keep their `Base`
# meaning.

# Julia name => arities provided. `atan` appears in both tables.
const _ARITIES = let d = Dict{Symbol,Vector{Int}}()
    for (tbl, nargs) in ((UNARY, 1), (BINARY, 2)), t in tbl
        push!(get!(d, first(t), Int[]), nargs)
    end
    d
end

_dispatcher(name::Symbol) = Symbol(:_simdmath_, name)

for (jlname, arities) in _ARITIES
    disp = _dispatcher(jlname)
    for nargs in arities
        args = [Symbol(:x, i) for i in 1:nargs]
        sig = [:($(a)::T) for a in args]
        @eval @inline $disp($(sig...)) where {T<:Union{Float32,Float64}} = $jlname($(args...))
    end
    # Anything else keeps its ordinary meaning. `pow` is what `x^y` rewrites to, so
    # its fallback is `^`; `nextafter` and `remainder` have no `Base` counterpart and
    # so get no fallback -- a `MethodError` naming them is the honest answer there.
    base = jlname === :pow ? :(Base.:^) :
           isdefined(Base, jlname) ? :(Base.$jlname) : nothing
    base === nothing || @eval @inline $disp(args...) = $base(args...)
end

_rewritable(f, nargs::Int) = f isa Symbol && nargs in get(_ARITIES, f, ())

# Positional, unsplatted arguments only: keywords and splats never reach a SIMDMath
# method, so such calls are left exactly as written.
_plain_args(args) = !any(a -> a isa Expr && a.head in (:parameters, :kw, :...), args)

_rewrite(x) = x
function _rewrite(ex::Expr)
    # Quoted code is data, not code that runs here.
    ex.head === :quote && return ex
    if ex.head === :call && _plain_args(ex.args[2:end])
        f, args = ex.args[1], ex.args[2:end]
        if f === :^ && length(args) == 2
            # `x^2` lowers to `Base.literal_pow`, i.e. `x*x`; a `pow` call would be
            # both slower and less accurate, so literal integer exponents stay put.
            args[2] isa Integer && return Expr(:call, :^, map(_rewrite, args)...)
            return Expr(:call, GlobalRef(SIMDMath, _dispatcher(:pow)), map(_rewrite, args)...)
        elseif _rewritable(f, length(args))
            return Expr(:call, GlobalRef(SIMDMath, _dispatcher(f)), map(_rewrite, args)...)
        end
    end
    # A method definition's signature is a name being bound, not a call being made.
    if ex.head in (:function, :->) ||
       (ex.head === :(=) && ex.args[1] isa Expr && ex.args[1].head in (:call, :where, :(::)))
        return Expr(ex.head, ex.args[1], map(_rewrite, ex.args[2:end])...)
    end
    return Expr(ex.head, map(_rewrite, ex.args)...)
end

"""
    SIMDMath.@simdmath expr

Rewrite the math calls inside `expr` to their `SIMDMath` equivalents, without
importing anything into the enclosing module:

```julia
using AppleAccelerate.SIMDMath: @simdmath

function weighted_logsum(X, W)
    u = zero(eltype(X))
    @simdmath @simd for i in eachindex(X, W)
        @inbounds u += W[i] * log(X[i])^W[i]
    end
    u
end
```

This is the alternative to `using AppleAccelerate.SIMDMath: log`, which replaces
`log` everywhere in the module rather than in one loop. `@simdmath` goes outside or
inside `@simd`; it does not add `@simd` or `@inbounds` for you, and the loop still
has to vectorise for any of this to matter (see [`SIMDMath`](@ref)).

What is rewritten:

  * unqualified calls to a name in the "Available functions" list of
    [`SIMDMath`](@ref), with a matching number of positional arguments;
  * `x^y`, which becomes `pow(x, y)` -- except when `y` is a literal integer
    (`x^2`), which Julia already lowers to multiplications.

What is left alone: qualified calls (`Base.log(x)`), broadcasts (`log.(x)`, `x .^ y`),
calls with keyword or splatted arguments, the signature of a method defined inside
`expr`, anything inside a quoted expression, and every function SIMDMath does not
provide.

A rewritten call only uses the SIMD routine when all its arguments are `Float32` or
all are `Float64`. For any other argument types it falls back to the `Base`
function, so index arithmetic such as `rem(i, 4)` and integer powers `x^n` keep
working unchanged. (`nextafter` and `remainder` have no `Base` counterpart to fall
back to.)

The rewrite is purely syntactic: a local variable or argument that happens to be
called `log` and is then *called* will be rewritten too.

The accuracy caveats of [`SIMDMath`](@ref) apply unchanged -- these routines are less
accurate than `Base`.
"""
macro simdmath(ex)
    esc(_rewrite(ex))
end

end # module SIMDMath
