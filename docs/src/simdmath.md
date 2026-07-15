# SIMD Math in `@simd` Loops

`AppleAccelerate.SIMDMath` exposes scalar math functions that LLVM replaces with
Apple's SIMD math routines when they appear inside a loop that vectorises.

This is the non-array counterpart to [Array Operations](array.md). Where
`AppleAccelerate.log!(out, X)` needs a whole array in and a whole array out,
`SIMDMath.log(x)` is an ordinary scalar call you can put in the middle of a loop.

## When to use it

!!! important
    **`SIMDMath` is not a faster replacement for the array API.** On an M-series
    machine, vForce's `vv*` routines (`AppleAccelerate.log!` and friends) are about
    **2x faster per element** than these SIMD routines at *every* array size measured
    (N = 16 through 500,000). There is no small-N crossover, and vForce wins even
    when it needs extra temporary arrays and extra passes over memory.

    If your data is a dense `Array` and you can call `AppleAccelerate.log!`, do that.

`SIMDMath` is for the cases where the array API does not apply:

  * strided or otherwise non-contiguous access;
  * values computed on the fly rather than read out of an array;
  * loops with control flow, or over iterables that are not dense `Array`s;
  * hot paths where you cannot allocate and cannot preallocate a temporary.

In those situations the realistic alternative is a scalar `Base` loop, and `SIMDMath`
is roughly **2–4x faster for `Float32`** and **1.2–2x for `Float64`**.

```julia
using AppleAccelerate.SIMDMath: log

# strided access: vForce would need a gather into a temporary first
function logsum_strided(X, stride)
    u = zero(eltype(X))
    @simd for i in 1:stride:length(X)
        @inbounds u += log(X[i])
    end
    u
end
```

## Accuracy

!!! warning
    These routines trade accuracy for speed and are **less accurate than `Base`**,
    whose math functions are correctly rounded (≤ 0.5 ULP). Worst case measured
    against `Base` is about **3 ULP**.

    Accuracy is also **not consistent across a single call site**: if the loop
    vectorises you get the SIMD routine, and if it does not — an unvectorisable
    loop, or the scalar remainder of one that did — you get the correctly-rounded
    libm call instead. Two runs over arrays of different length can therefore give
    slightly different answers. Do not use `SIMDMath` where bitwise reproducibility
    matters.

## Confirming you got the vectorised form

A call only becomes a SIMD call if the surrounding loop actually vectorises, which
in practice means `@simd` and usually `@inbounds`. Nothing warns you if it does not
— the code stays correct and simply runs at scalar speed. To check, look for the
symbol in the generated code:

```julia
using InteractiveUtils
@code_native logsum_strided(rand(1000), 3)   # expect a call to _simd_log_d2
```

`Float32` runs 4 lanes at a time (`_simd_*_f4`), `Float64` runs 2 (`_simd_*_d2`).

## Available functions

All are defined for `Float32` and `Float64` only.

| | |
|---|---|
| One argument | `acos`, `acosh`, `asin`, `asinh`, `atan`, `atanh`, `cbrt`, `cos`, `cosh`, `cospi`, `erf`, `erfc`, `exp`, `exp10`, `exp2`, `expm1`, `log`, `log10`, `log1p`, `log2`, `sin`, `sinh`, `sinpi`, `tan`, `tanh`, `tanpi`, `tgamma` |
| Two arguments | `atan(y, x)`, `hypot`, `nextafter`, `pow`, `rem`, `remainder` |

`erf`, `erfc`, `tgamma`, `nextafter`, `pow` and `remainder` have no `Base`
counterpart with matching semantics, so they are named after their C equivalents.
Use `pow(x, y)` rather than `x^y`.

Functions LLVM already lowers to native instructions — `sqrt`, `floor`, `ceil`,
`round`, `trunc`, `fma`, `abs` — are deliberately **absent**: routing those through
a library call would be slower than the instruction the compiler emits anyway.

`lgamma` is also absent, because it writes the global `signgam` and so cannot be
declared side-effect-free.

## How it works

LLVM's loop vectoriser can replace a scalar call inside a vectorising loop with a
call to an equivalent SIMD routine, if it is told which routine to use. `clang`
spells this `-fveclib=Accelerate`. The underlying mechanism is the LLVM *vector
function ABI* (VFABI), which is driven entirely by IR metadata, so Julia can reach
it through `llvmcall` — no compiler flag and no patched Julia required.

Each function is emitted as a scalar call to the real libm symbol (`logf`), carrying
a `"vector-function-abi-variant"` attribute naming the SIMD replacement:

```
_ZGV_LLVM_N4v_logf(_simd_log_f4)
```

which reads as "`_LLVM_` ISA, **N**ot masked, **4** lanes, one **v**ector operand,
mapping scalar `logf` onto vector `_simd_log_f4`". The vector routine must also be
declared in the module *and* anchored in `@llvm.compiler.used`, or it gets stripped
before the vectoriser runs and the mapping is silently dropped.

The routines come from libsystem_m's `<simd/math.h>` rather than Accelerate's own
`vfp.h`. Both are public Apple API, but `vfp.h` is `Float32`-only and covers fewer
functions, whereas `<simd/math.h>` covers both widths — and since `Float64` is
Julia's default, using one source for both keeps accuracy consistent across types.

```@docs
AppleAccelerate.SIMDMath
```
