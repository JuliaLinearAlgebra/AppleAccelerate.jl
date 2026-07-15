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
| One argument | `acos`, `acosh`, `asin`, `asinh`, `atan`, `atanh`, `cbrt`, `cos`, `cosh`, `cospi`, `exp`, `exp10`, `exp2`, `expm1`, `log`, `log10`, `log1p`, `log2`, `sin`, `sinh`, `sinpi`, `tan`, `tanh`, `tanpi` |
| Two arguments | `atan(y, x)`, `hypot`, `nextafter`, `pow`, `rem`, `remainder` |

`nextafter`, `pow` and `remainder` have no `Base` counterpart with matching
semantics, so they are named after their C equivalents. Use `pow(x, y)` rather than
`x^y`.

Deliberately **absent**:

  * `sqrt`, `floor`, `ceil`, `round`, `trunc`, `fma`, `abs` — LLVM already lowers
    these to native instructions, so a library call would be *slower*.
  * `erf`, `erfc`, `tgamma` — `<simd/math.h>` declares `_simd_erf_f4` and friends,
    but they measure at 0.97–0.99x of a plain scalar libm loop: they are scalar
    loops behind a vector signature and buy nothing.
  * `lgamma` — writes the global `signgam`, so it cannot be declared
    side-effect-free.

## Measured speedups

`out[i] = f(x[i])` over N = 100,000 on an M-series machine. **`SM/Base`** is
SIMDMath against a scalar `Base` loop (higher is better); **`SM/vForce`** is against
`AppleAccelerate.f!` (below 1.0 means vForce wins).

| func | Float32 SM/Base | Float32 SM/vForce | Float64 SM/Base | Float64 SM/vForce |
|---|---|---|---|---|
| `sinpi` / `cospi` | **26x** | 0.61x | **7.5x** | 0.45x |
| `tanpi` | **22x** | 0.78x | **5.9x** | 0.57x |
| `atan(y,x)` | **10.7x** | 0.73x | 3.9x | 0.62x |
| `atan` | **10.2x** | 0.69x | 3.4x | 0.53x |
| `asin` / `acos` | **8.6x** | 0.65x | 3.1–3.7x | 0.51–0.63x |
| `sin` / `cos` | **8.2x** | 0.57x | 2.5x | 0.51x |
| `rem` | 6.9x | 0.97x | 3.1x | 0.99x |
| `expm1` | 6.7x | 0.71x | 3.4x | 0.58x |
| `pow` | 6.4x | 0.78x | 2.6x | 0.66x |
| `tan` | 6.4x | 0.64x | 4.0x | 0.52x |
| `tanh` | 5.7x | 0.78x | 3.2x | 0.57x |
| `atanh` | 5.2x | 0.57x | 2.1x | 0.54x |
| `exp2` / `exp` / `exp10` | 4.2–5.0x | 0.57x | 1.5x | 0.68x |
| `log` / `log2` / `log10` / `log1p` | 3.4–4.5x | 0.55–0.67x | 1.8–2.2x | 0.59–0.65x |
| `sinh` | 4.4x | 0.66x | 1.6x | 0.57x |
| `acosh` | 4.1x | 0.56x | 1.1x | 0.52x |
| `asinh` | 2.5x | 0.65x | 1.1x | 0.60x |
| `cbrt` | 2.1x | 0.73x | 1.3x | 0.66x |
| `hypot` | 1.2x | — | 4.5x | — |
| `cosh` | 1.1x | 0.63x | 1.5x | 0.55x |
| `nextafter` | — | 0.18x | — | 0.33x |

Reading this:

  * **`SM/vForce` is below 1.0 in every single row.** vForce wins across the board,
    which is why the guidance above is to prefer it whenever it applies. `nextafter`
    is the extreme case at 0.18x — always prefer `AppleAccelerate.nextafter!`.
  * The `Float32` trigonometric functions are where this shines against `Base`,
    especially `sinpi`/`cospi`/`tanpi`, where `Base` is unusually slow.
  * `cosh` and `hypot` barely beat `Base` in `Float32` because `Base` implements
    them natively there — those loops vectorise with no library call at all, so
    there is little left to win.
  * `hypot` and `exp10` have no vForce equivalent, so `SIMDMath` is the only
    accelerated option for them.

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
