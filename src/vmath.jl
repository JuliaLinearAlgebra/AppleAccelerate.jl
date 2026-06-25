## vmath.jl — coverage notes for the vecLib math subframeworks ##
#
# This file documents the idiomatic-API status of the vecLib math subframeworks
# (vForce, vfp, vectorOps, vBasicOps, vBigNum). It deliberately contains no new
# wrappers: after auditing the SDK headers against the generated `LibAccelerate`
# raw layer, every function that has a *natural* Julia array surface is already
# wrapped (vForce), and everything else operates on SIMD/NEON register or
# big-integer vector types that have no clean idiomatic surface and is therefore
# intentionally left to the raw `LibAccelerate` layer.
#
# The notes live here (rather than only in the PR) so the rationale is
# discoverable from the source tree and from `?AppleAccelerate.VMATH_COVERAGE`.

"""
    AppleAccelerate.VMATH_COVERAGE

Documentation of idiomatic-wrapper coverage for the vecLib math subframeworks.

## vForce — fully covered (in `array.jl`)
The `vv*` transcendental / elementwise array functions all take plain
`const T *` pointers with a length, which maps directly onto Julia `Array{T}`
for `T <: Union{Float32,Float64}`. Every `vForce.h` entry has an idiomatic
allocating + mutating (`!`) wrapper in `array.jl`: `ceil`, `floor`, `trunc`,
`round`, `sqrt`, `rsqrt`, `rec`, `exp`, `exp2`, `expm1`, `log`, `log1p`,
`log2`, `log10`, `logb` (`exponent`), `sin`, `cos`, `tan`, `sinpi`, `cospi`,
`tanpi`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `asinh`,
`acosh`, `atanh`, `abs` (`fabs`), `copysign`, `rem` (`fmod`), `div_float`,
`pow`, `pows`, `cbrt`, `nextafter`, `remainder`, `sincos`, and `cis`
(`cosisin`). There are therefore no remaining vForce gaps to fill.

See [vForce documentation](https://developer.apple.com/documentation/accelerate/veclib).

## vfp — intentionally left to the raw layer
Every `vfp.h` function (`vsinf`, `vexpf`, `vsqrtf`, `vpowf`, `vremquof`,
`vscalbf`, `vsignbitf`, `vclassifyf`, `vtablelookup`, …) operates on the
128-bit SIMD type `vFloat` (4×`Float32` in a NEON/SSE register), taking and
returning *register* values rather than memory arrays. There is no idiomatic
Julia array surface for register-width SIMD operands, so these are accessible
only through the raw `LibAccelerate` layer.

## vectorOps — intentionally left to the raw layer
The `vS*` / `vIs*` BLAS-style routines (`vSgemm`, `vSdot`, `vSaxpy`,
`vIsamax`, `vSnorm2`, …) take `const vFloat *` operands and express their
length as a count of `vFloat` SIMD vectors (4 floats each), not as a plain
element count. That awkward, alignment-sensitive surface duplicates standard
BLAS, which Julia programs reach far more cleanly via `LinearAlgebra.BLAS`
(Accelerate is already wired in as the BLAS backend). Left to the raw layer.

## vBasicOps — intentionally left to the raw layer
`vBasicOps.h` provides integer arithmetic, shifts, and rotates on SIMD vector
types (`vUInt8`, `vSInt16`, … through the 128-bit `vU128`/`vS128`), again
operating on register-width vectors with no idiomatic array mapping. Raw layer
only.

## vBigNum — intentionally left to the raw layer
`vBigNum.h` implements fixed-width big-integer arithmetic over 256/512/1024-bit
SIMD struct types (`vU256`, `vU512`, `vU1024`, …). These have no natural Julia
counterpart (Julia uses `BigInt` for arbitrary precision), so they are left to
the raw layer.
"""
const VMATH_COVERAGE = nothing
