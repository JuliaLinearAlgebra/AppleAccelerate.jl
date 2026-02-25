# Array Operations

AppleAccelerate wraps Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) (`vv*`) and [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) (`vDSP_*`) functions to provide accelerated element-wise operations on `Array{Float32}` and `Array{Float64}`.

These functions are **not exported** to avoid conflicts with Base. Access them via the `AppleAccelerate.` prefix.

## Element-wise Math Functions

### One-argument functions

Each function `f` has an allocating variant `f(X)` and a mutating variant `f!(out, X)`:

| Function | Description |
|----------|-------------|
| `ceil`, `floor`, `trunc`, `round` | Rounding |
| `sqrt`, `rsqrt`, `rec` | Square root, reciprocal square root, reciprocal |
| `exp`, `exp2`, `expm1` | Exponentials |
| `log`, `log1p`, `log2`, `log10` | Logarithms |
| `sin`, `sinpi`, `cos`, `cospi`, `tan`, `tanpi` | Trigonometric |
| `asin`, `acos`, `atan` | Inverse trigonometric |
| `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` | Hyperbolic |
| `abs`, `exponent` | Miscellaneous |

### Two-argument functions

| Function | Description |
|----------|-------------|
| `copysign(X, Y)` | Copy sign of Y to X |
| `rem(X, Y)` | Element-wise remainder |
| `div_float(X, Y)` | Element-wise division (via vecLib) |
| `atan(X, Y)` | Two-argument arctangent |
| `pow(X, Y)` | Element-wise power |

### Special return types

| Function | Description |
|----------|-------------|
| `sincos(X)` | Returns `(sin(X), cos(X))` tuple |
| `cis(X)` | Returns `Complex` array `cos(X) + im*sin(X)` |

## Unary vDSP Operations

| Function | Description |
|----------|-------------|
| [`vneg`](@ref AppleAccelerate.vneg) | Negate each element: `result[i] = -X[i]` |
| [`vnabs`](@ref AppleAccelerate.vnabs) | Negative absolute value: `result[i] = -|X[i]|` |
| [`vabs`](@ref AppleAccelerate.vabs) | Absolute value: `result[i] = \|X[i]\|` |
| [`vsq`](@ref AppleAccelerate.vsq) | Square each element: `result[i] = X[i]^2` |
| [`vssq`](@ref AppleAccelerate.vssq) | Signed square: `result[i] = X[i] * \|X[i]\|` |
| [`vfrac`](@ref AppleAccelerate.vfrac) | Fractional part: `result[i] = X[i] - trunc(X[i])` |
| [`vreverse!`](@ref AppleAccelerate.vreverse!) | Reverse vector in-place |
| [`vreverse`](@ref AppleAccelerate.vreverse) | Return a reversed copy |

## Vector Reductions

| Function | Description |
|----------|-------------|
| `maximum(X)`, `minimum(X)` | Max/min value |
| `findmax(X)`, `findmin(X)` | Max/min value and index |
| `sum(X)`, `mean(X)` | Sum and mean |
| `meanmag(X)` | Mean of absolute values |
| `meansqr(X)` | Mean of squares |
| `meanssqr(X)` | Mean of signed squares |
| `summag(X)` | Sum of absolute values |
| `sumsqr(X)` | Sum of squares |
| `sumssqr(X)` | Sum of signed squares |
| [`dot`](@ref AppleAccelerate.dot) | Dot product: `sum(X .* Y)` |
| [`distancesq`](@ref AppleAccelerate.distancesq) | Squared Euclidean distance: `sum((X .- Y).^2)` |

## Vector-Vector Arithmetic

| Function | Description |
|----------|-------------|
| [`vadd`](@ref AppleAccelerate.vadd) / [`vadd!`](@ref AppleAccelerate.vadd!) | Element-wise addition |
| [`vsub`](@ref AppleAccelerate.vsub) / [`vsub!`](@ref AppleAccelerate.vsub!) | Element-wise subtraction |
| [`vmul`](@ref AppleAccelerate.vmul) / [`vmul!`](@ref AppleAccelerate.vmul!) | Element-wise multiplication |
| [`vdiv`](@ref AppleAccelerate.vdiv) / [`vdiv!`](@ref AppleAccelerate.vdiv!) | Element-wise division |

## Two-Vector Comparison & Distance

| Function | Description |
|----------|-------------|
| [`vmax`](@ref AppleAccelerate.vmax) | Element-wise maximum |
| [`vmin`](@ref AppleAccelerate.vmin) | Element-wise minimum |
| [`vmaxmg`](@ref AppleAccelerate.vmaxmg) | Element-wise maximum magnitude |
| [`vminmg`](@ref AppleAccelerate.vminmg) | Element-wise minimum magnitude |
| [`vdist`](@ref AppleAccelerate.vdist) | Element-wise Euclidean distance |
| [`vtmerg`](@ref AppleAccelerate.vtmerg) | Tapered merge of two vectors |

## Vector-Scalar Operations

| Function | Description |
|----------|-------------|
| [`vsadd`](@ref AppleAccelerate.vsadd) / [`vsadd!`](@ref AppleAccelerate.vsadd!) | Vector + scalar |
| [`vssub`](@ref AppleAccelerate.vssub) / [`vssub!`](@ref AppleAccelerate.vssub!) | Vector - scalar |
| [`svsub`](@ref AppleAccelerate.svsub) / [`svsub!`](@ref AppleAccelerate.svsub!) | Scalar - vector |
| [`vsmul`](@ref AppleAccelerate.vsmul) / [`vsmul!`](@ref AppleAccelerate.vsmul!) | Vector * scalar |
| [`vsdiv`](@ref AppleAccelerate.vsdiv) / [`vsdiv!`](@ref AppleAccelerate.vsdiv!) | Vector / scalar |
| [`svdiv`](@ref AppleAccelerate.svdiv) | Scalar / vector |

## Compound Arithmetic

These operations fuse multiple arithmetic steps into a single vDSP call for better performance.

### Three-vector operations

| Function | Description |
|----------|-------------|
| [`vam`](@ref AppleAccelerate.vam) | `(A + B) * C` |
| [`vsbm`](@ref AppleAccelerate.vsbm) | `(A - B) * C` |
| [`venvlp`](@ref AppleAccelerate.venvlp) | Signal envelope |

### Four-vector operations

| Function | Description |
|----------|-------------|
| [`vaam`](@ref AppleAccelerate.vaam) | `(A + B) * (C + D)` |
| [`vsbsbm`](@ref AppleAccelerate.vsbsbm) | `(A - B) * (C - D)` |
| [`vasbm`](@ref AppleAccelerate.vasbm) | `(A + B) * (C - D)` |
| [`vpythg`](@ref AppleAccelerate.vpythg) | Pythagorean distance |

### Vector-vector-scalar operations

| Function | Description |
|----------|-------------|
| [`vasm`](@ref AppleAccelerate.vasm) | `(A + B) * c` |
| [`vsbsm`](@ref AppleAccelerate.vsbsm) | `(A - B) * c` |
| [`vsma`](@ref AppleAccelerate.vsma) | `A * b + C` |
| [`vsmsa`](@ref AppleAccelerate.vsmsa) | `A * b + c` |

### Dual output

| Function | Description |
|----------|-------------|
| [`vaddsub`](@ref AppleAccelerate.vaddsub) | Simultaneous add and subtract: returns `(A .+ B, A .- B)` |

## Clipping & Thresholding

| Function | Description |
|----------|-------------|
| [`vclip`](@ref AppleAccelerate.vclip) | Clip values to `[low, high]` |
| [`viclip`](@ref AppleAccelerate.viclip) | Inverted clip: pass values outside `[low, high]` |
| [`vthr`](@ref AppleAccelerate.vthr) | Threshold: keep or clamp to threshold |
| [`vthres`](@ref AppleAccelerate.vthres) | Threshold to zero |
| [`vcmprs`](@ref AppleAccelerate.vcmprs) | Compress: gather elements where gate is nonzero |

## Type Conversion

| Function | Description |
|----------|-------------|
| [`vdouble`](@ref AppleAccelerate.vdouble) | Convert Float32 to Float64 |
| [`vsingle`](@ref AppleAccelerate.vsingle) | Convert Float64 to Float32 |

## Ramp Generation

| Function | Description |
|----------|-------------|
| [`vramp`](@ref AppleAccelerate.vramp) | Generate a ramp: `start + i * step` |
| [`vrampmul`](@ref AppleAccelerate.vrampmul) | Multiply vector by a generated ramp |
| [`vrampmul2`](@ref AppleAccelerate.vrampmul2) | Stereo ramp multiply (two outputs) |

## Linear Average

| Function | Description |
|----------|-------------|
| [`vavlin`](@ref AppleAccelerate.vavlin) | Weighted linear average of two vectors |

## Integration & Running Operations

| Function | Description |
|----------|-------------|
| [`vrsum`](@ref AppleAccelerate.vrsum) | Running sum scaled by `scale` |
| [`vsimps`](@ref AppleAccelerate.vsimps) | Simpson's rule integration |
| [`vtrapz`](@ref AppleAccelerate.vtrapz) | Trapezoidal integration |
| [`vswsum`](@ref AppleAccelerate.vswsum) | Sliding window sum |
| [`vswmax`](@ref AppleAccelerate.vswmax) | Sliding window maximum |

## Interpolation

| Function | Description |
|----------|-------------|
| [`vintb`](@ref AppleAccelerate.vintb) | Linear interpolation: `A + t * (B - A)` |
| [`vlint`](@ref AppleAccelerate.vlint) | Linear interpolation from lookup table |
| [`vqint`](@ref AppleAccelerate.vqint) | Quadratic interpolation from lookup table |

## Polynomial Evaluation

| Function | Description |
|----------|-------------|
| [`vpoly`](@ref AppleAccelerate.vpoly) | Evaluate polynomial at each point |

## Normalization

| Function | Description |
|----------|-------------|
| [`vnormalize`](@ref AppleAccelerate.vnormalize) | Normalize to zero mean and unit standard deviation |

## Zero Crossings

| Function | Description |
|----------|-------------|
| [`nzcros`](@ref AppleAccelerate.nzcros) | Find zero crossings |

## Decibel Conversion

| Function | Description |
|----------|-------------|
| [`vdbcon`](@ref AppleAccelerate.vdbcon) | Convert to decibels relative to a reference |

## Complex Array Operations

AppleAccelerate also wraps vDSP's split-complex functions for `Vector{Complex{Float32}}` and `Vector{Complex{Float64}}`. These extend existing function names (e.g., `vneg`, `vabs`, `vmul`) with methods that dispatch on complex element types — no naming conflicts with the real-valued versions above.

### Complex unary operations

| Function | Description |
|----------|-------------|
| `vneg(X)` / `vneg!(result, X)` | Negate: `-X` |
| `vabs(X)` / `vabs!(result, X)` | Modulus: `abs.(X)` |
| [`vconj`](@ref AppleAccelerate.vconj) | Complex conjugate |
| [`vcopy`](@ref AppleAccelerate.vcopy) | Copy via split-complex move |

### Complex → real operations

| Function | Description |
|----------|-------------|
| [`vphase`](@ref AppleAccelerate.vphase) | Complex phase (angle) |
| [`vmags`](@ref AppleAccelerate.vmags) | Squared magnitude (`abs2`) |
| [`vmagsa`](@ref AppleAccelerate.vmagsa) | Squared magnitude + accumulate |

### Complex binary operations

| Function | Description |
|----------|-------------|
| `vmul(X, Y)` / `vmul!(result, X, Y)` | Element-wise multiply: `X .* Y` |
| `vdiv(X, Y)` / `vdiv!(result, X, Y)` | Element-wise divide: `X ./ Y` |
| `vsmul(X, c)` / `vsmul!(result, X, c)` | Scalar multiply (complex scalar) |
| `dot(X, Y)` | Unconjugated dot product: `sum(X .* Y)` |

### Coordinate conversion

| Function | Description |
|----------|-------------|
| [`polar`](@ref AppleAccelerate.polar) | Cartesian to polar coordinates |
| [`rect`](@ref AppleAccelerate.rect) | Polar to Cartesian coordinates |

## Broadcasting

AppleAccelerate overrides `Base.copy` and `Base.copyto!` for `Broadcasted` objects, so that broadcasting syntax like `f.(X)` automatically uses the accelerated implementation.

## Benchmarks

Performance comparison of vDSP array operations vs Julia Base equivalents (`map(Base.f, X)` for unary, `@simd` loops for binary/compound). Run with `julia --project=test/bench test/bench/run_benchmarks.jl array` ([source](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl)).

### Unary Math Functions

Transcendental functions show the biggest gains — vDSP is 7–19× faster for `sin`/`cos` on Float32:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| exp | Float64 | 100,000 | 171 | 440 | 2.6× |
| log | Float64 | 100,000 | 191 | 652 | 3.4× |
| sin | Float64 | 100,000 | 152 | 1,051 | 6.9× |
| cos | Float64 | 100,000 | 160 | 1,083 | 6.8× |
| sqrt | Float64 | 100,000 | 42 | 97 | 2.3× |
| exp | Float32 | 100,000 | 60 | 464 | 7.8× |
| log | Float32 | 100,000 | 80 | 542 | 6.8× |
| sin | Float32 | 100,000 | 54 | 1,012 | 18.8× |
| cos | Float32 | 100,000 | 55 | 1,035 | 18.9× |
| sqrt | Float32 | 100,000 | 40 | 98 | 2.5× |

### Reductions

`sum`/`maximum`/`minimum` are 1.1–2× faster. Note: `dot` is slower via vDSP because Julia's `LinearAlgebra.dot` already uses the Accelerate-forwarded BLAS:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| sum | Float64 | 1,000,000 | 128 | 143 | 1.1× |
| maximum | Float64 | 1,000,000 | 128 | 136 | 1.1× |
| minimum | Float64 | 1,000,000 | 127 | 134 | 1.1× |
| sum | Float32 | 1,000,000 | 62 | 77 | 1.2× |
| maximum | Float32 | 1,000,000 | 62 | 71 | 1.1× |
| minimum | Float32 | 1,000,000 | 63 | 71 | 1.1× |

### Binary Element-wise Ops

Addition and multiplication are memory-bandwidth-bound for Float64; vDSP is faster for Float32 at large sizes:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| vadd | Float64 | 1,000,000 | 354 | 360 | 1.0× |
| vmul | Float64 | 1,000,000 | 343 | 340 | 1.0× |
| vadd | Float32 | 1,000,000 | 133 | 167 | 1.3× |
| vmul | Float32 | 1,000,000 | 87 | 167 | 1.9× |

!!! note "Benchmark environment"
    Apple M2 Max, macOS 26, single-threaded. Julia 1.12.5, AppleAccelerate v0.6.0, LinearAlgebra v1.12.0. Times are minimum of 5 trials. Julia reference uses `map(Base.f, X)` for unary ops and `@inbounds @simd` loops for binary/compound ops. Run [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl) to reproduce.
