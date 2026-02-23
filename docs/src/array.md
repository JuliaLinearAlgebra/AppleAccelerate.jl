# Array Operations

AppleAccelerate wraps Apple's vecLib (`vv*` and `vDSP_*` functions) to provide accelerated element-wise operations on `Array{Float32}` and `Array{Float64}`.

These functions are **not exported** to avoid conflicts with Base. Access them via the `AppleAccelerate.` prefix or use `@replaceBase` to override Base functions.

## Element-wise Math Functions (vecLib)

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

These functions wrap `vDSP_*` routines. Each has an allocating variant `f(X)` and a mutating variant `f!(result, X)`:

| Function | Description |
|----------|-------------|
| `vneg(X)` | Negate: `-X` |
| `vnabs(X)` | Negative absolute value: `-abs.(X)` |
| `vabs(X)` | Absolute value (vDSP) |
| `vsq(X)` | Square: `X.^2` |
| `vssq(X)` | Signed square: `X .* abs.(X)` |
| `vfrac(X)` | Fractional part: `X .- trunc.(X)` |
| `vreverse(X)` | Reverse vector (`vreverse!` operates in-place) |

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
| `dot(X, Y)` | Dot product: `sum(X .* Y)` |
| `distancesq(X, Y)` | Squared distance: `sum((X .- Y).^2)` |

## Vector-Vector Operations

### Arithmetic

| Function | Description |
|----------|-------------|
| `vadd(X, Y)` | Element-wise addition |
| `vsub(X, Y)` | Element-wise subtraction |
| `vmul(X, Y)` | Element-wise multiplication |
| `vdiv(X, Y)` | Element-wise division |

### Comparison & Distance

| Function | Description |
|----------|-------------|
| `vmax(X, Y)` | Element-wise maximum |
| `vmin(X, Y)` | Element-wise minimum |
| `vmaxmg(X, Y)` | Element-wise max of magnitudes: `max.(abs.(X), abs.(Y))` |
| `vminmg(X, Y)` | Element-wise min of magnitudes: `min.(abs.(X), abs.(Y))` |
| `vdist(X, Y)` | Element-wise distance: `hypot.(X, Y)` |
| `vtmerg(X, Y)` | Tapered merge of two vectors |

## Vector-Scalar Operations

| Function | Description |
|----------|-------------|
| `vsadd(X, c)` | `X .+ c` |
| `vssub(X, c)` | `X .- c` |
| `svsub(X, c)` | `c .- X` |
| `vsmul(X, c)` | `X .* c` |
| `vsdiv(X, c)` | `X ./ c` |
| `svdiv(X, c)` | `c ./ X` (scalar divided by vector) |

All vector operations have mutating `!` variants (e.g., `vadd!(result, X, Y)`).

## Compound Arithmetic

These operations fuse multiple arithmetic steps into a single vDSP call for better performance.

### Three-vector operations

| Function | Description |
|----------|-------------|
| `vam(A, B, C)` | `(A .+ B) .* C` |
| `vsbm(A, B, C)` | `(A .- B) .* C` |
| `venvlp(A, B, C)` | Signal envelope |

### Four-vector operations

| Function | Description |
|----------|-------------|
| `vaam(A, B, C, D)` | `(A .+ B) .* (C .+ D)` |
| `vsbsbm(A, B, C, D)` | `(A .- B) .* (C .- D)` |
| `vasbm(A, B, C, D)` | `(A .+ B) .* (C .- D)` |
| `vpythg(A, B, C, D)` | `sqrt.((A .- C).^2 .+ (B .- D).^2)` |

### Vector-vector-scalar operations

| Function | Description |
|----------|-------------|
| `vasm(A, B, c)` | `(A .+ B) .* c` |
| `vsbsm(A, B, c)` | `(A .- B) .* c` |
| `vsma(A, b, C)` | `A .* b .+ C` |
| `vsmsa(A, b, c)` | `A .* b .+ c` |

### Dual output

| Function | Description |
|----------|-------------|
| `vaddsub(A, B)` | Returns `(A .+ B, B .- A)` as two vectors |

## Clipping & Thresholding

| Function | Description |
|----------|-------------|
| `vclip(X, low, high)` | Clip values to `[low, high]` |
| `viclip(X, low, high)` | Inverted clip: pass values outside `[low, high]` |
| `vthr(X, threshold)` | `max.(X, threshold)` |
| `vthres(X, threshold)` | `X[n] >= threshold ? X[n] : 0` |
| `vcmprs(X, gate)` | Compress: keep `X[n]` where `gate[n] != 0` |

## Type Conversion

| Function | Description |
|----------|-------------|
| `vdouble(X::Vector{Float32})` | Convert to `Vector{Float64}` |
| `vsingle(X::Vector{Float64})` | Convert to `Vector{Float32}` |

## Ramp Generation

| Function | Description |
|----------|-------------|
| `vramp(start, step, n)` | Generate ramp: `start + i*step` for `i = 0, ..., n-1` |
| `vrampmul(X, start, step)` | Multiply vector by generated ramp |

## Integration & Running Operations

| Function | Description |
|----------|-------------|
| `vrsum(X, scale)` | Running sum scaled by `scale` |
| `vsimps(X, step)` | Simpson's rule integration |
| `vtrapz(X, step)` | Trapezoidal integration |
| `vswsum(X, window)` | Sliding window sum |
| `vswmax(X, window)` | Sliding window maximum |

## Interpolation

| Function | Description |
|----------|-------------|
| `vintb(A, B, t)` | Linear interpolation: `A .+ t .* (B .- A)` |
| `vlint(table, indices)` | Linear interpolation lookup from table |
| `vqint(table, indices)` | Quadratic interpolation lookup from table |

## Polynomial Evaluation

| Function | Description |
|----------|-------------|
| `vpoly(coeffs, X)` | Evaluate polynomial at each point in `X`. Coefficients are ordered highest degree first: `[a_P, a_{P-1}, ..., a_1, a_0]` |

## Normalization

| Function | Description |
|----------|-------------|
| `vnormalize(X)` | Returns `(normalized_X, mean, stddev)` where `normalized_X = (X .- mean) ./ stddev` |

## Zero Crossings

| Function | Description |
|----------|-------------|
| `nzcros(X, max_crossings=0)` | Find zero crossings. Returns `(indices, count)`. Default `max_crossings=0` means up to `length(X)`. |

## Decibel Conversion

| Function | Description |
|----------|-------------|
| `vdbcon(X, ref, power=true)` | Convert to decibels. `power=true`: `10*log10(X/ref)`. `power=false`: `20*log10(X/ref)`. |

## Broadcasting

AppleAccelerate overrides `Base.copy` and `Base.copyto!` for `Broadcasted` objects, so that broadcasting syntax like `f.(X)` automatically uses the accelerated implementation.

## Replacing Base Functions

```@docs
AppleAccelerate.@replaceBase
```
