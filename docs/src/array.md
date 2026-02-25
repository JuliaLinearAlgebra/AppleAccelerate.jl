# Array Operations

AppleAccelerate wraps Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) (`vv*`) and [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) (`vDSP_*`) functions to provide accelerated element-wise operations on `Array{Float32}` and `Array{Float64}`.

These functions are **not exported** to avoid conflicts with Base. Access them via the `AppleAccelerate.` prefix or use `@replaceBase` to override Base functions.

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

```@docs
AppleAccelerate.vneg
AppleAccelerate.vnabs
AppleAccelerate.vabs
AppleAccelerate.vsq
AppleAccelerate.vssq
AppleAccelerate.vfrac
AppleAccelerate.vreverse!
AppleAccelerate.vreverse
```

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

```@docs
AppleAccelerate.dot
AppleAccelerate.distancesq
```

## Vector-Vector Arithmetic

```@docs
AppleAccelerate.vadd
AppleAccelerate.vadd!
AppleAccelerate.vsub
AppleAccelerate.vsub!
AppleAccelerate.vmul
AppleAccelerate.vmul!
AppleAccelerate.vdiv
AppleAccelerate.vdiv!
```

## Two-Vector Comparison & Distance

```@docs
AppleAccelerate.vmax
AppleAccelerate.vmin
AppleAccelerate.vmaxmg
AppleAccelerate.vminmg
AppleAccelerate.vdist
AppleAccelerate.vtmerg
```

## Vector-Scalar Operations

```@docs
AppleAccelerate.vsadd
AppleAccelerate.vsadd!
AppleAccelerate.vssub
AppleAccelerate.vssub!
AppleAccelerate.svsub
AppleAccelerate.svsub!
AppleAccelerate.vsmul
AppleAccelerate.vsmul!
AppleAccelerate.vsdiv
AppleAccelerate.vsdiv!
AppleAccelerate.svdiv
```

## Compound Arithmetic

These operations fuse multiple arithmetic steps into a single vDSP call for better performance.

### Three-vector operations

```@docs
AppleAccelerate.vam
AppleAccelerate.vsbm
AppleAccelerate.venvlp
```

### Four-vector operations

```@docs
AppleAccelerate.vaam
AppleAccelerate.vsbsbm
AppleAccelerate.vasbm
AppleAccelerate.vpythg
```

### Vector-vector-scalar operations

```@docs
AppleAccelerate.vasm
AppleAccelerate.vsbsm
AppleAccelerate.vsma
AppleAccelerate.vsmsa
```

### Dual output

```@docs
AppleAccelerate.vaddsub
```

## Clipping & Thresholding

```@docs
AppleAccelerate.vclip
AppleAccelerate.viclip
AppleAccelerate.vthr
AppleAccelerate.vthres
AppleAccelerate.vcmprs
```

## Type Conversion

```@docs
AppleAccelerate.vdouble
AppleAccelerate.vsingle
```

## Ramp Generation

```@docs
AppleAccelerate.vramp
AppleAccelerate.vrampmul
```

## Integration & Running Operations

```@docs
AppleAccelerate.vrsum
AppleAccelerate.vsimps
AppleAccelerate.vtrapz
AppleAccelerate.vswsum
AppleAccelerate.vswmax
```

## Interpolation

```@docs
AppleAccelerate.vintb
AppleAccelerate.vlint
AppleAccelerate.vqint
```

## Polynomial Evaluation

```@docs
AppleAccelerate.vpoly
```

## Normalization

```@docs
AppleAccelerate.vnormalize
```

## Zero Crossings

```@docs
AppleAccelerate.nzcros
```

## Decibel Conversion

```@docs
AppleAccelerate.vdbcon
```

## Complex Array Operations

AppleAccelerate also wraps vDSP's split-complex functions for `Vector{Complex{Float32}}` and `Vector{Complex{Float64}}`. These extend existing function names (e.g., `vneg`, `vabs`, `vmul`) with methods that dispatch on complex element types — no naming conflicts with the real-valued versions above.

### Complex unary operations

| Function | Description | vDSP function |
|----------|-------------|---------------|
| `vneg(X)` / `vneg!(result, X)` | Negate: `-X` | `vDSP_zvneg` |
| `vabs(X)` / `vabs!(result, X)` | Modulus: `abs.(X)` | `vDSP_zvabs` |

```@docs
AppleAccelerate.vconj
AppleAccelerate.vcopy
```

### Complex → real operations

```@docs
AppleAccelerate.vphase
AppleAccelerate.vmags
AppleAccelerate.vmagsa
```

### Complex binary operations

| Function | Description | vDSP function |
|----------|-------------|---------------|
| `vmul(X, Y)` / `vmul!(result, X, Y)` | Element-wise multiply: `X .* Y` | `vDSP_zvmul` |
| `vdiv(X, Y)` / `vdiv!(result, X, Y)` | Element-wise divide: `X ./ Y` | `vDSP_zvdiv` |
| `vsmul(X, c)` / `vsmul!(result, X, c)` | Scalar multiply: `X .* c` (complex scalar) | `vDSP_zvzsml` |
| `dot(X, Y)` | Unconjugated dot product: `sum(X .* Y)` | `vDSP_zdotpr` |

### Coordinate conversion

```@docs
AppleAccelerate.polar
AppleAccelerate.rect
```

## Broadcasting

AppleAccelerate overrides `Base.copy` and `Base.copyto!` for `Broadcasted` objects, so that broadcasting syntax like `f.(X)` automatically uses the accelerated implementation.

## Replacing Base Functions

```@docs
AppleAccelerate.@replaceBase
```

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
