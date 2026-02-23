# Array Operations

AppleAccelerate wraps Apple's vecLib (`vv*` and `vDSP_*` functions) to provide accelerated element-wise operations on `Array{Float32}` and `Array{Float64}`.

These functions are **not exported** to avoid conflicts with Base. Access them via the `AppleAccelerate.` prefix or use [`@replaceBase`](@ref) to override Base functions.

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

## Vector-Vector Operations

| Function | Description |
|----------|-------------|
| `vadd(X, Y)` | Element-wise addition |
| `vsub(X, Y)` | Element-wise subtraction |
| `vmul(X, Y)` | Element-wise multiplication |
| `vdiv(X, Y)` | Element-wise division |

## Vector-Scalar Operations

| Function | Description |
|----------|-------------|
| `vsadd(X, c)` | `X .+ c` |
| `vssub(X, c)` | `X .- c` |
| `svsub(X, c)` | `c .- X` |
| `vsmul(X, c)` | `X .* c` |
| `vsdiv(X, c)` | `X ./ c` |

All vector operations have mutating `!` variants (e.g., `vadd!(result, X, Y)`).

## Broadcasting

AppleAccelerate overrides `Base.copy` and `Base.copyto!` for `Broadcasted` objects, so that broadcasting syntax like `f.(X)` automatically uses the accelerated implementation.

## Replacing Base Functions

```@docs
AppleAccelerate.@replaceBase
```
