# Array Operations

AppleAccelerate wraps Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) (`vv*`) and [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) (`vDSP_*`) functions to provide accelerated element-wise operations on `Array{Float32}` and `Array{Float64}`.

These functions are **not exported** to avoid conflicts with Base. Access them via the `AppleAccelerate.` prefix.

```@setup array
using AppleAccelerate
```

## Element-wise Math Functions

These functions wrap Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) `vv*` routines.

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

```@example array
X = randn(Float64, 1000)

# Element-wise math — 3–19× faster than Base
Y_exp = AppleAccelerate.exp(X)
Y_sin = AppleAccelerate.sin(X)
Y_log = AppleAccelerate.log(X .+ 10)  # shift to positive domain

# Mutating variant (pre-allocate output)
out = similar(X)
AppleAccelerate.exp!(out, X)

# Broadcasting works automatically
Y_broadcast = AppleAccelerate.sin.(X)
nothing # hide
```

```@docs
AppleAccelerate.sincos
AppleAccelerate.cis
```

## Unary vDSP Operations

Wraps [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) unary vector operations.

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

Wraps [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) reduction functions.

| Function | Description | Apple function |
|----------|-------------|----------------|
| `maximum(X)`, `minimum(X)` | Max/min value | [`vDSP_maxv`](https://developer.apple.com/documentation/accelerate/vdsp_maxv), [`vDSP_minv`](https://developer.apple.com/documentation/accelerate/vdsp_minv) |
| `findmax(X)`, `findmin(X)` | Max/min value and index | [`vDSP_maxvi`](https://developer.apple.com/documentation/accelerate/vdsp_maxvi), [`vDSP_minvi`](https://developer.apple.com/documentation/accelerate/vdsp_minvi) |
| `sum(X)`, `mean(X)` | Sum and mean | [`vDSP_sve`](https://developer.apple.com/documentation/accelerate/vdsp_sve), [`vDSP_meanv`](https://developer.apple.com/documentation/accelerate/vdsp_meanv) |
| `meanmag(X)` | Mean of absolute values | [`vDSP_meamgv`](https://developer.apple.com/documentation/accelerate/vdsp_meamgv) |
| `meansqr(X)` | Mean of squares | [`vDSP_measqv`](https://developer.apple.com/documentation/accelerate/vdsp_measqv) |
| `meanssqr(X)` | Mean of signed squares | [`vDSP_mvessq`](https://developer.apple.com/documentation/accelerate/vdsp_mvessq) |
| `summag(X)` | Sum of absolute values | [`vDSP_svemg`](https://developer.apple.com/documentation/accelerate/vdsp_svemg) |
| `sumsqr(X)` | Sum of squares | [`vDSP_svesq`](https://developer.apple.com/documentation/accelerate/vdsp_svesq) |
| `sumssqr(X)` | Sum of signed squares | [`vDSP_svs`](https://developer.apple.com/documentation/accelerate/vdsp_svs) |
| [`dot`](@ref AppleAccelerate.dot) | Dot product: `sum(X .* Y)` | [`vDSP_dotpr`](https://developer.apple.com/documentation/accelerate/vdsp_dotpr) |
| [`dotpr2`](@ref AppleAccelerate.dotpr2) | Dual dot product: one `B` dotted against two `A0`/`A1` | [`vDSP_dotpr2`](https://developer.apple.com/documentation/accelerate/vdsp_dotpr2) |
| [`distancesq`](@ref AppleAccelerate.distancesq) | Squared Euclidean distance: `sum((X .- Y).^2)` | [`vDSP_distancesq`](https://developer.apple.com/documentation/accelerate/vdsp_distancesq) |
| [`rmsqv`](@ref AppleAccelerate.rmsqv) | Root mean square: `sqrt(sum(X.^2)/N)` |
| [`sve_svesq`](@ref AppleAccelerate.sve_svesq) | Simultaneous sum and sum-of-squares |
| [`maxmgv`](@ref AppleAccelerate.maxmgv) | Maximum magnitude: `max(\|X\|)` |
| [`minmgv`](@ref AppleAccelerate.minmgv) | Minimum magnitude: `min(\|X\|)` |
| [`maxmgvi`](@ref AppleAccelerate.maxmgvi) | Maximum magnitude with index |
| [`minmgvi`](@ref AppleAccelerate.minmgvi) | Minimum magnitude with index |

```@example array
X = randn(Float64, 10_000)

# Reductions
s = AppleAccelerate.sum(X)
mx = AppleAccelerate.maximum(X)
val, idx = AppleAccelerate.findmax(X)
avg = AppleAccelerate.mean(X)
nothing # hide
```

```@docs
AppleAccelerate.maximum
AppleAccelerate.minimum
AppleAccelerate.sum
AppleAccelerate.mean
AppleAccelerate.findmax
AppleAccelerate.findmin
AppleAccelerate.meanmag
AppleAccelerate.meansqr
AppleAccelerate.meanssqr
AppleAccelerate.summag
AppleAccelerate.sumsqr
AppleAccelerate.sumssqr
AppleAccelerate.dot
AppleAccelerate.dotpr2
AppleAccelerate.distancesq
AppleAccelerate.rmsqv
AppleAccelerate.sve_svesq
AppleAccelerate.maxmgv
AppleAccelerate.minmgv
AppleAccelerate.maxmgvi
AppleAccelerate.minmgvi
```

## Vector-Vector Arithmetic

| Function | Description | Apple function |
|----------|-------------|----------------|
| [`vadd`](@ref AppleAccelerate.vadd) / [`vadd!`](@ref AppleAccelerate.vadd!) | Element-wise addition | [`vDSP_vadd`](https://developer.apple.com/documentation/accelerate/vdsp_vadd) |
| [`vsub`](@ref AppleAccelerate.vsub) / [`vsub!`](@ref AppleAccelerate.vsub!) | Element-wise subtraction | [`vDSP_vsub`](https://developer.apple.com/documentation/accelerate/vdsp_vsub) |
| [`vmul`](@ref AppleAccelerate.vmul) / [`vmul!`](@ref AppleAccelerate.vmul!) | Element-wise multiplication | [`vDSP_vmul`](https://developer.apple.com/documentation/accelerate/vdsp_vmul) |
| [`vdiv`](@ref AppleAccelerate.vdiv) / [`vdiv!`](@ref AppleAccelerate.vdiv!) | Element-wise division | [`vDSP_vdiv`](https://developer.apple.com/documentation/accelerate/vdsp_vdiv) |

```@example array
A = randn(Float64, 1000)
B = randn(Float64, 1000)

# Vector arithmetic
C = AppleAccelerate.vadd(A, B)   # A .+ B
D = AppleAccelerate.vmul(A, B)   # A .* B

# Compound operation: A * scalar + B
E = AppleAccelerate.vsma(A, 2.5, B)  # A .* 2.5 .+ B
nothing # hide
```

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

| Function | Description |
|----------|-------------|
| [`vmax`](@ref AppleAccelerate.vmax) | Element-wise maximum |
| [`vmin`](@ref AppleAccelerate.vmin) | Element-wise minimum |
| [`vmaxmg`](@ref AppleAccelerate.vmaxmg) | Element-wise maximum magnitude |
| [`vminmg`](@ref AppleAccelerate.vminmg) | Element-wise minimum magnitude |
| [`vdist`](@ref AppleAccelerate.vdist) | Element-wise Euclidean distance |
| [`vtmerg`](@ref AppleAccelerate.vtmerg) | Tapered merge of two vectors |

```@docs
AppleAccelerate.vmax
AppleAccelerate.vmin
AppleAccelerate.vmaxmg
AppleAccelerate.vminmg
AppleAccelerate.vdist
AppleAccelerate.vtmerg
```

## Vector-Scalar Operations

| Function | Description | Apple function |
|----------|-------------|----------------|
| [`vsadd`](@ref AppleAccelerate.vsadd) / [`vsadd!`](@ref AppleAccelerate.vsadd!) | Vector + scalar | [`vDSP_vsadd`](https://developer.apple.com/documentation/accelerate/vdsp_vsadd) |
| [`vssub`](@ref AppleAccelerate.vssub) / [`vssub!`](@ref AppleAccelerate.vssub!) | Vector - scalar | [`vDSP_vsadd`](https://developer.apple.com/documentation/accelerate/vdsp_vsadd) |
| [`svsub`](@ref AppleAccelerate.svsub) / [`svsub!`](@ref AppleAccelerate.svsub!) | Scalar - vector | [`vDSP_vsadd`](https://developer.apple.com/documentation/accelerate/vdsp_vsadd) |
| [`vsmul`](@ref AppleAccelerate.vsmul) / [`vsmul!`](@ref AppleAccelerate.vsmul!) | Vector * scalar | [`vDSP_vsmul`](https://developer.apple.com/documentation/accelerate/vdsp_vsmul) |
| [`vsdiv`](@ref AppleAccelerate.vsdiv) / [`vsdiv!`](@ref AppleAccelerate.vsdiv!) | Vector / scalar | [`vDSP_vsdiv`](https://developer.apple.com/documentation/accelerate/vdsp_vsdiv) |
| [`svdiv`](@ref AppleAccelerate.svdiv) | Scalar / vector | [`vDSP_svdiv`](https://developer.apple.com/documentation/accelerate/vdsp_svdiv) |

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

| Function | Description |
|----------|-------------|
| [`vam`](@ref AppleAccelerate.vam) | `(A + B) * C` |
| [`vsbm`](@ref AppleAccelerate.vsbm) | `(A - B) * C` |
| [`vma`](@ref AppleAccelerate.vma) | `A * B + C` |
| [`vmsb`](@ref AppleAccelerate.vmsb) | `A * B - C` |
| [`venvlp`](@ref AppleAccelerate.venvlp) | Signal envelope |

### Four-vector operations

| Function | Description |
|----------|-------------|
| [`vaam`](@ref AppleAccelerate.vaam) | `(A + B) * (C + D)` |
| [`vsbsbm`](@ref AppleAccelerate.vsbsbm) | `(A - B) * (C - D)` |
| [`vasbm`](@ref AppleAccelerate.vasbm) | `(A + B) * (C - D)` |
| [`vmma`](@ref AppleAccelerate.vmma) | `A * B + C * D` |
| [`vmmsb`](@ref AppleAccelerate.vmmsb) | `A * B - C * D` |
| [`vpythg`](@ref AppleAccelerate.vpythg) | Pythagorean distance |

### Vector-vector-scalar operations

| Function | Description |
|----------|-------------|
| [`vasm`](@ref AppleAccelerate.vasm) | `(A + B) * c` |
| [`vsbsm`](@ref AppleAccelerate.vsbsm) | `(A - B) * c` |
| [`vsma`](@ref AppleAccelerate.vsma) | `A * b + C` |
| [`vsmsa`](@ref AppleAccelerate.vsmsa) | `A * b + c` |
| [`vmsa`](@ref AppleAccelerate.vmsa) | `A * B + c` |
| [`vsmsb`](@ref AppleAccelerate.vsmsb) | `A * b - C` |
| [`vsmsma`](@ref AppleAccelerate.vsmsma) | `A * b + C * d` |

### Dual output

| Function | Description |
|----------|-------------|
| [`vaddsub`](@ref AppleAccelerate.vaddsub) | Simultaneous add and subtract: returns `(A .+ B, A .- B)` |

```@docs
AppleAccelerate.vam
AppleAccelerate.vsbm
AppleAccelerate.vma
AppleAccelerate.vmsb
AppleAccelerate.venvlp
AppleAccelerate.vaam
AppleAccelerate.vsbsbm
AppleAccelerate.vasbm
AppleAccelerate.vmma
AppleAccelerate.vmmsb
AppleAccelerate.vpythg
AppleAccelerate.vasm
AppleAccelerate.vsbsm
AppleAccelerate.vsma
AppleAccelerate.vsmsa
AppleAccelerate.vmsa
AppleAccelerate.vsmsb
AppleAccelerate.vsmsma
AppleAccelerate.vaddsub
```

## Clipping & Thresholding

| Function | Description |
|----------|-------------|
| [`vclip`](@ref AppleAccelerate.vclip) | Clip values to `[low, high]` |
| [`vclipc`](@ref AppleAccelerate.vclipc) | Clip with count: returns `(clipped, nlow, nhigh)` |
| [`viclip`](@ref AppleAccelerate.viclip) | Inverted clip: pass values outside `[low, high]` |
| [`vthr`](@ref AppleAccelerate.vthr) | Threshold: keep or clamp to threshold |
| [`vthres`](@ref AppleAccelerate.vthres) | Threshold to zero |
| [`vlim`](@ref AppleAccelerate.vlim) | Test limit: `(b <= A[i]) ? c : -c` |
| [`vthrsc`](@ref AppleAccelerate.vthrsc) | Threshold with signed constant |
| [`vcmprs`](@ref AppleAccelerate.vcmprs) | Compress: gather elements where gate is nonzero |

```@docs
AppleAccelerate.vclip
AppleAccelerate.vclipc
AppleAccelerate.viclip
AppleAccelerate.vthr
AppleAccelerate.vthres
AppleAccelerate.vlim
AppleAccelerate.vthrsc
AppleAccelerate.vcmprs
```

## Type Conversion

| Function | Description |
|----------|-------------|
| [`vdouble`](@ref AppleAccelerate.vdouble) | Convert Float32 to Float64 |
| [`vsingle`](@ref AppleAccelerate.vsingle) | Convert Float64 to Float32 |

```@docs
AppleAccelerate.vdouble
AppleAccelerate.vsingle
```

## Ramp Generation

| Function | Description |
|----------|-------------|
| [`vramp`](@ref AppleAccelerate.vramp) | Generate a ramp: `start + i * step` |
| [`vrampmul`](@ref AppleAccelerate.vrampmul) | Multiply vector by a generated ramp |
| [`vrampmul2`](@ref AppleAccelerate.vrampmul2) | Stereo ramp multiply (two outputs) |
| [`vrampmuladd`](@ref AppleAccelerate.vrampmuladd) | Ramp-multiply then accumulate into an existing vector |
| [`vrampmuladd2`](@ref AppleAccelerate.vrampmuladd2) | Stereo ramp-multiply then accumulate (two outputs) |

```@docs
AppleAccelerate.vramp
AppleAccelerate.vrampmul
AppleAccelerate.vrampmul2
AppleAccelerate.vrampmuladd
AppleAccelerate.vrampmuladd2
```

## Linear Average

| Function | Description |
|----------|-------------|
| [`vavlin`](@ref AppleAccelerate.vavlin) | Weighted linear average of two vectors |

```@docs
AppleAccelerate.vavlin
```

## Integration & Running Operations

| Function | Description |
|----------|-------------|
| [`vrsum`](@ref AppleAccelerate.vrsum) | Running sum scaled by `scale` |
| [`vsimps`](@ref AppleAccelerate.vsimps) | Simpson's rule integration |
| [`vtrapz`](@ref AppleAccelerate.vtrapz) | Trapezoidal integration |
| [`vswsum`](@ref AppleAccelerate.vswsum) | Sliding window sum |
| [`vswmax`](@ref AppleAccelerate.vswmax) | Sliding window maximum |

```@docs
AppleAccelerate.vrsum
AppleAccelerate.vsimps
AppleAccelerate.vtrapz
AppleAccelerate.vswsum
AppleAccelerate.vswmax
```

## Interpolation

| Function | Description |
|----------|-------------|
| [`vintb`](@ref AppleAccelerate.vintb) | Linear interpolation: `A + t * (B - A)` |
| [`vlint`](@ref AppleAccelerate.vlint) | Linear interpolation from lookup table |
| [`vqint`](@ref AppleAccelerate.vqint) | Quadratic interpolation from lookup table |

```@docs
AppleAccelerate.vintb
AppleAccelerate.vlint
AppleAccelerate.vqint
```

## Polynomial Evaluation

| Function | Description |
|----------|-------------|
| [`vpoly`](@ref AppleAccelerate.vpoly) | Evaluate polynomial at each point |

```@docs
AppleAccelerate.vpoly
```

## Normalization

| Function | Description |
|----------|-------------|
| [`vnormalize`](@ref AppleAccelerate.vnormalize) | Normalize to zero mean and unit standard deviation |

```@docs
AppleAccelerate.vnormalize
```

## Zero Crossings

| Function | Description |
|----------|-------------|
| [`nzcros`](@ref AppleAccelerate.nzcros) | Find zero crossings |

```@docs
AppleAccelerate.nzcros
```

## Decibel Conversion

| Function | Description |
|----------|-------------|
| [`vdbcon`](@ref AppleAccelerate.vdbcon) | Convert to decibels relative to a reference |

```@docs
AppleAccelerate.vdbcon
```

## Vector Fill, Swap & Sort

| Function | Description |
|----------|-------------|
| [`vclr!`](@ref AppleAccelerate.vclr!) | Fill vector with zeros |
| [`vfill!`](@ref AppleAccelerate.vfill!) | Fill vector with scalar value |
| [`vswap!`](@ref AppleAccelerate.vswap!) | Swap two vectors in-place |
| [`vsort!`](@ref AppleAccelerate.vsort!) | Sort vector in-place |
| [`vsorti`](@ref AppleAccelerate.vsorti) | Return sort permutation (indices) |

```@docs
AppleAccelerate.vclr!
AppleAccelerate.vfill!
AppleAccelerate.vswap!
AppleAccelerate.vsort!
AppleAccelerate.vsorti
```

## Gathering & Indexing

| Function | Description |
|----------|-------------|
| [`vgathr`](@ref AppleAccelerate.vgathr) | Gather by index: `C[i] = A[B[i]]` |
| [`vgathra`](@ref AppleAccelerate.vgathra) | Gather via an array of pointers: `C[i] = A[i][1]` |
| [`vindex`](@ref AppleAccelerate.vindex) | Index with float indices |
| [`vgen`](@ref AppleAccelerate.vgen) | Generate linear ramp between two values |
| [`vgenp`](@ref AppleAccelerate.vgenp) | Piecewise linear interpolation from breakpoints |
| [`vtabi`](@ref AppleAccelerate.vtabi) | Table lookup with interpolation |

```@docs
AppleAccelerate.vgathr
AppleAccelerate.vgathra
AppleAccelerate.vindex
AppleAccelerate.vgen
AppleAccelerate.vgenp
AppleAccelerate.vtabi
```

## Matrix Operations

| Function | Description |
|----------|-------------|
| [`mmul`](@ref AppleAccelerate.mmul) | Matrix multiply: `C = A * B` |
| [`mtrans`](@ref AppleAccelerate.mtrans) | Matrix transpose: `C = Aᵀ` |
| [`mmov`](@ref AppleAccelerate.mmov) | Matrix copy (submatrix move) |

```@docs
AppleAccelerate.mmul
AppleAccelerate.mtrans
AppleAccelerate.mmov
```

## Integer Operations (Int32)

| Function | Description |
|----------|-------------|
| [`vaddi`](@ref AppleAccelerate.vaddi) | Int32 vector addition |
| [`vabsi`](@ref AppleAccelerate.vabsi) | Int32 absolute value |
| [`vfilli!`](@ref AppleAccelerate.vfilli!) | Fill Int32 vector with scalar |
| [`veqvi`](@ref AppleAccelerate.veqvi) | Int32 bitwise XNOR |
| [`vdivi`](@ref AppleAccelerate.vdivi) | Int32 vector divide: `C[i] = div(A[i], B[i])` |
| [`vsaddi`](@ref AppleAccelerate.vsaddi) | Int32 scalar add: `C[i] = A[i] + b` |
| [`vsdivi`](@ref AppleAccelerate.vsdivi) | Int32 scalar divide: `C[i] = div(A[i], b)` |

```@docs
AppleAccelerate.vaddi
AppleAccelerate.vabsi
AppleAccelerate.vfilli!
AppleAccelerate.veqvi
AppleAccelerate.vdivi
AppleAccelerate.vsaddi
AppleAccelerate.vsdivi
```

## Fixed-Point (Q1.15 / Q8.24) Operations

vDSP fixed-point kernels operate directly on integer storage: an `Int16` Q1.15
value `v` represents the real number `v / 32768`, and an `Int32` Q8.24 value `v`
represents `v / 2^24`. Results are computed and re-encoded at the same scale.

| Function | Description |
|----------|-------------|
| [`dotpr_s1_15`](@ref AppleAccelerate.dotpr_s1_15) / [`dotpr_s8_24`](@ref AppleAccelerate.dotpr_s8_24) | Fixed-point dot product |
| [`dotpr2_s1_15`](@ref AppleAccelerate.dotpr2_s1_15) / [`dotpr2_s8_24`](@ref AppleAccelerate.dotpr2_s8_24) | Fixed-point dual dot product |
| [`vrampmul_s1_15`](@ref AppleAccelerate.vrampmul_s1_15) / [`vrampmul_s8_24`](@ref AppleAccelerate.vrampmul_s8_24) | Fixed-point ramp multiply |
| [`vrampmul2_s1_15`](@ref AppleAccelerate.vrampmul2_s1_15) / [`vrampmul2_s8_24`](@ref AppleAccelerate.vrampmul2_s8_24) | Fixed-point stereo ramp multiply |
| [`vrampmuladd_s1_15`](@ref AppleAccelerate.vrampmuladd_s1_15) / [`vrampmuladd_s8_24`](@ref AppleAccelerate.vrampmuladd_s8_24) | Fixed-point ramp-multiply then accumulate |
| [`vrampmuladd2_s1_15`](@ref AppleAccelerate.vrampmuladd2_s1_15) / [`vrampmuladd2_s8_24`](@ref AppleAccelerate.vrampmuladd2_s8_24) | Fixed-point stereo ramp-multiply then accumulate |

```@example array
A = Int16[16384, -8192, 4096]   # 0.5, -0.25, 0.125 in Q1.15
B = Int16[8192, 16384, -16384]  # 0.25, 0.5, -0.5 in Q1.15
c = AppleAccelerate.dotpr_s1_15(A, B)  # fixed-point dot product, Q1.15-encoded
nothing # hide
```

```@docs
AppleAccelerate.dotpr_s1_15
AppleAccelerate.dotpr_s8_24
AppleAccelerate.dotpr2_s1_15
AppleAccelerate.dotpr2_s8_24
AppleAccelerate.vrampmul_s1_15
AppleAccelerate.vrampmul_s8_24
AppleAccelerate.vrampmul2_s1_15
AppleAccelerate.vrampmul2_s8_24
AppleAccelerate.vrampmuladd_s1_15
AppleAccelerate.vrampmuladd_s8_24
AppleAccelerate.vrampmuladd2_s1_15
AppleAccelerate.vrampmuladd2_s8_24
```

## 24-bit Packed Integer Conversion

vDSP represents packed 24-bit integers with a 3-byte, unaligned in-memory layout.
These wrappers surface that as ordinary `Int32`/`UInt32` Julia vectors holding
values in the 24-bit range (`-8388608:8388607` signed, `0:16777215` unsigned),
packing/unpacking around each call.

| Function | Description |
|----------|-------------|
| [`vflt24`](@ref AppleAccelerate.vflt24) | Packed 24-bit signed int → `Float32` |
| [`vfltu24`](@ref AppleAccelerate.vfltu24) | Packed 24-bit unsigned int → `Float32` |
| [`vfltsm24`](@ref AppleAccelerate.vfltsm24) | Packed 24-bit signed int → `Float32`, scaled by `b` |
| [`vfltsmu24`](@ref AppleAccelerate.vfltsmu24) | Packed 24-bit unsigned int → `Float32`, scaled by `b` |
| [`vsmfix24`](@ref AppleAccelerate.vsmfix24) | `Float32` scaled by `b`, truncated to packed 24-bit signed int |
| [`vsmfixu24`](@ref AppleAccelerate.vsmfixu24) | `Float32` scaled by `b`, truncated to packed 24-bit unsigned int |

```@docs
AppleAccelerate.vflt24
AppleAccelerate.vfltu24
AppleAccelerate.vfltsm24
AppleAccelerate.vfltsmu24
AppleAccelerate.vsmfix24
AppleAccelerate.vsmfixu24
```

## Type Conversion (int ↔ float)

| Direction | Functions | Description |
|-----------|-----------|-------------|
| float → signed int (truncate) | `vfix8`, `vfix16`, `vfix32` | Truncating conversion |
| float → unsigned int (truncate) | `vfixu8`, `vfixu16`, `vfixu32` | Truncating conversion |
| float → signed int (round) | `vfixr8`, `vfixr16`, `vfixr32` | Rounding conversion |
| float → unsigned int (round) | `vfixru8`, `vfixru16`, `vfixru32` | Rounding conversion |
| signed int → float | `vflt8`, `vflt16`, `vflt32` | Signed integer to float |
| unsigned int → float | `vfltu8`, `vfltu16`, `vfltu32` | Unsigned integer to float |

Every function in this family has an allocating variant and a mutating variant
`f!(C, A)` that writes into a preallocated `C`. The number in the name is the integer
bit width (`8`, `16`, `32`), so `vfix32` truncates to `Int32`, `vfltu16` converts
`UInt16` to float, and so on. Both `Float32` and `Float64` are supported.

For the mutating `f!(C, A)` forms, `C` must satisfy `length(C) ≥ length(A)`;
otherwise a `DimensionMismatch` is thrown before any elements are written.

The two directions differ in how the output type is chosen. For **float → int** the
integer type is fixed by the function name, so the allocating form is `f(A)`. For
**int → float** the float width is ambiguous, so the allocating form takes it
explicitly as `f(A, Float64)` (or `Float32`).

```@example array
X = Float64[-1.7, 0.4, 2.9]

I32 = AppleAccelerate.vfix32(X)            # truncate toward zero → Int32[-1, 0, 2]
R32 = AppleAccelerate.vfixr32(X)           # round to nearest    → Int32[-2, 0, 3]
Xf  = AppleAccelerate.vflt32(I32, Float64) # back to Float64 (target type required)

# Mutating variant writes into a preallocated output
out = Vector{Int32}(undef, length(X))
AppleAccelerate.vfix32!(out, X)
nothing # hide
```

## Image Convolution

| Function | Description |
|----------|-------------|
| [`f3x3`](@ref AppleAccelerate.f3x3) | 2D convolution with 3×3 filter |
| [`f5x5`](@ref AppleAccelerate.f5x5) | 2D convolution with 5×5 filter |
| [`imgfir`](@ref AppleAccelerate.imgfir) | General 2D image convolution |

```@docs
AppleAccelerate.f3x3
AppleAccelerate.f5x5
AppleAccelerate.imgfir
```

## Broadcasting

AppleAccelerate overrides `Base.copy` and `Base.copyto!` for `Broadcasted` objects, so that broadcasting syntax like `f.(X)` automatically uses the accelerated implementation.
