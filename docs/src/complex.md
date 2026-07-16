# Complex Array Operations

AppleAccelerate wraps [vDSP](https://developer.apple.com/documentation/accelerate/vdsp)'s split-complex (`vDSP_z*`) functions for `Vector{Complex{Float32}}` and `Vector{Complex{Float64}}`. These extend existing function names (e.g., `vneg`, `vabs`, `vmul`) with methods that dispatch on complex element types — no naming conflicts with the real-valued versions on the [Array Operations](array.md) page.

Julia's interleaved complex storage is passed to the split-complex vDSP routines directly (using a stride-2 split-complex view of the same memory), so the mutating (`!`) variants are in-place and allocation-free.

These functions are **not exported**. Access them via the `AppleAccelerate.` prefix.

```@setup complex
using AppleAccelerate
```

## Complex unary operations

| Function | Description |
|----------|-------------|
| [`vneg`](@ref AppleAccelerate.vneg)`(X)` / `vneg!(result, X)` | Negate: `-X` |
| [`vabs`](@ref AppleAccelerate.vabs)`(X)` / `vabs!(result, X)` | Modulus: `abs.(X)` |
| [`vconj`](@ref AppleAccelerate.vconj) | Complex conjugate |
| [`vcopy`](@ref AppleAccelerate.vcopy) | Copy via split-complex move |

The complex-valued docstrings for `vneg` and `vabs` are rendered alongside the real-valued ones on the [Array Operations](array.md) page.

## Complex → real operations

| Function | Description |
|----------|-------------|
| [`vphase`](@ref AppleAccelerate.vphase) | Complex phase (angle) |
| [`vmags`](@ref AppleAccelerate.vmags) | Squared magnitude (`abs2`) |
| [`vmagsa`](@ref AppleAccelerate.vmagsa) | Squared magnitude + accumulate |

## Complex binary operations

| Function | Description |
|----------|-------------|
| [`vmul`](@ref AppleAccelerate.vmul)`(X, Y)` / `vmul!(result, X, Y)` | Element-wise multiply: `X .* Y` |
| [`vdiv`](@ref AppleAccelerate.vdiv)`(X, Y)` / `vdiv!(result, X, Y)` | Element-wise divide: `X ./ Y` |
| [`vsmul`](@ref AppleAccelerate.vsmul)`(X, c)` / `vsmul!(result, X, c)` | Scalar multiply (complex scalar) |
| [`dot`](@ref AppleAccelerate.dot)`(X, Y)` | Conjugated dot product matching `LinearAlgebra.dot`: `sum(conj(X) .* Y)` |
| [`dotu`](@ref AppleAccelerate.dotu)`(X, Y)` | Unconjugated (bilinear) dot product: `sum(X .* Y)` |
| [`zvadd`](@ref AppleAccelerate.zvadd) | Complex addition: `A + B` |
| [`zvsub`](@ref AppleAccelerate.zvsub) | Complex subtraction: `A - B` |
| [`zvcmul`](@ref AppleAccelerate.zvcmul) | Conjugate multiply: `conj(A) * B` |

## Complex-real operations

| Function | Description |
|----------|-------------|
| [`zrvmul`](@ref AppleAccelerate.zrvmul) | Complex × real |
| [`zrvdiv`](@ref AppleAccelerate.zrvdiv) | Complex / real |
| [`zrvadd`](@ref AppleAccelerate.zrvadd) | Complex + real (adds to real part) |
| [`zrvsub`](@ref AppleAccelerate.zrvsub) | Complex − real |

## Complex compound operations

| Function | Description |
|----------|-------------|
| [`zvcma`](@ref AppleAccelerate.zvcma) | `conj(A)*B + C` |
| [`zvma`](@ref AppleAccelerate.zvma) | `A*B + C` |
| [`zvsma`](@ref AppleAccelerate.zvsma) | `A*b + C` (b is complex scalar) |

## Complex dot products

| Function | Description |
|----------|-------------|
| [`zidotpr`](@ref AppleAccelerate.zidotpr) | Conjugate dot: `sum(conj(A) .* B)` |
| [`zrdotpr`](@ref AppleAccelerate.zrdotpr) | Complex-real dot: `sum(A .* B)` |

## Complex fill & convolution

| Function | Description |
|----------|-------------|
| [`zvfill!`](@ref AppleAccelerate.zvfill!) | Fill complex vector with scalar |
| [`zconv`](@ref AppleAccelerate.zconv) | Complex convolution |
| [`zmmul`](@ref AppleAccelerate.zmmul) | Complex matrix multiply |

## Coordinate conversion

| Function | Description |
|----------|-------------|
| [`polar`](@ref AppleAccelerate.polar) | Cartesian to polar coordinates |
| [`rect`](@ref AppleAccelerate.rect) | Polar to Cartesian coordinates |

## Format conversion

| Function | Description |
|----------|-------------|
| [`ctoz`](@ref AppleAccelerate.ctoz) | Interleaved complex → split (real, imag) vectors |
| [`ztoc`](@ref AppleAccelerate.ztoc) | Split (real, imag) vectors → interleaved complex |

```@example complex
Z = AppleAccelerate.cis(randn(Float64, 100))  # complex array

# Complex operations
conj_Z = AppleAccelerate.vconj(Z)
phases = AppleAccelerate.vphase(Z)
mags = AppleAccelerate.vmags(Z)

# Coordinate conversion
r, θ = AppleAccelerate.polar(Z)
Z_back = AppleAccelerate.rect(r, θ)
nothing # hide
```

```@docs
AppleAccelerate.vconj
AppleAccelerate.vcopy
AppleAccelerate.vphase
AppleAccelerate.vmags
AppleAccelerate.vmagsa
AppleAccelerate.dotu
AppleAccelerate.polar
AppleAccelerate.rect
AppleAccelerate.zvadd
AppleAccelerate.zvsub
AppleAccelerate.zvcmul
AppleAccelerate.zrvmul
AppleAccelerate.zrvdiv
AppleAccelerate.zrvadd
AppleAccelerate.zrvsub
AppleAccelerate.zvcma
AppleAccelerate.zvma
AppleAccelerate.zvsma
AppleAccelerate.zidotpr
AppleAccelerate.zrdotpr
AppleAccelerate.zvfill!
AppleAccelerate.zconv
AppleAccelerate.zmmul
AppleAccelerate.ctoz
AppleAccelerate.ztoc
```
