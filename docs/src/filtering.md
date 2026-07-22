# Filtering & Spectral Analysis

AppleAccelerate wraps [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) functions for convolution, correlation, filtering, spectral analysis, and window generation.

```@setup filtering
using AppleAccelerate
```

## Convolution

Wraps [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv). `conv(X, K)` computes the convolution of vectors `X` and `K`. `conv!(result, X, K)` stores the result in a preallocated vector.

```@example filtering
X = randn(Float64, 100)
K = randn(Float64, 10)
result = AppleAccelerate.conv(X, K)  # length = length(X) + length(K) - 1
nothing # hide
```

```@docs
AppleAccelerate.conv
AppleAccelerate.conv!
```

## Cross-Correlation

Wraps [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv) (with reversed kernel). `xcorr(X, Y)` computes the cross-correlation. `xcorr(X)` computes auto-correlation.

```@example filtering
X = randn(Float64, 100)
Y = randn(Float64, 100)
result = AppleAccelerate.xcorr(X, Y)
nothing # hide
```

```@docs
AppleAccelerate.xcorr
AppleAccelerate.xcorr!
```

## Biquad Filtering

Wraps [`vDSP_biquad`](https://developer.apple.com/documentation/accelerate/vdsp_biquad). IIR filtering via cascaded biquad sections. Both `Float64` and `Float32` processing are supported; coefficients are always `Float64`.

Each section takes 5 coefficients in the order `[b0, b1, b2, a1, a2]` (feedforward
then feedback). Use real filter-design coefficients rather than random values — an
arbitrary `a1`/`a2` can place the poles outside the unit circle and make the filter
unstable. The example below is a stable second-order lowpass:

```@example filtering
X = randn(Float64, 100)
coefficients = [0.2929, 0.5858, 0.2929, 0.0, 0.1716]  # [b0, b1, b2, a1, a2]
bq = AppleAccelerate.biquadcreate(coefficients, 1)  # defaults to Float64
delays = zeros(4)
result = AppleAccelerate.biquad(X, delays, length(X), bq)

# Float32 processing
X32 = randn(Float32, 100)
bq32 = AppleAccelerate.biquadcreate(coefficients, 1, Float32)
delays32 = zeros(Float32, 4)
result32 = AppleAccelerate.biquad(X32, delays32, length(X32), bq32)
nothing # hide
```

A single-precision setup's coefficients can be updated in place with
[`biquad_setcoefficients!`](@ref AppleAccelerate.biquad_setcoefficients!), avoiding
a full recreate when re-tuning a filter (vDSP provides this setter only for the
`Float32` setup).

```@docs
AppleAccelerate.biquadcreate
AppleAccelerate.biquad
AppleAccelerate.biquad_setcoefficients!
AppleAccelerate.biquaddestroy
```

## Multi-Channel Biquad Filtering

Wraps [`vDSP_biquadm`](https://developer.apple.com/documentation/accelerate/vdsp_biquadm). Applies a multi-channel IIR filter in a single call. Both `Float32` (default) and `Float64` are supported.

```@example filtering
channels = 2
sections = 1
coefficients = randn(5 * channels * sections)

setup = AppleAccelerate.biquadm_create(coefficients, channels, sections)  # Float32
inputs = [randn(Float32, 64) for _ in 1:channels]
outputs = AppleAccelerate.biquadm(inputs, 64, setup)

# Float64 variant
setup64 = AppleAccelerate.biquadm_create(coefficients, channels, sections, Float64)
inputs64 = [randn(Float64, 64) for _ in 1:channels]
outputs64 = AppleAccelerate.biquadm(inputs64, 64, setup64)
nothing # hide
```

### Live-state controls

A `BiquadMulti` setup keeps internal per-channel/section state and coefficients
that can be manipulated between processing calls — for gain automation,
click-free parameter changes, resetting or transplanting filter memory, and
bypassing sections:

| Function | Effect |
|----------|--------|
| [`biquadm_setcoefficients!`](@ref AppleAccelerate.biquadm_setcoefficients!) | Immediately set a block of `(section, channel)` coefficients |
| [`biquadm_settargets!`](@ref AppleAccelerate.biquadm_settargets!) | Set target coefficients the filter interpolates toward (smooth changes) |
| [`biquadm_resetstate!`](@ref AppleAccelerate.biquadm_resetstate!) | Zero all internal delays |
| [`biquadm_copystate!`](@ref AppleAccelerate.biquadm_copystate!) | Copy internal state from one setup to another |
| [`biquadm_setactivefilters!`](@ref AppleAccelerate.biquadm_setactivefilters!) | Enable/disable (bypass) individual sections |

```@example filtering
setup = AppleAccelerate.biquadm_create([1.0,0,0,0,0, 1.0,0,0,0,0], 2, 1, Float64)
# Retune channel 1 (0-based) to a gain of 5 without recreating the setup:
AppleAccelerate.biquadm_setcoefficients!(setup, [5.0,0,0,0,0], 0, 1, 1, 1)
y = AppleAccelerate.biquadm([ones(8), ones(8)], 8, setup)
@assert y[2] ≈ fill(5.0, 8)
AppleAccelerate.biquadm_resetstate!(setup)  # clear filter memory
nothing # hide
```

```@docs
AppleAccelerate.biquadm_create
AppleAccelerate.biquadm
AppleAccelerate.biquadm_setcoefficients!
AppleAccelerate.biquadm_settargets!
AppleAccelerate.biquadm_resetstate!
AppleAccelerate.biquadm_copystate!
AppleAccelerate.biquadm_setactivefilters!
```

## Recursive Filters

### Second-Order Recursive Filter (deq22)

Wraps [`vDSP_deq22`](https://developer.apple.com/documentation/accelerate/vdsp_deq22) for second-order (two-pole two-zero) recursive filtering. Both `Float32` and `Float64` are supported.

```@example filtering
A = randn(Float64, 64)
B = Float64[0.5, -0.3, 0.2, 0.1, -0.05]  # 5 filter coefficients
C = AppleAccelerate.deq22(A, B)  # returns N output samples
nothing # hide
```

The `deq22!` variant operates on pre-padded arrays with explicit initial state.

### FIR Decimation Filter (desamp)

Wraps [`vDSP_desamp`](https://developer.apple.com/documentation/accelerate/vdsp_desamp) for FIR filtering with decimation. Both `Float32` and `Float64` are supported.

```@example filtering
A = randn(Float64, 100)
F = randn(Float64, 5)    # 5-tap FIR filter
DF = 3                    # decimation factor
C = AppleAccelerate.desamp(A, DF, F)  # output length = div(100 - 5, 3) + 1
nothing # hide
```

### Wiener-Levinson Filter (wiener)

Wraps [`vDSP_wiener`](https://developer.apple.com/documentation/accelerate/vdsp_wiener) for solving the Wiener-Hopf equation via Levinson-Durbin recursion. Both `Float32` and `Float64` are supported.

```@example filtering
L = 8
autocorr = Float64[1.0 / (1 + abs(k)) for k in 0:(L-1)]
crosscorr = randn(Float64, L)
F, err = AppleAccelerate.wiener(autocorr, crosscorr)  # err == 0 on success
nothing # hide
```

```@docs
AppleAccelerate.deq22
AppleAccelerate.deq22!
AppleAccelerate.desamp
AppleAccelerate.desamp!
AppleAccelerate.wiener
AppleAccelerate.wiener!
```

## Spectral Analysis

Wraps Apple's [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) spectral analysis functions for computing power spectra, cross-spectra, coherence, and transfer functions. All functions support both `Float32` and `Float64`.

### Autospectrum ([`vDSP_zaspec`](https://developer.apple.com/documentation/accelerate/vdsp_zaspec))

```@example filtering
x = randn(ComplexF64, 256)
power = AppleAccelerate.zaspec(x)  # |x[n]|^2
nothing # hide
```

The `zaspec!` variant accumulates into an existing vector (`C[n] += |A[n]|^2`), useful for averaging multiple frames.

### Cross-Spectrum ([`vDSP_zcspec`](https://developer.apple.com/documentation/accelerate/vdsp_zcspec))

```@example filtering
x = randn(ComplexF64, 256)
y = randn(ComplexF64, 256)
csd = AppleAccelerate.zcspec(x, y)  # conj(x[n]) * y[n], accumulated
nothing # hide
```

### Coherence ([`vDSP_zcoher`](https://developer.apple.com/documentation/accelerate/vdsp_zcoher))

```@example filtering
n = 256
Sxx = abs2.(randn(ComplexF64, n))  # autospectrum of x
Syy = abs2.(randn(ComplexF64, n))  # autospectrum of y
Sxy = randn(ComplexF64, n)          # cross-spectrum
coh = AppleAccelerate.zcoher(Sxx, Syy, Sxy)  # |Sxy|^2 / (Sxx * Syy)
nothing # hide
```

### Transfer Function ([`vDSP_ztrans`](https://developer.apple.com/documentation/accelerate/vdsp_ztrans))

```@example filtering
Sxx = abs2.(randn(ComplexF64, 256))
Sxy = randn(ComplexF64, 256)
H = AppleAccelerate.ztrans(Sxx, Sxy)  # Sxy[n] / Sxx[n]
nothing # hide
```

```@docs
AppleAccelerate.zaspec
AppleAccelerate.zaspec!
AppleAccelerate.zcspec
AppleAccelerate.zcspec!
AppleAccelerate.zcoher
AppleAccelerate.zcoher!
AppleAccelerate.ztrans
AppleAccelerate.ztrans!
```

## Window Functions

Wraps Apple's [vDSP window generation functions](https://developer.apple.com/documentation/accelerate/vdsp/vector_generation).

| Function | Description |
|----------|-------------|
| [`blackman`](@ref AppleAccelerate.blackman) | Generate a Blackman window |
| [`hamming`](@ref AppleAccelerate.hamming) | Generate a Hamming window |
| [`hanning`](@ref AppleAccelerate.hanning) | Generate a Hanning window |

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
```
