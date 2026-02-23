# DSP & FFT

AppleAccelerate wraps vDSP functions for signal processing, including convolution, correlation, filtering, windowing, DCT, and FFT.

```@setup dsp
using AppleAccelerate, AbstractFFTs
```

## FFT

The FFT API follows the same naming conventions as [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl). Both 1D vectors and 2D matrices are supported with `ComplexF64` and `ComplexF32` inputs. All dimensions must be powers of 2.

```@example dsp
x = randn(ComplexF64, 1024)

# Forward FFT (auto-creates plan)
X = AppleAccelerate.fft(x)

# Normalized inverse FFT
x_recovered = AppleAccelerate.ifft(X)

# Unnormalized inverse (backward) FFT
X_back = AppleAccelerate.bfft(X)

# Reusable plan for repeated transforms of the same size
setup = AppleAccelerate.plan_fft(x)
X = AppleAccelerate.fft(x, setup)
x_recovered = AppleAccelerate.ifft(X, setup)
nothing # hide
```

2D FFT works the same way — `fft` dispatches on the input shape:

```@example dsp
x2 = randn(ComplexF64, 16, 32)
X2 = AppleAccelerate.fft(x2)
x2_recovered = AppleAccelerate.ifft(X2)
nothing # hide
```

### In-place FFT

In-place variants avoid allocating the output array:

```@example dsp
x = randn(ComplexF64, 1024)
y = copy(x)
AppleAccelerate.fft!(y)       # forward, modifies y in-place
AppleAccelerate.ifft!(y)      # inverse, modifies y in-place
@assert y ≈ x
nothing # hide
```

### Real FFT

`rfft` computes the FFT of real input, returning only the non-redundant complex coefficients (length `N÷2+1`). `irfft` inverts it back to real output.

```@example dsp
x = randn(Float64, 1024)
X = AppleAccelerate.rfft(x)          # Complex vector of length 513
x_recovered = AppleAccelerate.irfft(X, 1024)  # Back to real, length 1024
@assert x_recovered ≈ x
nothing # hide
```

### AbstractFFTs Integration

AppleAccelerate provides a package extension for [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl), allowing you to use the standard Julia FFT interface backed by Apple's vDSP library. Load both packages to activate the extension:

```@example dsp
x = randn(ComplexF64, 1024)

# Out-of-place
X = fft(x)
x_recovered = ifft(X)

# In-place
y = copy(x)
fft!(y)
ifft!(y)
@assert y ≈ x

# Plan-based API (works for all variants)
p = plan_fft(x)
X = p * x
x_recovered = p \ X
nothing # hide
```

The extension supports the full AbstractFFTs interface for 1D and 2D:

| Function | Description |
|----------|-------------|
| `fft`, `ifft`, `bfft` | Out-of-place complex FFT |
| `fft!`, `ifft!`, `bfft!` | In-place complex FFT |
| `plan_fft`, `plan_ifft`, `plan_bfft` | Out-of-place plans |
| `plan_fft!`, `plan_bfft!` | In-place plans |
| `inv(plan)`, `p \\ x`, `mul!` | Plan operations |

Input must be `Vector` or `Matrix` with `Complex{T}` elements and power-of-2 dimensions. Non-power-of-2 inputs, empty arrays, and 3D+ arrays will throw an `ArgumentError` with a message suggesting to use FFTW.jl.

!!! note
    Real FFT (`rfft`, `irfft`, `brfft`) is available only via the direct `AppleAccelerate.*` API and is not wired into the AbstractFFTs extension, to avoid dispatch conflicts with other packages that call `rfft` internally on non-power-of-2 inputs.

## DCT (Discrete Cosine Transform)

Float32 only. Supports DCT types II, III, and IV.

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
```

## Convolution

`conv(X, K)` computes the convolution of vectors `X` and `K`. `conv!(result, X, K)` stores the result in a preallocated vector.

```@example dsp
X = randn(Float64, 100)
K = randn(Float64, 10)
result = AppleAccelerate.conv(X, K)  # length = length(X) + length(K) - 1
nothing # hide
```

## Cross-Correlation

`xcorr(X, Y)` computes the cross-correlation. `xcorr(X)` computes auto-correlation.

```@example dsp
X = randn(Float64, 100)
Y = randn(Float64, 100)
result = AppleAccelerate.xcorr(X, Y)
nothing # hide
```

## Biquad Filtering

IIR filtering via cascaded biquad sections (Float64 only).

```@example dsp
X = randn(Float64, 100)
coefficients = randn(5)  # 5 coefficients per section
bq = AppleAccelerate.biquadcreate(coefficients, 1)
delays = zeros(4)
result = AppleAccelerate.biquad(X, delays, length(X), bq)
nothing # hide
```

## Window Functions

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
AppleAccelerate.hann
```
