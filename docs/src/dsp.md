# DSP & FFT

AppleAccelerate wraps vDSP functions for signal processing, including convolution, correlation, filtering, windowing, DCT, and FFT.

```@setup dsp
using AppleAccelerate, AbstractFFTs
```

## FFT

The FFT API follows the same naming conventions as [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl): `fft`, `ifft`, `bfft`, and `plan_fft`. Both 1D vectors and 2D matrices are supported with `ComplexF64` and `ComplexF32` inputs. All dimensions must be powers of 2.

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

### AbstractFFTs Integration

AppleAccelerate provides a package extension for [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl), allowing you to use the standard Julia FFT interface backed by Apple's vDSP library. Load both packages to activate the extension:

```@example dsp
# 1D FFT
x = randn(ComplexF64, 1024)
X = fft(x)
x_recovered = ifft(X)

# 2D FFT
x2 = randn(ComplexF64, 16, 32)
X2 = fft(x2)
x2_recovered = ifft(X2)

# Plan-based API (works for both 1D and 2D)
p = plan_fft(x)
X = p * x
x_recovered = p \ X
nothing # hide
```

The extension supports `fft`, `ifft`, `bfft`, `plan_fft`, `plan_ifft`, `plan_bfft`, plan inversion via `inv(plan)`, and `mul!` for both 1D vectors and 2D matrices. Input must be `Vector{Complex{T}}` or `Matrix{Complex{T}}` with power-of-2 dimensions.

!!! note
    In-place plans (`plan_fft!`, `plan_bfft!`) and real FFTs (`rfft`, `irfft`) are not supported since vDSP only provides out-of-place complex FFTs. 3D and higher-dimensional arrays are not supported by vDSP and will fall through to FFTW if available.

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
