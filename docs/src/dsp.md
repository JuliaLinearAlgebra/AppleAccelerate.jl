# DSP & FFT

AppleAccelerate wraps vDSP functions for signal processing, including convolution, correlation, filtering, windowing, DCT, and FFT.

```@setup dsp
using AppleAccelerate
```

## FFT

```@example dsp
# Create an FFT plan (reusable for same-size transforms)
setup = AppleAccelerate.plan_fft(1024, Float64)

# Forward FFT
x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x, setup)

# Inverse FFT
x_recovered = AppleAccelerate.fft(X, setup, AppleAccelerate.FFT_INVERSE)
nothing # hide
```

Both `ComplexF64` and `ComplexF32` inputs are supported. Input length must be a power of 2.

## 2D FFT

```@example dsp
# Create an FFT plan (must accommodate the largest dimension)
setup = AppleAccelerate.plan_fft(32, Float64)

# Forward 2D FFT
x = randn(ComplexF64, 16, 32)
X = AppleAccelerate.fft2d(x, setup)

# Inverse 2D FFT (unnormalized — divide by nrows*ncols to recover original)
x_recovered = AppleAccelerate.fft2d(X, setup, AppleAccelerate.FFT_INVERSE) ./ (16 * 32)
nothing # hide
```

Both dimensions must be powers of 2. The `FFTSetup` must be created with `n >= max(nrows, ncols)`.

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
