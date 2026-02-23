# DSP & FFT

AppleAccelerate wraps vDSP functions for signal processing, including convolution, correlation, filtering, windowing, DCT, and FFT.

## FFT

```julia
# Create an FFT plan (reusable for same-size transforms)
setup = AppleAccelerate.plan_fft(1024, Float64)

# Forward FFT
x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x, setup)

# Inverse FFT
x_recovered = AppleAccelerate.fft(X, setup, AppleAccelerate.FFT_INVERSE)
```

Both `ComplexF64` and `ComplexF32` inputs are supported. Input length must be a power of 2.

## DCT (Discrete Cosine Transform)

Float32 only. Supports DCT types II, III, and IV.

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
```

## Convolution

`conv(X, K)` computes the convolution of vectors `X` and `K`. `conv!(result, X, K)` stores the result in a preallocated vector.

```julia
X = randn(Float64, 100)
K = randn(Float64, 10)
result = AppleAccelerate.conv(X, K)  # length = length(X) + length(K) - 1
```

## Cross-Correlation

`xcorr(X, Y)` computes the cross-correlation. `xcorr(X)` computes auto-correlation.

```julia
X = randn(Float64, 100)
Y = randn(Float64, 100)
result = AppleAccelerate.xcorr(X, Y)
```

## Biquad Filtering

IIR filtering via cascaded biquad sections (Float64 only).

```julia
coefficients = randn(5)  # 5 coefficients per section
bq = AppleAccelerate.biquadcreate(coefficients, 1)
delays = zeros(4)
result = AppleAccelerate.biquad(X, delays, length(X), bq)
```

## Window Functions

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
AppleAccelerate.hann
```
