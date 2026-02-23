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

Both `ComplexF64` and `ComplexF32` inputs are supported.

```@docs
AppleAccelerate.plan_fft
AppleAccelerate.fft
AppleAccelerate.destroy_fftsetup
AppleAccelerate.FFTSetup
```

## DCT (Discrete Cosine Transform)

Float32 only. Supports DCT types II, III, and IV.

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
```

## Convolution

```@docs
AppleAccelerate.conv
AppleAccelerate.conv!
```

## Cross-Correlation

```@docs
AppleAccelerate.xcorr
AppleAccelerate.xcorr!
```

## Biquad Filtering

IIR filtering via cascaded biquad sections (Float64 only).

```@docs
AppleAccelerate.biquadcreate
AppleAccelerate.biquad
```

## Window Functions

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
AppleAccelerate.hann
```
