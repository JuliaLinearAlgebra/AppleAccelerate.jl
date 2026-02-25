# Signal Processing

AppleAccelerate wraps [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) functions for signal processing, including convolution, correlation, filtering, windowing, DCT, and FFT.

```@setup dsp
using AppleAccelerate, AbstractFFTs
```

## FFT

The FFT API wraps Apple's [vDSP FFT functions](https://developer.apple.com/documentation/accelerate/fast_fourier_transforms) and follows the same naming conventions as [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl). Both 1D vectors and 2D matrices are supported with `ComplexF64` and `ComplexF32` inputs. All dimensions must be powers of 2.

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
X = AppleAccelerate.fft(x)
x_recovered = AppleAccelerate.ifft(X)
@assert x_recovered ≈ x

# In-place
y = copy(x)
AppleAccelerate.fft!(y)
AppleAccelerate.ifft!(y)
@assert y ≈ x
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

## FFT Benchmarks

Pre-planned FFT performance comparing Apple vDSP against FFTW, both single-threaded. Run the full benchmark suite with `julia --project=test/bench test/bench/run_benchmarks.jl fft` ([source](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl)).

**Complex 1D FFT** — vDSP and FFTW are closely matched. vDSP is notably faster at N=4096; at larger sizes FFTW has a slight edge:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| ComplexF64 | 1,024 | 7.0 | 7.1 | 1.01× |
| ComplexF64 | 4,096 | 20.5 | 34.1 | 1.67× |
| ComplexF64 | 65,536 | 718 | 653 | 0.91× |
| ComplexF64 | 1,048,576 | 18,767 | 17,208 | 0.92× |
| ComplexF32 | 1,024 | 3.1 | 3.1 | 1.01× |
| ComplexF32 | 4,096 | 11.5 | 18.4 | 1.61× |
| ComplexF32 | 65,536 | 329 | 323 | 0.98× |
| ComplexF32 | 1,048,576 | 4,891 | 4,753 | 0.97× |

**Real FFT** — For Float32, vDSP is faster at larger sizes (up to 2.2×). Float64 rfft favors FFTW:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| Float64 | 1,024 | 3.8 | 1.6 | 0.42× |
| Float64 | 65,536 | 360 | 209 | 0.58× |
| Float32 | 4,096 | 7.0 | 7.6 | 1.09× |
| Float32 | 65,536 | 115 | 169 | 1.47× |
| Float32 | 262,144 | 539 | 1,166 | 2.16× |

**Complex 2D FFT** — Nearly identical performance:

| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |
|------|------|-----------|-----------|---------|
| ComplexF64 | 64×64 | 47 | 47 | 1.00× |
| ComplexF64 | 256×256 | 4,659 | 4,708 | 1.01× |
| ComplexF32 | 64×64 | 18 | 22 | 1.19× |
| ComplexF32 | 256×256 | 321 | 325 | 1.01× |

!!! note "Benchmark environment"
    Apple M2 Max, macOS 26, single-threaded. Julia 1.12.5, AppleAccelerate v0.6.0, FFTW v1.10.0. Times are minimum of 5 trials. Both use pre-planned transforms. Run [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl) to reproduce.

## DCT (Discrete Cosine Transform)

Wraps Apple's [vDSP DCT functions](https://developer.apple.com/documentation/accelerate/discrete_cosine_transforms). Float32 only. Supports DCT types II, III, and IV.

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
AppleAccelerate.plan_destroy
```

## Convolution

Wraps [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv). `conv(X, K)` computes the convolution of vectors `X` and `K`. `conv!(result, X, K)` stores the result in a preallocated vector.

```@example dsp
X = randn(Float64, 100)
K = randn(Float64, 10)
result = AppleAccelerate.conv(X, K)  # length = length(X) + length(K) - 1
nothing # hide
```

## Cross-Correlation

Wraps [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv) (with reversed kernel). `xcorr(X, Y)` computes the cross-correlation. `xcorr(X)` computes auto-correlation.

```@example dsp
X = randn(Float64, 100)
Y = randn(Float64, 100)
result = AppleAccelerate.xcorr(X, Y)
nothing # hide
```

## Biquad Filtering

Wraps [`vDSP_biquad`](https://developer.apple.com/documentation/accelerate/vdsp_biquad). IIR filtering via cascaded biquad sections (Float64 only).

```@example dsp
X = randn(Float64, 100)
coefficients = randn(5)  # 5 coefficients per section
bq = AppleAccelerate.biquadcreate(coefficients, 1)
delays = zeros(4)
result = AppleAccelerate.biquad(X, delays, length(X), bq)
nothing # hide
```

## Window Functions

Wraps Apple's [vDSP window generation functions](https://developer.apple.com/documentation/accelerate/vdsp/vector_generation).

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
AppleAccelerate.hann
```
