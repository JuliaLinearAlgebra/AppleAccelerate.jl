# FFT & Transforms

AppleAccelerate wraps Apple's [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) transform functions for FFT, real FFT, DFT, and DCT.

```@setup fft
using AppleAccelerate, AbstractFFTs
```

## FFT

The FFT API wraps Apple's [vDSP FFT functions](https://developer.apple.com/documentation/accelerate/fast_fourier_transforms) and follows the same naming conventions as [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl). Both 1D vectors and 2D matrices are supported with `ComplexF64` and `ComplexF32` inputs. All dimensions must be powers of 2.

```@example fft
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

```@example fft
x2 = randn(ComplexF64, 16, 32)
X2 = AppleAccelerate.fft(x2)
x2_recovered = AppleAccelerate.ifft(X2)
nothing # hide
```

### In-place FFT

In-place variants avoid allocating the output array:

```@example fft
x = randn(ComplexF64, 1024)
y = copy(x)
AppleAccelerate.fft!(y)       # forward, modifies y in-place
AppleAccelerate.ifft!(y)      # inverse, modifies y in-place
@assert y ≈ x
nothing # hide
```

```@docs
AppleAccelerate.plan_fft
AppleAccelerate.fft
AppleAccelerate.ifft
AppleAccelerate.bfft
AppleAccelerate.fft!
AppleAccelerate.ifft!
AppleAccelerate.bfft!
```

## Real FFT

`rfft` computes the FFT of real input, returning only the non-redundant complex coefficients (length `N÷2+1`). `irfft` inverts it back to real output.

```@example fft
x = randn(Float64, 1024)
X = AppleAccelerate.rfft(x)          # Complex vector of length 513
x_recovered = AppleAccelerate.irfft(X, 1024)  # Back to real, length 1024
@assert x_recovered ≈ x
nothing # hide
```

```@docs
AppleAccelerate.plan_rfft
AppleAccelerate.rfft
AppleAccelerate.irfft
AppleAccelerate.brfft
```

## AbstractFFTs Integration

AppleAccelerate provides a package extension for [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl), allowing you to use the standard Julia FFT interface backed by Apple's vDSP library. Load both packages to activate the extension:

```@example fft
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
| `inv(plan)`, `p \ x`, `mul!` | Plan operations |

Input must be `Vector` or `Matrix` with `Complex{T}` elements and power-of-2 dimensions. Non-power-of-2 inputs, empty arrays, and 3D+ arrays will throw an `ArgumentError` with a message suggesting to use FFTW.jl.

### Region argument

For 2D transforms, the `region` argument selects which dimensions to transform — matching the AbstractFFTs convention:

```@example fft
x = randn(ComplexF64, 16, 32)

p_cols = plan_fft(x, (1,))   # transform along columns only (dim 1)
p_rows = plan_fft(x, (2,))   # transform along rows only (dim 2)
p_both = plan_fft(x, (1,2))  # transform along both (default)

# Equivalent to two 1D FFTs composed:
X_both = p_both * x
X_seq  = p_rows * (p_cols * x)
@assert X_both ≈ X_seq
nothing # hide
```

!!! note
    Real FFT (`rfft`, `irfft`, `brfft`) is available only via the direct `AppleAccelerate.*` API and is not wired into the AbstractFFTs extension, to avoid dispatch conflicts with other packages that call `rfft` internally on non-power-of-2 inputs.

## DFT (Complex Discrete Fourier Transform)

Wraps Apple's [`vDSP_DFT_zop`](https://developer.apple.com/documentation/accelerate/vdsp_dft_execute) for complex-to-complex DFT. Unlike the FFT functions, DFT supports non-power-of-2 lengths of the form `f * 2^n` where `f ∈ {1, 3, 5, 15}` and `n ≥ 3`. Both `Float32` and `Float64` are supported.

```@example fft
x = randn(ComplexF64, 120)  # 120 = 15 * 2^3, non-power-of-2

# Forward DFT (auto-creates setup)
X = AppleAccelerate.dft(x)

# Normalized inverse DFT
x_recovered = AppleAccelerate.idft(X)
@assert x_recovered ≈ x

# Reusable setup for repeated transforms
setup_fwd = AppleAccelerate.plan_dft(120, AppleAccelerate.DFT_FORWARD, Float64)
setup_inv = AppleAccelerate.plan_dft(120, AppleAccelerate.DFT_INVERSE, Float64)
X = AppleAccelerate.dft(x, setup_fwd)
x_recovered = AppleAccelerate.idft(X, setup_inv)
nothing # hide
```

```@docs
AppleAccelerate.plan_dft
AppleAccelerate.dft
AppleAccelerate.idft
```

## DCT (Discrete Cosine Transform)

Wraps Apple's [vDSP DCT functions](https://developer.apple.com/documentation/accelerate/discrete_cosine_transforms). Float32 only. Supports DCT types II, III, and IV.

| Function | Description |
|----------|-------------|
| [`plan_dct`](@ref AppleAccelerate.plan_dct) | Create a DCT setup object |
| [`dct`](@ref AppleAccelerate.dct) | Compute the Discrete Cosine Transform |
| [`plan_destroy`](@ref AppleAccelerate.plan_destroy) | Destroy a DCT setup object |

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
AppleAccelerate.plan_destroy
```
