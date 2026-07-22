# FFT & Transforms

AppleAccelerate wraps Apple's [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) transform functions for FFT, real FFT, DFT, and DCT.

```@setup fft
using AppleAccelerate
```

## FFT

The FFT API wraps Apple's [vDSP FFT functions](https://developer.apple.com/documentation/accelerate/fast_fourier_transforms) and follows the same naming conventions as [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl). Both 1D vectors and 2D matrices are supported with `ComplexF64` and `ComplexF32` inputs. A **1D** `fft`/`ifft`/`bfft`/`rfft`/`irfft` (called without an explicit `plan_fft` setup) accepts any power-of-two length, plus non-power-of-two lengths of the form `f * 2^k` with `f ∈ {3, 5, 15}` — for complex transforms `k ≥ 3` (so the smallest such lengths are 24, 40, 120), and for real transforms `k ≥ 4`. Powers of two use vDSP's fast FFT path, while the mixed-radix lengths transparently use Apple's DFT (see [`is_supported_fft_length`](@ref AppleAccelerate.is_supported_fft_length)). **2D transforms remain power-of-2 only.** Unsupported lengths throw an `ArgumentError`.

!!! note "Choosing between AppleAccelerate and FFTW"
    AppleAccelerate's FFT is a direct vDSP binding: 1-D transforms support powers of two plus `f * 2^k` lengths (`f ∈ {3, 5, 15}`, with `k ≥ 3` for complex and `k ≥ 4` for real transforms), while 2-D transforms are power-of-2 only. Only 1-D vectors and 2-D matrices are supported, and it is reached exclusively through the `AppleAccelerate.` prefix. It does **not** plug into the AbstractFFTs.jl interface — loading AppleAccelerate never changes what `fft`/`plan_fft` do in your session, and will never override FFTW.jl.

    For general sizes, 3-D and higher arrays, dimension-selective (`region`) transforms, or the standard Julia FFT interface, use [FFTW.jl](https://github.com/JuliaMath/FFTW.jl). vDSP is mainly competitive at small transform sizes; FFTW is typically faster at larger ones.

    Real FFT (`rfft`, `irfft`, `brfft`) is likewise available only through the `AppleAccelerate.` prefix.

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

### Batched 1D FFT

Passing a `dims` argument (1 = each column, 2 = each row) computes a batch of
independent 1D transforms over a matrix in a single call, like FFTW's
`fft(A, dims)`. This wraps [`vDSP_fftm_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fftm_zop). The transform length (`size(A, dims)`) must be a power of 2:

```@example fft
A = randn(ComplexF64, 16, 8)
F = AppleAccelerate.fft(A, 1)     # FFT of each column
A1 = AppleAccelerate.ifft(F, 1)   # normalized inverse
@assert A1 ≈ A
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

### Temp-buffer variants

For repeated transforms of the same size, an [`FFTWorkspace`](@ref AppleAccelerate.FFTWorkspace) supplies vDSP with a reusable scratch buffer (the `*t` vDSP variants), avoiding its internal per-call allocation. Pass it as the last argument to `fft`/`fft!`/`bfft`/`ifft` (1D, 2D, or batched) and to the real transforms. Results are identical to the non-workspace calls.

```@example fft
x = randn(ComplexF64, 1024)
setup = AppleAccelerate.plan_fft(x)
ws = AppleAccelerate.FFTWorkspace(x)
X = AppleAccelerate.fft(x, setup, ws)
@assert AppleAccelerate.ifft(X, setup, ws) ≈ x
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
AppleAccelerate.FFTWorkspace
AppleAccelerate.is_supported_fft_length
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

1D real transforms called without an explicit setup also accept non-power-of-2 lengths of the form `f * 2^k` with `f ∈ {3, 5, 15}` and `k ≥ 4` (e.g. 48, 160, 240), via Apple's real-input mixed-radix DFT ([`vDSP_DFT_zrop_CreateSetup`](https://developer.apple.com/documentation/accelerate/vdsp_dft_zrop_createsetup)).

2D real FFT is supported for power-of-2 dimensions and matches FFTW's `rfft` layout — an `n1×n2` real matrix transforms to an `(n1÷2+1)×n2` complex matrix:

```@example fft
x2 = randn(Float64, 16, 32)
X2 = AppleAccelerate.rfft(x2)                # 9×32 complex matrix
x2_recovered = AppleAccelerate.irfft(X2, 16) # back to real, 16×32
@assert x2_recovered ≈ x2
nothing # hide
```

Batched real FFT (`rfft(x, dims)` / `irfft(X, n, dims)`) transforms each column (`dims = 1`) or row (`dims = 2`) of a real matrix independently, like FFTW's `rfft(x, dims)`. In-place real variants (`rfft!`) reuse the input buffer as vDSP's scratch (overwriting it), and workspace variants accept an [`FFTWorkspace`](@ref AppleAccelerate.FFTWorkspace).

```@docs
AppleAccelerate.plan_rfft
AppleAccelerate.rfft
AppleAccelerate.rfft!
AppleAccelerate.irfft
AppleAccelerate.brfft
```

## Small-radix, fixed-size, and interleaved transforms

Beyond the general FFT/DFT paths, vDSP provides specialized complex-transform kernels:

- **Radix-3 / radix-5 FFTs** ([`fftradix3`](@ref AppleAccelerate.fftradix3), [`fftradix5`](@ref AppleAccelerate.fftradix5), and their unnormalized inverses `bfftradix3`/`bfftradix5`) handle lengths `3·2^k` and `5·2^k` — including `12`, `20`, `24`, `40`, which the mixed-radix `fft`/`dft` path cannot.
- **Fixed-size 16-/32-point FFTs** ([`fft16`](@ref AppleAccelerate.fft16), [`fft32`](@ref AppleAccelerate.fft32), inverses `bfft16`/`bfft32`) are dedicated `Float32` kernels that need no setup.
- **Interleaved-complex DFT** ([`dft_interleaved`](@ref AppleAccelerate.dft_interleaved) / [`idft_interleaved`](@ref AppleAccelerate.idft_interleaved), planned with [`plan_dft_interleaved`](@ref AppleAccelerate.plan_dft_interleaved)) operates directly on `Vector{Complex}` buffers rather than the split-complex layout of [`dft`](@ref AppleAccelerate.dft).

```@example fft
x = randn(ComplexF64, 12)           # 12 = 3·2^2, not an fft() length
@assert AppleAccelerate.bfftradix3(AppleAccelerate.fftradix3(x)) ≈ 12 .* x

y = randn(ComplexF32, 16)
@assert AppleAccelerate.bfft16(AppleAccelerate.fft16(y)) ≈ 16 .* y

z = randn(ComplexF64, 24)
@assert AppleAccelerate.idft_interleaved(AppleAccelerate.dft_interleaved(z)) ≈ z
nothing # hide
```

```@docs
AppleAccelerate.fftradix3
AppleAccelerate.fft16
AppleAccelerate.plan_dft_interleaved
AppleAccelerate.dft_interleaved
AppleAccelerate.idft_interleaved
```

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
| [`idct`](@ref AppleAccelerate.idct) | Compute the inverse Discrete Cosine Transform |
| [`plan_destroy`](@ref AppleAccelerate.plan_destroy) | Destroy a DCT setup object |

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
AppleAccelerate.idct
AppleAccelerate.plan_destroy
```
