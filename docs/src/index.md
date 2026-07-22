# Introduction

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate)
for macOS, providing accelerated:

- **dense linear algebra** — BLAS/LAPACK, forwarded automatically via LBT,
- **sparse linear algebra** (`libSparse`) — direct and iterative solvers, real and complex,
- **vectorized elementwise math** (vForce / vDSP),
- **signal processing** (vDSP) — FFT, DCT, convolution, biquad filtering, window functions,
- **complex-vector operations** (split-complex vDSP),
- **neural-network primitives** (BNNS), and
- **image processing** (vImage) — geometry, convolution, morphology, histogram, alpha
  compositing, and format/colorspace conversion,

plus **numerical integration** (Quadrature).

AppleAccelerate is built as **two layers**: an auto-generated raw ABI mirror of
Accelerate's C API (`AppleAccelerate.LibAccelerate`) and a hand-written idiomatic
API on top of it. AppleAccelerate ships no package extensions and defines no
methods on other packages' functions — loading it never changes what code
elsewhere in your session does. See [Architecture](@ref architecture).

!!! note "Namespace"
    Most functions are **not exported** to avoid conflicts with Base and LinearAlgebra. Access them via the `AppleAccelerate.` prefix (e.g., `AppleAccelerate.exp(X)`). The exception is dense linear algebra, which is activated automatically via LBT on package load.

!!! warning "Key Limitations"
    - **macOS only** — this package uses Apple's Accelerate framework, which is not available on Linux or Windows.
    - **FFT size limits** — 1D complex `fft`/`ifft`/`bfft` support lengths `f * 2^k` with `f ∈ {1, 3, 5, 15}`; 2D transforms and the real FFT (`rfft`/`brfft`) are power-of-2 only. DFT also supports `f * 2^n` where `f ∈ {1, 3, 5, 15}`. Use FFTW.jl for other sizes or N-D transforms.
    - **Precision** — the vectorized array/DSP ops are `Float32`/`Float64` only; BNNS is `Float32`-centric; sparse solvers also accept `ComplexF32`/`ComplexF64`; vImage covers the usual 8-/16-bit and float pixel formats.

## Installation

Requires macOS 13.4+ and Julia 1.10+.

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

## Quick Start

One example per major subsystem — the same set as the
[README](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl). Every function lives
under the `AppleAccelerate.` prefix; the package intentionally exports nothing, so it never
shadows `Base`/`LinearAlgebra`.

### Dense linear algebra (BLAS/LAPACK via LBT)

```@example quickstart-dense
using AppleAccelerate, LinearAlgebra
A = randn(1000, 1000)
F = lu(A)                                       # BLAS/LAPACK routed to Accelerate
nothing # hide
```

### Vectorized elementwise math (vForce / vDSP)

```@example quickstart-math
using AppleAccelerate
X = randn(10_000)
Y = AppleAccelerate.exp(X)                      # also sin, cos, log, sqrt, tanh, …
AppleAccelerate.sincos(X)                       # fused, both results in one pass
nothing # hide
```

### Signal processing — FFT / DCT / convolution / biquad (vDSP)

```@example quickstart-fft
using AppleAccelerate
x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x)                      # cached setup; also rfft, fft2d, dct
nothing # hide
```

### Complex vector operations (split-complex vDSP)

```@example quickstart-complex
using AppleAccelerate
z = randn(ComplexF64, 1000)
mags = AppleAccelerate.vmags(z)                 # squared magnitudes (abs2)
ang  = AppleAccelerate.vphase(z)                # phase angles
nothing # hide
```

### Sparse solvers — direct & iterative (libSparse)

```@example quickstart-sparse
using AppleAccelerate, LinearAlgebra, SparseArrays
S = sprandn(500, 500, 0.01); S = S*S' + 500I    # symmetric positive-definite
As = AppleAccelerate.AASparseMatrix(SparseMatrixCSC{Float64,Int64}(S))
xs = AppleAccelerate.solve(AppleAccelerate.cholesky(As), randn(500))
nothing # hide
```

### Neural-network primitives — a tiny 2-layer MLP (BNNS)

```@example quickstart-bnns
using AppleAccelerate
W1, b1 = randn(Float32, 16, 8), randn(Float32, 16)   # layer 1: 8 → 16
W2, b2 = randn(Float32, 4, 16), randn(Float32, 4)    # layer 2: 16 → 4
x = randn(Float32, 8)
h = AppleAccelerate.bnns_activation(:relu, AppleAccelerate.bnns_matmul(W1, reshape(x, :, 1)) .+ b1)
logits = AppleAccelerate.bnns_matmul(W2, h) .+ b2     # 4-class output
nothing # hide
```

### Image processing (vImage)

```@example quickstart-vimage
using AppleAccelerate
img   = rand(Float32, 64, 48)                        # a 64×48 planar (grayscale) image
small = AppleAccelerate.scale_PlanarF(img, 32, 24)   # resize to 32×24
flip  = AppleAccelerate.horizontalReflect_PlanarF(img)
nothing # hide
```
