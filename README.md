<p align="center">
  <img src="docs/src/assets/logo.png" alt="AppleAccelerate.jl" width="400"/>
</p>

<p align="center">
  <a href="https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml/badge.svg" alt="CI"/></a>
  <a href="https://codecov.io/gh/JuliaLinearAlgebra/AppleAccelerate.jl"><img src="https://codecov.io/gh/JuliaLinearAlgebra/AppleAccelerate.jl/branch/master/graph/badge.svg" alt="Coverage"/></a>
  <a href="https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Docs"/></a>
  <a href="https://juliahub.com/ui/Packages/General/AppleAccelerate"><img src="https://juliahub.com/docs/General/AppleAccelerate/stable/version.svg" alt="JuliaHub"/></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"/></a>
  <a href="https://julialang.org/downloads/"><img src="https://img.shields.io/badge/Julia-≥1.10-blue.svg" alt="Julia compat"/></a>
</p>

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate), providing:

- **Vectorized array operations** via vDSP and vForce — element-wise math, reductions, compound arithmetic, clipping, interpolation — **2–19× faster** than Base Julia for transcendentals (`sin`, `cos`, `exp`, `log`)
- **Dense linear algebra** — all of `LinearAlgebra` (`lu`, `qr`, `svd`, `cholesky`, `eigen`, …) accelerated transparently via [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) — **6–13× faster** single-threaded GEMM than OpenBLAS on Apple Silicon (SME/AMX co-processor), plus **2–4× faster** factorizations and solves
- **Sparse linear algebra** via `libSparse` — direct (Cholesky / LDLᵀ / LU / QR) and iterative (CG / GMRES / LSMR) solvers, real and complex
- **Signal processing** — 1D/2D real & complex FFT (batched, mixed-radix), DCT, convolution, biquad filtering, window functions; cached setups make no-plan `fft(x)` competitive with FFTW and drop the FFTW dependency
- **Neural-network primitives** via BNNS — `Float32` matrix multiply and pointwise activations
- **Image processing** via vImage — geometry (scale, rotate, affine warp), convolution, morphology, histogram, alpha compositing, and format/colorspace conversion (incl. Y′CbCr)

See the [benchmarks](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/benchmarks/) for full performance comparisons and methodology.

## Installation

Requires macOS 13.4+ and Julia 1.10+.

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

## Quick start

One self-contained, copy-pasteable example per subsystem. Every function lives under the
`AppleAccelerate.` prefix — the package intentionally exports nothing, so it never shadows
`Base`/`LinearAlgebra`.

### Dense linear algebra — all of `LinearAlgebra` accelerated transparently via LBT

```julia
using AppleAccelerate, LinearAlgebra
A = randn(1000, 1000)
F = lu(A)                                       # BLAS/LAPACK routed to Accelerate
```

### Vectorized elementwise math (vForce / vDSP)

```julia
using AppleAccelerate
X = randn(10_000)
Y = AppleAccelerate.exp(X)                      # also sin, cos, log, sqrt, tanh, …
AppleAccelerate.sincos(X)                       # fused, both results in one pass
```

### Signal processing — FFT / DCT / convolution / biquad filtering

```julia
using AppleAccelerate
x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x)                      # cached setup; also rfft, fft2d, dct
```

### Complex vector operations (split-complex vDSP)

```julia
using AppleAccelerate
z = randn(ComplexF64, 1000)
mags = AppleAccelerate.vmags(z)                 # squared magnitudes (abs2)
ang  = AppleAccelerate.vphase(z)                # phase angles
```

### Sparse solvers (libSparse) — direct Cholesky / LDLᵀ / LU / QR and iterative CG / GMRES / LSMR

```julia
using AppleAccelerate, LinearAlgebra, SparseArrays
S = sprandn(500, 500, 0.01); S = S*S' + 500I    # symmetric positive-definite
As = AppleAccelerate.AASparseMatrix(SparseMatrixCSC{Float64,Int64}(S))
xs = AppleAccelerate.solve(AppleAccelerate.cholesky(As), randn(500))
```

### Neural-network primitives (BNNS) — a tiny 2-layer MLP forward pass

```julia
using AppleAccelerate
W1, b1 = randn(Float32, 16, 8), randn(Float32, 16)   # layer 1: 8 → 16
W2, b2 = randn(Float32, 4, 16), randn(Float32, 4)    # layer 2: 16 → 4
x = randn(Float32, 8)
h = AppleAccelerate.bnns_activation(:relu, AppleAccelerate.bnns_matmul(W1, reshape(x, :, 1)) .+ b1)
logits = AppleAccelerate.bnns_matmul(W2, h) .+ b2     # 4-class output
```

### Image processing (vImage)

```julia
using AppleAccelerate
img   = rand(Float32, 64, 48)                        # a 64×48 planar (grayscale) image
small = AppleAccelerate.scale_PlanarF(img, 32, 24)   # resize to 32×24
flip  = AppleAccelerate.horizontalReflect_PlanarF(img)
```

## Documentation

See the [full documentation](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/) for the complete API reference.
