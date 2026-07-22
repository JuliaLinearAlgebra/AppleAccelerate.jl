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

- **Vectorized array operations** via vecLib (`vv*`) and vDSP (`vDSP_*`) — element-wise math, reductions, compound arithmetic, clipping, interpolation, integration, and more — **2–19× faster** than Base Julia for transcendental functions (`sin`, `cos`, `exp`, `log`)
- **Dense linear algebra** via BLAS/LAPACK forwarding through [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) — all standard `LinearAlgebra` operations (`lu`, `qr`, `svd`, `cholesky`, `eigen`, etc.) are accelerated transparently — **6–13× faster** GEMM than OpenBLAS on Apple Silicon when both are single-threaded; OpenBLAS narrows the gap with multiple threads, since Accelerate offloads GEMM to the SME/AMX matrix co-processor whose throughput is largely thread-independent while OpenBLAS scales across CPU cores — plus **2–4× faster** factorizations and solves
- **Sparse linear algebra** via `libSparse` — sparse matrix operations and direct solvers (QR, Cholesky, LDLT) — **faster for Float32 QR** and **Cholesky at N=5000** vs SuiteSparse
- **Signal Processing** — 1D/2D complex and real FFT (including batched and mixed-radix lengths), DCT, convolution, cross-correlation, biquad filtering, window functions — with pre-planned transforms FFTW is generally faster (roughly 1.5–14× across complex 1D, real, and 2D FFT), but for no-plan convenience calls AppleAccelerate's cached FFT setups make `fft(x)` **faster than FFTW at small sizes and competitive for Float32 throughout**; vDSP's FFT also avoids an FFTW dependency when Accelerate is already loaded
- **Neural-network primitives** via BNNS — `Float32` matrix multiply and pointwise activations (`:relu`, `:sigmoid`, `:tanh`, `:abs`, `:identity`)

## Installation

Requires macOS 13.4+ and Julia 1.10+.

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

## Quick start

An example from each major subsystem (every function lives under the `AppleAccelerate.`
prefix — the package intentionally exports nothing, so it never shadows `Base`/`LinearAlgebra`):

```julia
using AppleAccelerate
using LinearAlgebra, SparseArrays

# --- Dense linear algebra: all of LinearAlgebra is accelerated transparently via LBT ---
A = randn(1000, 1000)
F = lu(A)                                      # BLAS/LAPACK routed to Accelerate

# --- Vectorized elementwise math (vForce / vDSP) ---
X = randn(10_000)
Y = AppleAccelerate.exp(X)                     # also sin, cos, log, sqrt, tanh, …
AppleAccelerate.sincos(X)                      # fused, both results in one pass

# --- Signal processing: FFT / DCT / convolution / biquad filtering ---
x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x)                      # cached setup; also rfft, fft2d, dct

# --- Complex vector operations (split-complex vDSP) ---
z = randn(ComplexF64, 1000)
mags = AppleAccelerate.vmags(z)                 # squared magnitudes (abs2)
ang  = AppleAccelerate.vphase(z)                # phase angles

# --- Sparse direct & iterative solvers (libSparse): Cholesky / LDLᵀ / LU / QR / CG / GMRES ---
S = sprandn(500, 500, 0.01); S = S*S' + 500I    # symmetric positive-definite
As = AppleAccelerate.AASparseMatrix(SparseMatrixCSC{Float64,Int64}(S))
xs = AppleAccelerate.solve(AppleAccelerate.cholesky(As), randn(500))

# --- Neural-network primitives (BNNS) ---
W = rand(Float32, 4, 3); v = rand(Float32, 3, 8)
C = AppleAccelerate.bnns_matmul(W, v)           # α·(W*v), Float32
```

## Documentation

See the [full documentation](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/) for the complete API reference.
