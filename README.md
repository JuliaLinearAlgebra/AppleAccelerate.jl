<p align="center">
  <img src="docs/src/assets/logo.png" alt="AppleAccelerate.jl" width="400"/>
</p>

<p align="center">
  <a href="https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml"><img src="https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml/badge.svg" alt="CI"/></a>
  <a href="https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Docs"/></a>
  <a href="https://juliahub.com/ui/Packages/General/AppleAccelerate"><img src="https://juliahub.com/docs/General/AppleAccelerate/stable/version.svg" alt="JuliaHub"/></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"/></a>
  <a href="https://julialang.org/downloads/"><img src="https://img.shields.io/badge/Julia-≥1.10-blue.svg" alt="Julia compat"/></a>
</p>

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate), providing:

- **Vectorized array operations** via vecLib (`vv*`) and vDSP (`vDSP_*`) — element-wise math, reductions, compound arithmetic, clipping, interpolation, integration, and more — **2–19× faster** than Base Julia for transcendental functions (`sin`, `cos`, `exp`, `log`)
- **Dense linear algebra** via BLAS/LAPACK forwarding through [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) — all standard `LinearAlgebra` operations (`lu`, `qr`, `svd`, `cholesky`, `eigen`, etc.) are accelerated transparently — **6–14× faster** GEMM than OpenBLAS on Apple Silicon, **2–4× faster** factorizations and solves
- **Sparse linear algebra** via `libSparse` — sparse matrix operations and direct solvers (QR, Cholesky, LDLT) — **faster for Float32 QR** and **Cholesky at N=5000** vs SuiteSparse
- **Signal Processing** — FFT, DCT, convolution, cross-correlation, biquad filtering, window functions — **on par with FFTW** for complex FFT, **up to 2× faster** for Float32 real FFT
- **SIMD math inside `@simd` loops** via `AppleAccelerate.SIMDMath` — scalar math functions that LLVM turns into SIMD calls, for loops the array API can't express (strided access, values computed on the fly) — **2–4× faster** than a scalar Base loop

## Installation

Requires macOS 13.4+ and Julia 1.10+.

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

## Quick start

```julia
using AppleAccelerate
using LinearAlgebra

# All LinearAlgebra is accelerated through LBT
A = randn(1000, 1000)
F = lu(A)

# Accelerated vectorized math
X = randn(10_000)
Y = AppleAccelerate.exp(X)
Y = AppleAccelerate.sin(X)

# FFT
x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x)

# SIMD math inside a loop the array API above can't express.
# Prefer AppleAccelerate.exp/log on whole arrays when you can -- they're faster.
using AppleAccelerate.SIMDMath: log
function logsum_strided(X, stride)
    u = zero(eltype(X))
    @simd for i in 1:stride:length(X)
        @inbounds u += log(X[i])
    end
    u
end
```

## Documentation

See the [full documentation](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/) for the complete API reference.
