# AppleAccelerate.jl

[![CI](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/)
[![JuliaHub](https://juliahub.com/docs/General/AppleAccelerate/stable/version.svg)](https://juliahub.com/ui/Packages/General/AppleAccelerate)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia compat](https://img.shields.io/badge/Julia-≥1.10-blue.svg)](https://julialang.org/downloads/)

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate), providing:

- **Vectorized array operations** via vecLib (`vv*`) and vDSP (`vDSP_*`) — element-wise math, reductions, compound arithmetic, clipping, interpolation, integration, and more — **2–19× faster** than Base Julia for transcendental functions (`sin`, `cos`, `exp`, `log`)
- **Dense linear algebra** via BLAS/LAPACK forwarding through [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) — all standard `LinearAlgebra` operations (`lu`, `qr`, `svd`, `cholesky`, `eigen`, etc.) are accelerated transparently — **6–14× faster** GEMM than OpenBLAS on Apple Silicon, **2–4× faster** factorizations and solves
- **Sparse linear algebra** via `libSparse` — sparse matrix operations and direct solvers (QR, Cholesky, LDLT) — **faster for Float32 QR** and **Cholesky at N=5000** vs SuiteSparse
- **Signal Processing** — FFT, DCT, convolution, cross-correlation, biquad filtering, window functions — **on par with FFTW** for complex FFT, **up to 2× faster** for Float32 real FFT

## Requirements

- macOS 13.4 or later
- Julia 1.10 or later

## Installation

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
```

## Documentation

See the [full documentation](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/) for the complete API reference.
