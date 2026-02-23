# AppleAccelerate.jl

[![CI](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/)

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate), providing:

- **Vectorized array operations** via vecLib (`vv*`) and vDSP (`vDSP_*`) — element-wise math, reductions, compound arithmetic, clipping, interpolation, integration, and more
- **Dense linear algebra** via BLAS/LAPACK forwarding through [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) — all standard `LinearAlgebra` operations (`lu`, `qr`, `svd`, `cholesky`, `eigen`, etc.) are accelerated transparently
- **Sparse linear algebra** via `libSparse` — sparse matrix operations and direct solvers (QR, Cholesky, LDLT)
- **Signal Processing** — FFT, DCT, convolution, cross-correlation, biquad filtering, window functions

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
