# AppleAccelerate.jl

[![CI](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/)

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate), providing:

- **Dense linear algebra** via BLAS/LAPACK forwarding through [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) — all standard `LinearAlgebra` operations (`lu`, `qr`, `svd`, `cholesky`, `eigen`, etc.) are accelerated transparently
- **Vectorized array operations** via vecLib (`vv*`) and vDSP (`vDSP_*`) — element-wise math, reductions, compound arithmetic, clipping, interpolation, integration, and more
- **DSP & FFT** — FFT, DCT, convolution, cross-correlation, biquad filtering, window functions
- **Sparse linear algebra** via `libSparse` — sparse matrix operations and direct solvers (QR, Cholesky, LDLT)

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

# All BLAS/LAPACK operations now use Accelerate
A = randn(1000, 1000)
F = lu(A)

# Accelerated vectorized math
X = randn(10_000)
Y = AppleAccelerate.exp(X)
Y = AppleAccelerate.sin(X)

# Or replace Base functions for transparent acceleration
AppleAccelerate.@replaceBase sin cos exp log sqrt
sin.(X)  # Now uses Accelerate
```

## Documentation

See the [full documentation](https://JuliaLinearAlgebra.github.io/AppleAccelerate.jl/dev/) for the complete API reference, including all available vDSP functions, DSP operations, and sparse linear algebra.
