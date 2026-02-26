# Introduction

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/accelerate/), providing high-performance BLAS/LAPACK, vectorized math operations, DSP/FFT, and sparse linear algebra on macOS.

## Installation

Requires macOS 13.4+ and Julia 1.10+.

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

## Feature Overview

| Feature | Backend | Details |
|---------|---------|---------|
| Dense Linear Algebra | [BLAS](https://developer.apple.com/documentation/accelerate/blas) / [LAPACK](https://developer.apple.com/documentation/accelerate/solving-systems-of-linear-equations-with-lapack) | Automatic forwarding via LBT — `lu`, `qr`, `svd`, `mul!`, etc. |
| Array Operations | [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) / [vecLib](https://developer.apple.com/documentation/accelerate/veclib) | Element-wise math (`exp`, `sin`, `log`), reductions, vector arithmetic |
| FFT & Transforms | [vDSP FFT](https://developer.apple.com/documentation/accelerate/fast_fourier_transforms) | Complex/real FFT (1D/2D), DFT, DCT with AbstractFFTs integration |
| Filtering & Spectral | [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) | Convolution, biquad IIR, spectral analysis, window functions |
| Sparse Linear Algebra | [Sparse Solvers](https://developer.apple.com/documentation/accelerate/sparse_solvers) | SpMV, direct solvers (QR, Cholesky, LDLT) |

!!! note "Namespace"
    Most functions are **not exported** to avoid conflicts with Base and LinearAlgebra. Access them via the `AppleAccelerate.` prefix (e.g., `AppleAccelerate.exp(X)`). The exception is dense linear algebra, which is activated automatically via LBT on package load.

!!! warning "Key Limitations"
    - **macOS only** — this package uses Apple's Accelerate framework, which is not available on Linux or Windows.
    - **FFT requires power-of-2 sizes** — use FFTW.jl for arbitrary sizes. DFT supports `f * 2^n` where `f ∈ {1, 3, 5, 15}`.
    - **Float32/Float64 only** — no support for integer or extended-precision types.

## Quick Start

### Array Operations

AppleAccelerate provides accelerated element-wise math operations via Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib):

```@example
using AppleAccelerate

X = randn(Float64, 10_000)

# Accelerated math functions (not exported to avoid conflicts)
Y = AppleAccelerate.exp(X)
Y = AppleAccelerate.sin(X)
Y = AppleAccelerate.log(X)

# Broadcasting works automatically
Y = AppleAccelerate.sin.(X)
nothing # hide
```

### Dense Linear Algebra

AppleAccelerate automatically forwards BLAS and LAPACK calls to Apple's Accelerate framework on package load:

```@example
using AppleAccelerate
using LinearAlgebra

# All BLAS/LAPACK operations now use Accelerate
A = randn(1000, 1000)
F = lu(A)  # Uses Accelerate LAPACK
nothing # hide
```

### Sparse Linear Algebra

```@example
using AppleAccelerate, SparseArrays, LinearAlgebra

A = sprandn(100, 100, 0.1)
A = A + A' + 20I  # Make symmetric positive definite
b = randn(100)

import AppleAccelerate: AAFactorization, solve
f = AAFactorization(A)
x = solve(f, b)
nothing # hide
```

### Signal Processing

```@example
using AppleAccelerate

x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x)
x_recovered = AppleAccelerate.ifft(X)
nothing # hide
```
