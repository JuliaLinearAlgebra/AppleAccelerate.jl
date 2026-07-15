# Introduction

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate), providing high-performance BLAS/LAPACK, vectorized math operations, DSP/FFT, and sparse linear algebra on macOS.

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
| FFT & Transforms | [vDSP FFT](https://developer.apple.com/documentation/accelerate/fast_fourier_transforms) | Complex/real FFT (1D/2D, power-of-2), DFT, DCT — direct vDSP bindings under the `AppleAccelerate.` namespace |
| Filtering & Spectral | [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) | Convolution, biquad IIR, spectral analysis, window functions |
| Sparse Linear Algebra | [Sparse Solvers](https://developer.apple.com/documentation/accelerate/sparse_solvers) | SpMV, direct solvers (QR, Cholesky, LDLT), factorization reuse (`refactor!`) |
| Neural Network Primitives | [BNNS](https://developer.apple.com/documentation/accelerate/bnns) | `Float32` matmul and activations — see [BNNS](@ref "Neural Network Primitives (BNNS)") |
| Numerical Integration | [Quadrature](https://developer.apple.com/documentation/accelerate/quadrature) | Adaptive QUADPACK integration — see [Numerical Integration](@ref) |

AppleAccelerate is built as **two layers**: an auto-generated raw ABI mirror of
Accelerate's C API (`AppleAccelerate.LibAccelerate`) and a hand-written idiomatic
API on top of it. AppleAccelerate ships no package extensions and defines no
methods on other packages' functions — loading it never changes what code
elsewhere in your session does. See [Architecture](@ref architecture).

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

### Neural Network Primitives

BNNS-backed `Float32` matrix multiply and activations (see [Neural Network Primitives (BNNS)](@ref)):

```@example
using AppleAccelerate

A = rand(Float32, 4, 3); B = rand(Float32, 3, 5)
C = AppleAccelerate.bnns_matmul(A, B)
y = AppleAccelerate.bnns_activation(:relu, Float32[-1, 0, 1])
nothing # hide
```

### Integration

Adaptive QUADPACK integration (see [Numerical Integration](@ref)):

```@example
using AppleAccelerate

AppleAccelerate.integrate(x -> exp(-x^2), -Inf, Inf).value  # ≈ √π
```
