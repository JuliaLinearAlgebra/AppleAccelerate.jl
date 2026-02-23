# AppleAccelerate.jl

A Julia interface to Apple's [Accelerate framework](https://developer.apple.com/accelerate/), providing high-performance BLAS/LAPACK, vectorized math operations, DSP/FFT, and sparse linear algebra on macOS.

## Requirements

- macOS 13.4 or later
- Julia 1.10 or later

## Installation

```julia
using Pkg
Pkg.add("AppleAccelerate")
```

## Quick Start

### BLAS/LAPACK Forwarding

AppleAccelerate automatically forwards BLAS and LAPACK calls to Apple's Accelerate framework on package load:

```julia
using AppleAccelerate
using LinearAlgebra

# All BLAS/LAPACK operations now use Accelerate
A = randn(1000, 1000)
F = lu(A)  # Uses Accelerate LAPACK
```

### Vectorized Math

AppleAccelerate provides accelerated element-wise math operations via Apple's vecLib:

```julia
using AppleAccelerate

X = randn(Float64, 10_000)

# Accelerated math functions (not exported to avoid conflicts)
Y = AppleAccelerate.exp(X)
Y = AppleAccelerate.sin(X)
Y = AppleAccelerate.log(X)

# Or replace Base functions for transparent acceleration
AppleAccelerate.@replaceBase sin cos exp log sqrt
sin.(X)  # Now uses Accelerate
```

### FFT

```julia
using AppleAccelerate

x = randn(ComplexF64, 1024)
setup = AppleAccelerate.plan_fft(length(x))
y = AppleAccelerate.fft(x, setup)
```

### Sparse Linear Algebra

```julia
using AppleAccelerate, SparseArrays

A = sprandn(100, 100, 0.1)
A = A + A' + 20I  # Make symmetric positive definite
b = randn(100)

import AppleAccelerate: AAFactorization, solve
f = AAFactorization(A)
x = solve(f, b)
```
