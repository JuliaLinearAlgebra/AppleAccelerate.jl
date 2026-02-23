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

```@example
using AppleAccelerate
using LinearAlgebra

# All BLAS/LAPACK operations now use Accelerate
A = randn(1000, 1000)
F = lu(A)  # Uses Accelerate LAPACK
nothing # hide
```

### Vectorized Math

AppleAccelerate provides accelerated element-wise math operations via Apple's vecLib:

```@example
using AppleAccelerate

X = randn(Float64, 10_000)

# Accelerated math functions (not exported to avoid conflicts)
Y = AppleAccelerate.exp(X)
Y = AppleAccelerate.sin(X)
Y = AppleAccelerate.log(X)

# Or replace Base functions for transparent acceleration
AppleAccelerate.@replaceBase sin cos exp log sqrt
sin.(X)  # Now uses Accelerate
nothing # hide
```

### FFT

```@example
using AppleAccelerate

x = randn(ComplexF64, 1024)
X = AppleAccelerate.fft(x)
x_recovered = AppleAccelerate.ifft(X)
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
