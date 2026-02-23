# Dense Linear Algebra

AppleAccelerate forwards [BLAS](https://developer.apple.com/documentation/accelerate/blas) and [LAPACK](https://developer.apple.com/documentation/accelerate/lapack) calls to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate) via Julia's [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) (LBT) mechanism. This happens automatically when the package is loaded.

## How it works

On `__init__`, AppleAccelerate loads both LP64 and ILP64 BLAS/LAPACK interfaces from Accelerate. OpenBLAS remains as a fallback for operations not provided by Accelerate (e.g., `gemmt`).

Since Accelerate provides a full [BLAS](https://developer.apple.com/documentation/accelerate/blas) and [LAPACK](https://developer.apple.com/documentation/accelerate/lapack) implementation, all standard Julia [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) operations are accelerated transparently. This includes:

**Factorizations:** `lu`, `qr`, `cholesky`, `svd`, `eigen`, `schur`, `ldlt`, `hessenberg`, and their in-place `!` variants

**Solvers:** `\`, `ldiv!`, `rdiv!`

**Matrix operations:** `mul!`, `*`, `det`, `tr`, `inv`, `pinv`, `rank`, `norm`, `cond`, `opnorm`

**BLAS routines:** All Level 1 (vector), Level 2 (matrix-vector), and Level 3 (matrix-matrix) operations via `LinearAlgebra.BLAS`

For the complete list of available operations, see the [Julia LinearAlgebra documentation](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/).

## Loading

```@example
using AppleAccelerate
```

```@docs
AppleAccelerate.load_accelerate
```

## Threading

On macOS 26+, you can control BLAS threading:

```@docs
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
```

## Utilities

```@docs
AppleAccelerate.get_macos_version
```

## Benchmarks

Performance comparison of Apple Accelerate vs OpenBLAS. The dense benchmark script loads OpenBLAS first, then switches to Accelerate via LBT, so both are measured in the same process. Run with `julia --project=test/bench test/bench/run_benchmarks.jl dense` ([source](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl)).

### GEMM (`mul!`) — GFLOPS (higher is better)

Accelerate is 6–14× faster for matrix multiply, with the largest gains for Float32 due to AMX/NEON support:

| Type | N | OpenBLAS | Accelerate | Speedup |
|------|---|----------|------------|---------|
| Float64 | 64 | 17 | 180 | 10.6× |
| Float64 | 256 | 31 | 251 | 8.0× |
| Float64 | 1,024 | 36 | 258 | 7.2× |
| Float64 | 4,096 | 36 | 241 | 6.6× |
| Float32 | 64 | 62 | 406 | 6.5× |
| Float32 | 256 | 70 | 889 | 12.7× |
| Float32 | 1,024 | 73 | 1,029 | 14.1× |
| Float32 | 4,096 | 73 | 906 | 12.4× |

### Factorizations (Float64) — time in μs (lower is better)

| Operation | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |
|-----------|---|---------------|-----------------|---------|
| LU | 256 | 639 | 285 | 2.2× |
| LU | 1,024 | 24,639 | 15,471 | 1.6× |
| LU | 2,048 | 182,654 | 92,604 | 2.0× |
| QR | 512 | 8,193 | 4,559 | 1.8× |
| QR | 2,048 | 407,071 | 200,736 | 2.0× |
| Cholesky | 256 | 359 | 114 | 3.1× |
| Cholesky | 1,024 | 13,647 | 3,801 | 3.6× |
| SVD | 256 | 16,459 | 8,423 | 2.0× |
| SVD | 512 | 100,127 | 40,355 | 2.5× |
| SVD | 1,024 | 625,263 | 208,055 | 3.0× |

### Linear Solve (`A\b`) — time in μs (lower is better)

| Type | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |
|------|---|---------------|-----------------|---------|
| Float64 | 256 | 647 | 323 | 2.0× |
| Float64 | 1,024 | 24,789 | 7,310 | 3.4× |
| Float64 | 2,048 | 175,949 | 50,714 | 3.5× |
| Float32 | 256 | 448 | 176 | 2.5× |
| Float32 | 1,024 | 14,337 | 3,626 | 4.0× |
| Float32 | 2,048 | 95,772 | 21,360 | 4.5× |

!!! note "Benchmark environment"
    Apple M2 Max, macOS 26, single-threaded. Julia 1.12.5, AppleAccelerate v0.6.0, LinearAlgebra v1.12.0 (OpenBLAS 0.3.29). Times are minimum of 5 trials. OpenBLAS benchmarked before loading AppleAccelerate; Accelerate benchmarked after `using AppleAccelerate` forwards BLAS via LBT. Run [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl) to reproduce.
