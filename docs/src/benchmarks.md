# Benchmarks

All benchmarks were run on an Apple M2 Max, macOS 26, single-threaded, with Julia 1.12.5 and AppleAccelerate.jl v0.6.0. Times are the minimum of 5 trials.

Run the full benchmark suite with:

```
julia --project=test/bench test/bench/run_benchmarks.jl
```

## Array Operations

Performance comparison of vDSP array operations vs Julia Base equivalents (`map(Base.f, X)` for unary, `@simd` loops for binary/compound). Source: [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl).

### Unary Math Functions

Transcendental functions show the biggest gains — vDSP is 7–19× faster for `sin`/`cos` on Float32:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| exp | Float64 | 100,000 | 171 | 440 | 2.6× |
| log | Float64 | 100,000 | 191 | 652 | 3.4× |
| sin | Float64 | 100,000 | 152 | 1,051 | 6.9× |
| cos | Float64 | 100,000 | 160 | 1,083 | 6.8× |
| sqrt | Float64 | 100,000 | 42 | 97 | 2.3× |
| exp | Float32 | 100,000 | 60 | 464 | 7.8× |
| log | Float32 | 100,000 | 80 | 542 | 6.8× |
| sin | Float32 | 100,000 | 54 | 1,012 | 18.8× |
| cos | Float32 | 100,000 | 55 | 1,035 | 18.9× |
| sqrt | Float32 | 100,000 | 40 | 98 | 2.5× |

### Reductions

`sum`/`maximum`/`minimum` are 1.1–2× faster. Note: `dot` is slower via vDSP because Julia's `LinearAlgebra.dot` already uses the Accelerate-forwarded BLAS:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| sum | Float64 | 1,000,000 | 128 | 143 | 1.1× |
| maximum | Float64 | 1,000,000 | 128 | 136 | 1.1× |
| minimum | Float64 | 1,000,000 | 127 | 134 | 1.1× |
| sum | Float32 | 1,000,000 | 62 | 77 | 1.2× |
| maximum | Float32 | 1,000,000 | 62 | 71 | 1.1× |
| minimum | Float32 | 1,000,000 | 63 | 71 | 1.1× |

### Binary Element-wise Ops

Addition and multiplication are memory-bandwidth-bound for Float64; vDSP is faster for Float32 at large sizes:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| vadd | Float64 | 1,000,000 | 354 | 360 | 1.0× |
| vmul | Float64 | 1,000,000 | 343 | 340 | 1.0× |
| vadd | Float32 | 1,000,000 | 133 | 167 | 1.3× |
| vmul | Float32 | 1,000,000 | 87 | 167 | 1.9× |

!!! note "Benchmark environment"
    Julia reference uses `map(Base.f, X)` for unary ops and `@inbounds @simd` loops for binary/compound ops. Source: [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl).

## Dense Linear Algebra

Performance comparison of Apple Accelerate vs OpenBLAS. The dense benchmark script loads OpenBLAS first, then switches to Accelerate via LBT, so both are measured in the same process. Source: [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl).

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
    LinearAlgebra.jl v1.12.0 (OpenBLAS 0.3.29). OpenBLAS benchmarked before loading AppleAccelerate; Accelerate benchmarked after `using AppleAccelerate` forwards BLAS via LBT. Source: [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl).

## Sparse Linear Algebra

Performance comparison of Apple Sparse Solvers vs SuiteSparse (CHOLMOD/UMFPACK). The sparse benchmark script runs SuiteSparse first (before loading AppleAccelerate), so SuiteSparse uses OpenBLAS internally, then loads AppleAccelerate and re-runs with Apple's sparse solvers. Source: [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl).

### Sparse Matrix-Vector Multiply (density=0.01)

SuiteSparse CSC SpMV is 2.5–4.6× faster due to its simpler data layout:

| Type | N | Apple (μs) | SuiteSparse (μs) | Ratio |
|------|---|-----------|-------------------|-------|
| Float64 | 1,000 | 35 | 10 | 0.30× |
| Float64 | 10,000 | 2,998 | 728 | 0.24× |
| Float64 | 50,000 | 76,061 | 30,615 | 0.40× |
| Float32 | 1,000 | 35 | 10 | 0.29× |
| Float32 | 10,000 | 2,938 | 644 | 0.22× |

### QR Factorize + Solve

SuiteSparse LU (`\`) is faster for Float64. Apple QR wins for Float32 at N≥500 (up to 1.8×):

| Type | N | Apple (μs) | SuiteSparse (μs) | Speedup |
|------|---|-----------|-------------------|---------|
| Float64 | 500 | 3,057 | 2,610 | 0.85× |
| Float64 | 2,000 | 148,372 | 98,608 | 0.66× |
| Float64 | 5,000 | 2,481,285 | 1,834,568 | 0.74× |
| Float32 | 1,000 | 11,798 | 16,460 | 1.40× |
| Float32 | 2,000 | 70,343 | 124,074 | 1.76× |
| Float32 | 5,000 | 1,307,779 | 1,800,863 | 1.38× |

### Cholesky Factorize + Solve

Apple Cholesky is faster at N=5000, SuiteSparse faster at smaller sizes:

| Type | N | Apple (μs) | SuiteSparse (μs) | Speedup |
|------|---|-----------|-------------------|---------|
| Float64 | 500 | 1,957 | 1,546 | 0.79× |
| Float64 | 2,000 | 52,412 | 46,820 | 0.89× |
| Float64 | 5,000 | 344,766 | 505,765 | 1.47× |
| Float32 | 2,000 | 42,128 | 33,195 | 0.79× |
| Float32 | 5,000 | 273,444 | 339,105 | 1.24× |

!!! note "Benchmark environment"
    SparseArrays.jl v1.12.0 (SuiteSparse 7.8.3). Matrices have density 0.01. Source: [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl).

!!! warning "Benchmark limitations"
    These benchmarks use random sparse matrices (`sprandn`) which lack the structure found in real-world problems (e.g., banded, block-diagonal, or mesh-derived sparsity patterns). The matrix sizes tested (N up to 5,000–50,000) are also modest by sparse solver standards. Performance on structured problems from applications like FEM, circuit simulation, or graph analysis may differ significantly.

## FFT

Pre-planned FFT performance comparing Apple vDSP against FFTW, both single-threaded. Source: [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl).

### Complex 1D FFT

vDSP and FFTW are closely matched. vDSP is notably faster at N=4096; at larger sizes FFTW has a slight edge:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| ComplexF64 | 1,024 | 7.0 | 7.1 | 1.01× |
| ComplexF64 | 4,096 | 20.5 | 34.1 | 1.67× |
| ComplexF64 | 65,536 | 718 | 653 | 0.91× |
| ComplexF64 | 1,048,576 | 18,767 | 17,208 | 0.92× |
| ComplexF32 | 1,024 | 3.1 | 3.1 | 1.01× |
| ComplexF32 | 4,096 | 11.5 | 18.4 | 1.61× |
| ComplexF32 | 65,536 | 329 | 323 | 0.98× |
| ComplexF32 | 1,048,576 | 4,891 | 4,753 | 0.97× |

### Real FFT

For Float32, vDSP is faster at larger sizes (up to 2.2×). Float64 rfft favors FFTW:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| Float64 | 1,024 | 3.8 | 1.6 | 0.42× |
| Float64 | 65,536 | 360 | 209 | 0.58× |
| Float32 | 4,096 | 7.0 | 7.6 | 1.09× |
| Float32 | 65,536 | 115 | 169 | 1.47× |
| Float32 | 262,144 | 539 | 1,166 | 2.16× |

### Complex 2D FFT

Nearly identical performance:

| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |
|------|------|-----------|-----------|---------|
| ComplexF64 | 64×64 | 47 | 47 | 1.00× |
| ComplexF64 | 256×256 | 4,659 | 4,708 | 1.01× |
| ComplexF32 | 64×64 | 18 | 22 | 1.19× |
| ComplexF32 | 256×256 | 321 | 325 | 1.01× |

!!! note "Benchmark environment"
    FFTW.jl v1.10.0. Both use pre-planned transforms. Source: [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl).
