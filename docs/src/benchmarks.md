# Benchmarks

All benchmarks were run on an Apple M2 Max, macOS 26, single-threaded, with Julia 1.12.6 and AppleAccelerate.jl v0.8.0 (all tables measured July 2026). Times are the minimum of 5 trials.

Single-threaded execution is used to compare **kernel quality** on a level footing. This is *not* representative of how OpenBLAS is normally deployed — OpenBLAS scales GEMM across CPU cores, whereas Accelerate offloads large GEMM to a shared SME co-processor that a single thread already saturates. See the GEMM note below for a multi-threaded comparison.

## How these benchmarks are generated

Every table on this page is produced by the scripts in
[`test/bench/`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/tree/master/test/bench),
driven by
[`run_benchmarks.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/run_benchmarks.jl):

```
julia --project=test/bench test/bench/run_benchmarks.jl          # all suites
julia --project=test/bench test/bench/run_benchmarks.jl dense    # one suite
```

Each suite runs in its own fresh Julia process so that the OpenBLAS/SuiteSparse
baselines are measured *before* AppleAccelerate is loaded (loading it forwards
BLAS/LAPACK through libblastrampoline for the rest of the process). All suites
default to a single thread; the dense suite additionally honors a
`BENCH_THREADS` environment variable for thread-scaling sweeps like the GEMM
table below:

```
BENCH_THREADS=4 julia --project=test/bench test/bench/run_benchmarks.jl dense
```

## Array Operations

Performance comparison of vDSP array operations vs Julia Base equivalents (`map(Base.f, X)` for unary, `@simd` loops for binary/compound). Source: [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl).

### Unary Math Functions

Transcendental functions show the biggest gains — vDSP is 7–19× faster for `sin`/`cos` on Float32:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| exp | Float64 | 100,000 | 176 | 356 | 2.0× |
| log | Float64 | 100,000 | 174 | 457 | 2.6× |
| sin | Float64 | 100,000 | 106 | 740 | 7.0× |
| cos | Float64 | 100,000 | 113 | 763 | 6.8× |
| sqrt | Float64 | 100,000 | 30 | 59 | 2.0× |
| exp | Float32 | 100,000 | 42 | 330 | 8.0× |
| log | Float32 | 100,000 | 56 | 381 | 6.7× |
| sin | Float32 | 100,000 | 38 | 707 | 18.8× |
| cos | Float32 | 100,000 | 38 | 717 | 18.7× |
| sqrt | Float32 | 100,000 | 28 | 69 | 2.5× |

### Reductions

`sum`/`maximum`/`minimum` are 1.1–2× faster. Note: `dot` is slower via vDSP because Julia's `LinearAlgebra.dot` already uses the Accelerate-forwarded BLAS:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| sum | Float64 | 1,000,000 | 90 | 100 | 1.1× |
| maximum | Float64 | 1,000,000 | 90 | 94 | 1.1× |
| minimum | Float64 | 1,000,000 | 88 | 93 | 1.1× |
| sum | Float32 | 1,000,000 | 44 | 55 | 1.2× |
| maximum | Float32 | 1,000,000 | 45 | 50 | 1.1× |
| minimum | Float32 | 1,000,000 | 45 | 50 | 1.1× |

### Binary Element-wise Ops

Addition and multiplication are memory-bandwidth-bound, so results hover near parity — Julia's `@simd` loop even edges out vDSP for Float64 `vadd`, while vDSP is modestly faster elsewhere:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| vadd | Float64 | 1,000,000 | 349 | 253 | 1.4× slower |
| vmul | Float64 | 1,000,000 | 239 | 280 | 1.2× |
| vadd | Float32 | 1,000,000 | 111 | 120 | 1.1× |
| vmul | Float32 | 1,000,000 | 111 | 117 | 1.1× |

!!! note "Benchmark environment"
    Julia reference uses `map(Base.f, X)` for unary ops and `@inbounds @simd` loops for binary/compound ops. Source: [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl).

## Dense Linear Algebra

Performance comparison of Apple Accelerate vs OpenBLAS. The dense benchmark script loads OpenBLAS first, then switches to Accelerate via LBT, so both are measured in the same process. Source: [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl).

### GEMM (`mul!`) — GFLOPS (higher is better)

Single-threaded, Accelerate is 6–13× faster for matrix multiply, with the largest gains for Float32:

| Type | N | OpenBLAS | Accelerate | Speedup |
|------|---|----------|------------|---------|
| Float64 | 64 | 24 | 225 | 9.3× |
| Float64 | 256 | 46 | 350 | 7.6× |
| Float64 | 1,024 | 52 | 359 | 6.9× |
| Float64 | 4,096 | 52 | 300 | 5.7× |
| Float32 | 64 | 88 | 547 | 6.2× |
| Float32 | 256 | 100 | 1,207 | 12.1× |
| Float32 | 1,024 | 105 | 1,394 | 13.3× |
| Float32 | 4,096 | 105 | 1,127 | 10.8× |

!!! warning "These are single-threaded numbers — read this before comparing"
    The table above pins **both** libraries to a single thread (`BLAS.set_num_threads(1)`),
    which is *not* how OpenBLAS is normally used and understates it substantially.

    The two libraries use different execution units on Apple Silicon:

    - **OpenBLAS** runs GEMM on the CPU's P-cores (NEON), and **scales with thread count**.
    - **Accelerate** dispatches large GEMM to the **SME/AMX matrix co-processor** — a
      dedicated unit shared *per core cluster*. This is a CPU-side co-processor, **not** the
      GPU (Apple GPUs do not natively support Float64, and
      [Apple documents Accelerate as CPU-based](https://developer.apple.com/documentation/accelerate#overview)).
      Because each unit is shared and a single thread on
      that cluster already saturates it, Accelerate GEMM throughput scales only with the
      number of co-processors (i.e. the number of core clusters), **not** with
      `BLAS.set_num_threads` beyond that.

    Measured on this page's reference machine, an **Apple M2 Max** (8 P-cores in two
    clusters + 4 E-cores), Float64/Float32 GEMM at N=4096:

    | Threads | OpenBLAS F64 | OB speedup | Accelerate F64 | OpenBLAS F32 | Accelerate F32 |
    |---------|-------------:|:----------:|---------------:|-------------:|---------------:|
    | 1 |  52.5 | 1.0× | 314.8 | 105.1 | 1158.3 |
    | 2 | 104.4 | 2.0× | 616.2 | 206.7 | 2251.8 |
    | 3 | 152.1 | 2.9× | 615.8 | 299.4 | 2236.7 |
    | 4 | 205.8 | 3.9× | 617.4 | 398.8 | 2254.9 |
    | 6 | 276.2 | 5.3× | 611.7 | 561.1 | 2248.7 |
    | 8 | 354.1 | 6.7× | 618.1 | 715.6 | 2240.3 |

    (GFLOPS; higher is better. This is a separate thread-scaling sweep, so its 1-thread
    column differs slightly from the single-threaded table above due to run-to-run variance.)
    Two effects are visible. Accelerate **doubles** from 1→2
    threads — the M2 Max has *two* P-core clusters, so a second thread engages the second
    SME unit — then stays flat, since there are no more co-processors to recruit. OpenBLAS
    keeps scaling across all 8 P-cores (~6.7×). On a single-P-cluster chip such as the M4,
    Accelerate is flat from thread 1 (see the M4 data in
    [issue #132](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/issues/132)).

    So the multiple by which Accelerate leads shrinks as OpenBLAS is given more cores.
    Accelerate still wins here, but the single-threaded speedups above should be read as a
    kernel-quality comparison, not as the gap you will see against a fully-threaded OpenBLAS.
    See Apple's [CPU optimization guide](https://developer.apple.com/documentation/apple-silicon/cpu-optimization-guide)
    for details on the SME engine, and [arXiv:2409.18779](https://arxiv.org/abs/2409.18779)
    for independent benchmarks of matrix-multiply throughput on Apple Silicon's
    matrix co-processor.

### Factorizations (Float64) — time in μs (lower is better)

| Operation | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |
|-----------|---|---------------|-----------------|---------|
| LU | 256 | 440 | 217 | 2.0× |
| LU | 1,024 | 17,115 | 11,274 | 1.5× |
| LU | 2,048 | 123,414 | 68,209 | 1.8× |
| QR | 512 | 5,686 | 3,428 | 1.7× |
| QR | 2,048 | 281,766 | 145,400 | 1.9× |
| Cholesky | 256 | 245 | 86 | 2.8× |
| Cholesky | 1,024 | 9,333 | 2,835 | 3.3× |
| SVD | 256 | 11,489 | 6,590 | 1.7× |
| SVD | 512 | 69,465 | 31,085 | 2.2× |
| SVD | 1,024 | 438,359 | 161,072 | 2.7× |

### Linear Solve (`A\b`) — time in μs (lower is better)

| Type | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |
|------|---|---------------|-----------------|---------|
| Float64 | 256 | 444 | 201 | 2.2× |
| Float64 | 1,024 | 16,791 | 5,191 | 3.2× |
| Float64 | 2,048 | 121,863 | 37,958 | 3.2× |
| Float32 | 256 | 303 | 123 | 2.5× |
| Float32 | 1,024 | 9,859 | 2,578 | 3.8× |
| Float32 | 2,048 | 66,486 | 14,703 | 4.5× |

!!! note "Benchmark environment"
    LinearAlgebra.jl v1.12.0 (OpenBLAS 0.3.29). OpenBLAS benchmarked before loading AppleAccelerate; Accelerate benchmarked after `using AppleAccelerate` forwards BLAS via LBT. Source: [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl).

## Sparse Linear Algebra

Performance comparison of Apple Sparse Solvers vs SuiteSparse (CHOLMOD/UMFPACK). The sparse benchmark script runs SuiteSparse first (before loading AppleAccelerate), so SuiteSparse uses OpenBLAS internally, then loads AppleAccelerate and re-runs with Apple's sparse solvers. Source: [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl).

### Sparse Matrix-Vector Multiply (density=0.01)

SuiteSparse CSC SpMV is 2.4–4.6× faster due to its simpler data layout:

| Type | N | Apple (μs) | SuiteSparse (μs) | Ratio |
|------|---|-----------|-------------------|-------|
| Float64 | 1,000 | 23 | 8 | 0.33× |
| Float64 | 10,000 | 2,103 | 506 | 0.24× |
| Float64 | 50,000 | 52,442 | 21,756 | 0.41× |
| Float32 | 1,000 | 24 | 7 | 0.30× |
| Float32 | 10,000 | 2,081 | 451 | 0.22× |

### QR Factorize + Solve

SuiteSparse LU (`\`) is faster for Float64. Apple QR wins for Float32 at N≥1000 (up to 2.3×):

| Type | N | Apple (μs) | SuiteSparse (μs) | Speedup |
|------|---|-----------|-------------------|---------|
| Float64 | 500 | 2,471 | 1,937 | 0.78× |
| Float64 | 2,000 | 120,362 | 102,256 | 0.85× |
| Float64 | 5,000 | 1,996,559 | 1,804,746 | 0.90× |
| Float32 | 1,000 | 9,397 | 16,314 | 1.74× |
| Float32 | 2,000 | 54,188 | 122,508 | 2.26× |
| Float32 | 5,000 | 960,933 | 1,833,739 | 1.91× |

### Cholesky Factorize + Solve

Apple Cholesky is around parity at N=500 and pulls ahead as N grows, reaching ~4× at N=5000:

| Type | N | Apple (μs) | SuiteSparse (μs) | Speedup |
|------|---|-----------|-------------------|---------|
| Float64 | 500 | 1,453 | 1,521 | 1.05× |
| Float64 | 2,000 | 45,324 | 77,211 | 1.70× |
| Float64 | 5,000 | 257,123 | 1,059,278 | 4.12× |
| Float32 | 2,000 | 30,099 | 51,089 | 1.70× |
| Float32 | 5,000 | 187,216 | 605,418 | 3.23× |

!!! note "Benchmark environment"
    SparseArrays.jl v1.12.0 (SuiteSparse 7.8.3). Matrices have density 0.01. Source: [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl).

!!! warning "Benchmark limitations"
    These benchmarks use random sparse matrices (`sprandn`) which lack the structure found in real-world problems (e.g., banded, block-diagonal, or mesh-derived sparsity patterns). The matrix sizes tested (N up to 5,000–50,000) are also modest by sparse solver standards. Performance on structured problems from applications like FEM, circuit simulation, or graph analysis may differ significantly.

## FFT

Pre-planned FFT performance comparing Apple vDSP against FFTW (`FFTW.MEASURE` plans), both single-threaded. On this machine FFTW leads across sizes and dimensions; vDSP only reaches parity for Float32 at the largest 1D sizes. vDSP's FFT is still useful when you want to avoid an FFTW dependency or already have Accelerate loaded. Source: [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl).

### Complex 1D FFT

FFTW leads for 1D complex transforms, typically by 1.4–2.5×; vDSP only reaches parity for Float32 at the largest sizes:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| ComplexF64 | 1,024 | 5.5 | 2.5 | 2.16× slower |
| ComplexF64 | 4,096 | 19.7 | 9.4 | 2.1× slower |
| ComplexF64 | 65,536 | 603 | 244 | 2.47× slower |
| ComplexF64 | 1,048,576 | 12,965 | 9,309 | 1.39× slower |
| ComplexF32 | 1,024 | 2.1 | 1.2 | 1.82× slower |
| ComplexF32 | 4,096 | 7.7 | 5.0 | 1.55× slower |
| ComplexF32 | 65,536 | 149 | 159 | 1.06× |
| ComplexF32 | 1,048,576 | 7,218 | 5,905 | 1.22× slower |

### Real FFT

FFTW leads for real transforms too, though vDSP closes to within ~1.1–1.3× for Float32 at large sizes:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| Float64 | 1,024 | 2.5 | 1.0 | 2.44× slower |
| Float64 | 65,536 | 244 | 124 | 1.96× slower |
| Float32 | 4,096 | 7.9 | 3.8 | 2.07× slower |
| Float32 | 65,536 | 93 | 71 | 1.31× slower |
| Float32 | 262,144 | 386 | 340 | 1.14× slower |

### Complex 2D FFT

FFTW is substantially faster for 2D complex transforms — up to ~14× for large Float64 grids — so prefer FFTW when 2D FFTs dominate:

| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |
|------|------|-----------|-----------|---------|
| ComplexF64 | 64×64 | 33.0 | 8.0 | 4.14× slower |
| ComplexF64 | 256×256 | 3,204 | 234 | 13.69× slower |
| ComplexF32 | 64×64 | 13.2 | 4.6 | 2.86× slower |
| ComplexF32 | 256×256 | 228 | 124 | 1.84× slower |

!!! note "Benchmark environment"
    FFTW.jl v1.10.0. Both use pre-planned transforms. Source: [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl).
