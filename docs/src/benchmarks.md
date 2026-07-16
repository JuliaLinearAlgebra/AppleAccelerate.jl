# Benchmarks

All benchmarks were run on an Apple M2 Max, macOS 26, single-threaded, with Julia 1.12.6 and AppleAccelerate.jl v0.8.0 master (all tables measured July 2026, after the 2D real FFT / batched FFT / FFT setup-caching work). Times are the minimum of 5 trials.

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
BLAS/LAPACK through libblastrampoline for the rest of the process). The tables
on this page were produced by running the suites **one at a time** so no two
benchmarks compete for the machine. All suites default to a single thread; the
dense suite additionally honors a `BENCH_THREADS` environment variable for
thread-scaling sweeps like the GEMM table below:

```
BENCH_THREADS=4 julia --project=test/bench test/bench/run_benchmarks.jl dense
```

## Array Operations

Performance comparison of vDSP array operations vs Julia Base equivalents (`map(Base.f, X)` for unary, `@simd` loops for binary/compound). Source: [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl).

### Unary Math Functions

Transcendental functions show the biggest gains — vDSP is 7–19× faster for `sin`/`cos` on Float32:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| exp | Float64 | 100,000 | 135 | 301 | 2.2× |
| log | Float64 | 100,000 | 174 | 500 | 2.9× |
| sin | Float64 | 100,000 | 147 | 731 | 5.0× |
| cos | Float64 | 100,000 | 113 | 759 | 6.7× |
| sqrt | Float64 | 100,000 | 30 | 59 | 2.0× |
| exp | Float32 | 100,000 | 42 | 319 | 7.7× |
| log | Float32 | 100,000 | 56 | 381 | 6.8× |
| sin | Float32 | 100,000 | 38 | 704 | 18.8× |
| cos | Float32 | 100,000 | 38 | 715 | 18.7× |
| sqrt | Float32 | 100,000 | 28 | 69 | 2.5× |

### Reductions

`sum`/`maximum`/`minimum` are 1.0–1.2× faster. Note: `dot` is slower via vDSP because Julia's `LinearAlgebra.dot` already uses the Accelerate-forwarded BLAS:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| sum | Float64 | 1,000,000 | 90 | 101 | 1.1× |
| maximum | Float64 | 1,000,000 | 90 | 97 | 1.1× |
| minimum | Float64 | 1,000,000 | 93 | 94 | 1.0× |
| sum | Float32 | 1,000,000 | 44 | 54 | 1.2× |
| maximum | Float32 | 1,000,000 | 45 | 49 | 1.1× |
| minimum | Float32 | 1,000,000 | 44 | 49 | 1.1× |

### Binary Element-wise Ops

Addition and multiplication are memory-bandwidth-bound, so Float64 results hover at parity, while vDSP shows a modest edge for Float32:

| Op | Type | N | vDSP (μs) | Julia (μs) | Speedup |
|----|------|---|-----------|------------|---------|
| vadd | Float64 | 1,000,000 | 254 | 252 | 1.0× |
| vmul | Float64 | 1,000,000 | 255 | 240 | 1.1× slower |
| vadd | Float32 | 1,000,000 | 89 | 141 | 1.6× |
| vmul | Float32 | 1,000,000 | 135 | 140 | 1.0× |

!!! note "Benchmark environment"
    Julia reference uses `map(Base.f, X)` for unary ops and `@inbounds @simd` loops for binary/compound ops. Source: [`bench_array.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_array.jl).

## Dense Linear Algebra

Performance comparison of Apple Accelerate vs OpenBLAS. The dense benchmark script loads OpenBLAS first, then switches to Accelerate via LBT, so both are measured in the same process. Source: [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl).

### GEMM (`mul!`) — GFLOPS (higher is better)

Single-threaded, Accelerate is 6–13× faster for matrix multiply, with the largest gains for Float32:

| Type | N | OpenBLAS | Accelerate | Speedup |
|------|---|----------|------------|---------|
| Float64 | 64 | 24 | 237 | 9.8× |
| Float64 | 256 | 45 | 349 | 7.7× |
| Float64 | 1,024 | 51 | 360 | 7.1× |
| Float64 | 4,096 | 52 | 291 | 5.6× |
| Float32 | 64 | 88 | 599 | 6.8× |
| Float32 | 256 | 100 | 1,217 | 12.2× |
| Float32 | 1,024 | 104 | 1,396 | 13.5× |
| Float32 | 4,096 | 106 | 1,136 | 10.8× |

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
| LU | 256 | 429 | 220 | 2.0× |
| LU | 1,024 | 16,956 | 11,190 | 1.5× |
| LU | 2,048 | 123,071 | 68,307 | 1.8× |
| QR | 512 | 5,571 | 3,422 | 1.6× |
| QR | 2,048 | 280,167 | 145,190 | 1.9× |
| Cholesky | 256 | 252 | 88 | 2.9× |
| Cholesky | 1,024 | 9,313 | 2,688 | 3.5× |
| SVD | 256 | 11,533 | 6,609 | 1.7× |
| SVD | 512 | 68,990 | 30,973 | 2.2× |
| SVD | 1,024 | 435,259 | 157,894 | 2.8× |

### Linear Solve (`A\b`) — time in μs (lower is better)

| Type | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |
|------|---|---------------|-----------------|---------|
| Float64 | 256 | 458 | 193 | 2.4× |
| Float64 | 1,024 | 17,156 | 5,341 | 3.2× |
| Float64 | 2,048 | 123,081 | 37,501 | 3.3× |
| Float32 | 256 | 312 | 124 | 2.5× |
| Float32 | 1,024 | 9,799 | 2,576 | 3.8× |
| Float32 | 2,048 | 65,855 | 14,773 | 4.5× |

!!! note "Benchmark environment"
    LinearAlgebra.jl v1.12.0 (OpenBLAS 0.3.29). OpenBLAS benchmarked before loading AppleAccelerate; Accelerate benchmarked after `using AppleAccelerate` forwards BLAS via LBT. Source: [`bench_dense.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_dense.jl).

## Sparse Linear Algebra

Performance comparison of Apple Sparse Solvers vs SuiteSparse (CHOLMOD/UMFPACK). The sparse benchmark script runs SuiteSparse first (before loading AppleAccelerate), so SuiteSparse uses OpenBLAS internally, then loads AppleAccelerate and re-runs with Apple's sparse solvers. Source: [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl).

### Sparse Matrix-Vector Multiply (density=0.01)

SuiteSparse CSC SpMV is 2.4–4.5× faster due to its simpler data layout:

| Type | N | Apple (μs) | SuiteSparse (μs) | Ratio |
|------|---|-----------|-------------------|-------|
| Float64 | 1,000 | 23 | 7 | 0.32× |
| Float64 | 10,000 | 2,139 | 512 | 0.24× |
| Float64 | 50,000 | 52,121 | 22,038 | 0.42× |
| Float32 | 1,000 | 23 | 7 | 0.31× |
| Float32 | 10,000 | 2,106 | 472 | 0.22× |

### QR Factorize + Solve

SuiteSparse LU (`\`) is modestly faster for Float64 (near parity by N=2000). Apple QR wins for Float32 at N≥1000 (up to 2×):

| Type | N | Apple (μs) | SuiteSparse (μs) | Speedup |
|------|---|-----------|-------------------|---------|
| Float64 | 500 | 2,412 | 1,990 | 0.83× |
| Float64 | 2,000 | 117,360 | 119,056 | 1.01× |
| Float64 | 5,000 | 2,038,889 | 1,845,175 | 0.91× |
| Float32 | 1,000 | 9,521 | 13,286 | 1.40× |
| Float32 | 2,000 | 57,290 | 116,995 | 2.04× |
| Float32 | 5,000 | 966,020 | 1,820,923 | 1.88× |

### Cholesky Factorize + Solve

Apple Cholesky is around parity at N=500 and pulls ahead as N grows, reaching ~4× at N=5000:

| Type | N | Apple (μs) | SuiteSparse (μs) | Speedup |
|------|---|-----------|-------------------|---------|
| Float64 | 500 | 1,488 | 1,513 | 1.02× |
| Float64 | 2,000 | 38,322 | 78,226 | 2.04× |
| Float64 | 5,000 | 258,357 | 1,042,905 | 4.04× |
| Float32 | 2,000 | 30,550 | 51,237 | 1.68× |
| Float32 | 5,000 | 185,574 | 613,517 | 3.31× |

!!! note "Benchmark environment"
    SparseArrays.jl v1.12.0 (SuiteSparse 7.8.3). Matrices have density 0.01. Source: [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl).

!!! warning "Benchmark limitations"
    These benchmarks use random sparse matrices (`sprandn`) which lack the structure found in real-world problems (e.g., banded, block-diagonal, or mesh-derived sparsity patterns). The matrix sizes tested (N up to 5,000–50,000) are also modest by sparse solver standards. Performance on structured problems from applications like FEM, circuit simulation, or graph analysis may differ significantly.

## FFT

FFT performance comparing Apple vDSP against FFTW, both single-threaded. With **pre-planned** transforms (`FFTW.MEASURE` plans) FFTW leads across sizes and dimensions — typically 1.5–2.4× for 1D and up to ~14× for large 2D Float64 grids — with vDSP reaching parity only for Float32 at the largest 1D sizes. For **no-plan convenience calls**, however, the picture changed with the FFT setup cache: `AppleAccelerate.fft(x)` reuses a cached setup after the first call at each size, and now beats FFTW's unplanned `fft(x)` at small sizes and wins or ties for Float32 at every size tested (see the last table). vDSP's FFT also remains useful to avoid an FFTW dependency or when Accelerate is already loaded. Source: [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl).

### Complex 1D FFT (pre-planned)

FFTW leads for 1D complex transforms, typically by 1.6–2.4×; vDSP reaches parity for Float32 at some large sizes:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| ComplexF64 | 1,024 | 4.6 | 2.0 | 2.36× slower |
| ComplexF64 | 4,096 | 19.5 | 9.7 | 2.0× slower |
| ComplexF64 | 65,536 | 611 | 263 | 2.32× slower |
| ComplexF64 | 1,048,576 | 13,423 | 8,249 | 1.63× slower |
| ComplexF32 | 1,024 | 2.2 | 1.2 | 1.79× slower |
| ComplexF32 | 4,096 | 10.9 | 4.6 | 2.36× slower |
| ComplexF32 | 65,536 | 185 | 164 | 1.12× slower |
| ComplexF32 | 262,144 | 653 | 675 | 1.03× |
| ComplexF32 | 1,048,576 | 7,143 | 5,850 | 1.22× slower |

### Real FFT (pre-planned)

FFTW leads for real transforms too, though vDSP closes to near parity for Float32 at large sizes:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| Float64 | 1,024 | 3.3 | 1.0 | 3.16× slower |
| Float64 | 65,536 | 241 | 128 | 1.87× slower |
| Float32 | 4,096 | 6.2 | 4.1 | 1.51× slower |
| Float32 | 65,536 | 93 | 71 | 1.3× slower |
| Float32 | 262,144 | 369 | 355 | 1.04× slower |

### Complex 2D FFT (pre-planned)

FFTW is substantially faster for 2D complex transforms — up to ~14× for large Float64 grids — so prefer FFTW when 2D FFTs dominate:

| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |
|------|------|-----------|-----------|---------|
| ComplexF64 | 64×64 | 32.5 | 7.9 | 4.12× slower |
| ComplexF64 | 256×256 | 3,255 | 231 | 14.06× slower |
| ComplexF32 | 64×64 | 13.6 | 4.7 | 2.92× slower |
| ComplexF32 | 256×256 | 228 | 124 | 1.84× slower |

### Real 2D FFT (pre-planned)

The new 2D real FFT (`rfft(A::Matrix)`) follows the same pattern — FFTW leads by roughly 2–3.5× at moderate sizes, more for large Float64 grids:

| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |
|------|------|-----------|-----------|---------|
| Float64 | 64×64 | 17.3 | 5.0 | 3.5× slower |
| Float64 | 128×128 | 69.4 | 23.2 | 2.99× slower |
| Float64 | 256×256 | 1,181 | 110 | 10.74× slower |
| Float32 | 64×64 | 7.0 | 3.8 | 1.86× slower |
| Float32 | 128×128 | 30.4 | 15.6 | 1.95× slower |
| Float32 | 256×256 | 156 | 69 | 2.25× slower |

### Batched 1D FFT (columns of a matrix)

The new `fft(A, 1)` transforms each column via `vDSP_fftm_*`. FFTW's planned batched transforms lead by 2–4.4×:

| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |
|------|------|-----------|-----------|---------|
| ComplexF64 | 256×64 | 57.4 | 23.0 | 2.5× slower |
| ComplexF64 | 1024×64 | 257 | 124 | 2.07× slower |
| ComplexF64 | 4096×16 | 294 | 152 | 1.94× slower |
| ComplexF32 | 256×64 | 50.6 | 11.4 | 4.43× slower |
| ComplexF32 | 1024×64 | 165 | 59 | 2.77× slower |
| ComplexF32 | 4096×16 | 181 | 80 | 2.25× slower |

### No-plan `fft(x)` — cached setups

This table compares the *convenience* entry points: `AppleAccelerate.fft(x)` vs `FFTW.fft(x)`, neither given a pre-built plan. AppleAccelerate reuses a cached vDSP setup after the first call at each size (FFTW builds an `ESTIMATE` plan per call), so for the common "just call `fft`" pattern vDSP now leads at small sizes and wins or ties for Float32 at every size tested:

| Type | N | vDSP (μs) | FFTW (μs) | Speedup |
|------|---|-----------|-----------|---------|
| ComplexF64 | 1,024 | 3.8 | 8.3 | 2.16× |
| ComplexF64 | 65,536 | 465 | 263 | 1.77× slower |
| ComplexF64 | 1,048,576 | 13,531 | 8,883 | 1.52× slower |
| ComplexF32 | 1,024 | 2.3 | 4.2 | 1.79× |
| ComplexF32 | 65,536 | 145 | 146 | 1.01× |
| ComplexF32 | 1,048,576 | 4,704 | 5,403 | 1.15× |

!!! note "Benchmark environment"
    FFTW.jl v1.10.0. Pre-planned tables use `FFTW.MEASURE` plans on both sides; the no-plan table calls the convenience functions directly. Source: [`bench_fft.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_fft.jl).
