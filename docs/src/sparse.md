# Sparse Linear Algebra

AppleAccelerate wraps Apple's [Sparse Solvers](https://developer.apple.com/documentation/accelerate/sparse_solvers) library for sparse matrix operations and direct solvers.

```@setup sparse
using AppleAccelerate, SparseArrays, LinearAlgebra
```

## AASparseMatrix

A wrapper around Apple's sparse matrix format, constructed from Julia's `SparseMatrixCSC`.

```@example sparse
using AppleAccelerate, SparseArrays
import AppleAccelerate: AASparseMatrix, muladd!

A_jl = sprandn(100, 100, 0.05)
A = AASparseMatrix(A_jl)

x = randn(100)
y = A * x  # Sparse matrix-vector multiply
nothing # hide
```

The constructor automatically detects symmetric and triangular structure and sets
the appropriate Apple Accelerate attributes.

### Matrix operations

| Function | Description |
|----------|-------------|
| [`AASparseMatrix`](@ref AppleAccelerate.AASparseMatrix) | Construct from Julia sparse matrix |
| `A * x` | Sparse matrix-vector or matrix-matrix multiply |
| `alpha * A * x` | Scaled sparse multiply |
| [`muladd!`](@ref AppleAccelerate.muladd!) | Multiply-add: `y += A * x` or `y += alpha * A * x` |
| `transpose(A)` | Transpose (sets flag, no copy) |

### Query functions

| Function | Description |
|----------|-------------|
| `size(A)` | Matrix dimensions |
| `eltype(A)` | Element type |
| `issymmetric(A)` | Check if symmetric |
| `istriu(A)` | Check if upper triangular |
| `istril(A)` | Check if lower triangular |
| `A[i, j]` | Element access |

## AAFactorization

Lazy factorization wrapper: the factorization is computed on the first call to `solve`
or by explicitly calling `factor!`.

```@example sparse
import AppleAccelerate: AAFactorization, solve, solve!, factor!

A = sprandn(100, 100, 0.1) + 20I
f = AAFactorization(A)

# Factorization computed lazily on first solve
b = randn(100)
x = solve(f, b)

# Or explicitly
factor!(f)
x = solve(f, b)
nothing # hide
```

### Factorization types

| Type | Use case |
|------|----------|
| `SparseFactorizationQR` | Default for non-symmetric matrices |
| `SparseFactorizationCholesky` | Default for symmetric positive definite |
| `SparseFactorizationLDLT` | Symmetric indefinite (default LDLT) |
| `SparseFactorizationLDLTUnpivoted` | Symmetric indefinite, no pivoting |
| `SparseFactorizationLDLTSBK` | Symmetric indefinite, Bunch-Kaufman |
| `SparseFactorizationLDLTTPP` | Symmetric indefinite, threshold partial pivoting |
| `SparseFactorizationCholeskyAtA` | Cholesky of A'A (for least squares) |

### Solve functions

| Function | Description |
|----------|-------------|
| [`AAFactorization`](@ref AppleAccelerate.AAFactorization) | Lazy factorization wrapper |
| `solve(f, b)` | Solve `Ax = b`, returns new vector/matrix |
| `solve!(f, xb)` | Solve in-place (`xb` is overwritten with solution) |
| `f \ b` | Equivalent to `solve(f, b)` |
| `ldiv!(f, xb)` | Equivalent to `solve!(f, xb)` |
| `ldiv!(x, f, b)` | Solve `Ax = b`, store result in `x` |
| `factor!(f)` | Explicitly compute the factorization |
| `factor!(f, type)` | Compute factorization with specific type |
| `factorize(A::AASparseMatrix)` | Create an `AAFactorization` from a sparse matrix |

## Benchmarks

Performance comparison of Apple Sparse Solvers vs SuiteSparse (CHOLMOD/UMFPACK). The sparse benchmark script runs SuiteSparse first (before loading AppleAccelerate), so SuiteSparse uses OpenBLAS internally, then loads AppleAccelerate and re-runs with Apple's sparse solvers. Run with `julia --project=test/bench test/bench/run_benchmarks.jl sparse` ([source](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl)).

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
    Apple M2 Max, macOS 26, single-threaded. Julia 1.12.5, AppleAccelerate v0.6.0, SparseArrays v1.12.0 (SuiteSparse 7.8.3). Times are minimum of 5 trials. Matrices have density 0.01. Run [`bench_sparse.jl`](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/blob/master/test/bench/bench_sparse.jl) to reproduce.

!!! warning "Benchmark limitations"
    These benchmarks use random sparse matrices (`sprandn`) which lack the structure found in real-world problems (e.g., banded, block-diagonal, or mesh-derived sparsity patterns). The matrix sizes tested (N up to 5,000–50,000) are also modest by sparse solver standards. Performance on structured problems from applications like FEM, circuit simulation, or graph analysis may differ significantly.
