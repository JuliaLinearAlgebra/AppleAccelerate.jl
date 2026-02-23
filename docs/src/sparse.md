# Sparse Linear Algebra

AppleAccelerate wraps Apple's `libSparse` library for sparse matrix operations and direct solvers.

## AASparseMatrix

A wrapper around Apple's sparse matrix format, constructed from Julia's `SparseMatrixCSC`.

```julia
using AppleAccelerate, SparseArrays
import AppleAccelerate: AASparseMatrix, muladd!

A_jl = sprandn(100, 100, 0.05)
A = AASparseMatrix(A_jl)

x = randn(100)
y = A * x  # Sparse matrix-vector multiply
```

The constructor automatically detects symmetric and triangular structure and sets
the appropriate Apple Accelerate attributes.

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.muladd!
```

## AAFactorization

Lazy factorization wrapper: the factorization is computed on the first call to `solve`
or by explicitly calling `factor!`.

```julia
import AppleAccelerate: AAFactorization, solve, solve!, factor!

A = sprandn(100, 100, 0.1) + 20I
f = AAFactorization(A)

# Factorization computed lazily on first solve
b = randn(100)
x = solve(f, b)

# Or explicitly
factor!(f)
x = solve(f, b)
```

Supported factorization types:
- `SparseFactorizationQR` (default for non-symmetric)
- `SparseFactorizationCholesky` (default for symmetric positive definite)
- `SparseFactorizationLDLT`, `SparseFactorizationLDLTUnpivoted`, `SparseFactorizationLDLTSBK`, `SparseFactorizationLDLTTPP`
- `SparseFactorizationCholeskyAtA`

```@docs
AppleAccelerate.AAFactorization
```
