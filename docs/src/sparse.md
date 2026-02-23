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
| `AASparseMatrix(M::SparseMatrixCSC)` | Construct from Julia sparse matrix |
| `A * x` | Sparse matrix-vector or matrix-matrix multiply |
| `alpha * A * x` | Scaled sparse multiply |
| `muladd!(A, x, y)` | Multiply-add: `y += A * x` |
| `muladd!(alpha, A, x, y)` | Scaled multiply-add: `y += alpha * A * x` |
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

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.muladd!
```

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
| `solve(f, b)` | Solve `Ax = b`, returns new vector/matrix |
| `solve!(f, xb)` | Solve in-place (`xb` is overwritten with solution) |
| `f \\ b` | Equivalent to `solve(f, b)` |
| `ldiv!(f, xb)` | Equivalent to `solve!(f, xb)` |
| `ldiv!(x, f, b)` | Solve `Ax = b`, store result in `x` |
| `factor!(f)` | Explicitly compute the factorization |
| `factor!(f, type)` | Compute factorization with specific type |
| `factorize(A::AASparseMatrix)` | Create an `AAFactorization` from a sparse matrix |

```@docs
AppleAccelerate.AAFactorization
```
