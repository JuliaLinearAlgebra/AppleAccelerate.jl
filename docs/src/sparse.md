# Sparse Linear Algebra

AppleAccelerate wraps Apple's [Sparse Solvers](https://developer.apple.com/documentation/accelerate/sparse_solvers) library for sparse matrix operations and direct solvers.

```@setup sparse
using AppleAccelerate, SparseArrays, LinearAlgebra
```

## AASparseMatrix

A wrapper around Apple's [`SparseMatrix`](https://developer.apple.com/documentation/accelerate/sparsematrix_double) format, constructed from Julia's `SparseMatrixCSC`.

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

Wraps Apple's [`SparseOpaqueFactorization`](https://developer.apple.com/documentation/accelerate/sparseopaquefactorization_double).
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

## Factorization conveniences

For matrices wrapped as an `AASparseMatrix`, the usual `LinearAlgebra`
factorization spellings build an `AAFactorization` of the requested kind:

```@example sparse
import AppleAccelerate: AASparseMatrix

A = AASparseMatrix(sprandn(100, 100, 0.1) + 20I)

F = lu(A)        # AAFactorization (LU; requires macOS 15.5+)
G = qr(A)        # QR (also works for rectangular / least-squares)
nothing # hide
```

`cholesky(A)` and `ldlt(A)` are available for symmetric/Hermitian matrices.

!!! note "No type piracy"
    These methods dispatch on `AppleAccelerate`'s own `AASparseMatrix`, **not** on
    `SparseMatrixCSC`, so they do not override Julia's built-in
    `lu(::SparseMatrixCSC)` / `cholesky` / `qr` / `ldlt` (UMFPACK/CHOLMOD/SPQR). You
    opt in to the Accelerate solvers by wrapping your matrix in `AASparseMatrix`
    first, e.g. `lu(AASparseMatrix(A))`.

## Reusing a factorization: `refactor!`

When a matrix changes its **values but not its sparsity pattern** — the common case
in Newton iterations, implicit time stepping, and parameter sweeps —
[`refactor!`](@ref AppleAccelerate.refactor!) recomputes the numeric factorization
in place while **reusing the existing symbolic factorization** (the fill-reducing
ordering and sparsity analysis). This is substantially cheaper than building a fresh
`AAFactorization`.

```@example sparse
import AppleAccelerate: AAFactorization, refactor!, solve

P = sprandn(100, 100, 0.05) + 20I
F = AAFactorization(P)
solve(F, randn(100))             # forces the first (full) factorization

# Same sparsity pattern, new values:
P2 = copy(P); P2.nzval .*= 1.5
refactor!(F, P2)                 # reuses the symbolic factorization
x = solve(F, randn(100))
nothing # hide
```

`F` must already hold a completed factorization (call `factor!`/`solve` once first),
and `A` must have the same number of stored nonzeros as the original matrix.

For LU/Cholesky/LDLᵀ factorizations you can equivalently use the `LinearAlgebra`-style
spellings that mirror how `SparseArrays` exposes symbolic reuse — `lu!(F, A)`,
`cholesky!(F, A)`, `ldlt!(F, A)`. Each requires `F` to already hold a factorization
of the matching kind and delegates to `refactor!`. QR reuse has no stdlib spelling,
so use `refactor!` directly for it.

## Round-tripping to `SparseMatrixCSC`

An `AASparseMatrix` can be materialized back to a standard Julia `SparseMatrixCSC`.
The transpose/adjoint/symmetric/Hermitian/triangular attribute bits are honored, so
the result equals the logical matrix the wrapper represents:

```@example sparse
A = AASparseMatrix(sprandn(50, 50, 0.1))
B = SparseMatrixCSC(A)
@assert B ≈ SparseMatrixCSC(A)   # round-trips the logical matrix
nothing # hide
```

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.AAFactorization
AppleAccelerate.muladd!
AppleAccelerate.factor!
AppleAccelerate.solve
AppleAccelerate.solve!
AppleAccelerate.refactor!
```
