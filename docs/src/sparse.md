# Sparse Linear Algebra (libSparse)

AppleAccelerate wraps Apple's [Sparse Solvers](https://developer.apple.com/documentation/accelerate/sparse_solvers) library for sparse matrix operations and direct solvers.

```@setup sparse
using AppleAccelerate, SparseArrays, LinearAlgebra
```

## AASparseMatrix

A wrapper around Apple's [`SparseMatrix`](https://developer.apple.com/documentation/accelerate/sparsematrix_double) format. Construct it either from Julia's `SparseMatrixCSC` or directly from coordinate (COO) triplets.

```@example sparse
using AppleAccelerate, SparseArrays
import AppleAccelerate: AASparseMatrix, muladd!

A_jl = sprandn(100, 100, 0.05)
A = AASparseMatrix(A_jl)

x = randn(100)
y = A * x  # Sparse matrix-vector multiply
nothing # hide
```

You can also build straight from coordinate (COO) triplets with 1-based indices,
like `SparseArrays.sparse(I, J, V, m, n)` (duplicate coordinates are summed):

```@example sparse_coo
using AppleAccelerate, SparseArrays
import AppleAccelerate: AASparseMatrix

I = [1, 2, 3, 1]; J = [1, 2, 3, 3]; V = [10.0, 20.0, 30.0, 5.0]
B = AASparseMatrix(I, J, V, 3, 3)
@assert SparseMatrixCSC(B) ≈ sparse(I, J, V, 3, 3)
nothing # hide
```

The `SparseMatrixCSC` constructor automatically detects symmetric/Hermitian and
triangular structure and sets the appropriate Apple Accelerate attributes. The COO
constructor uses Accelerate's `SparseConvertFromCoordinate` and accepts
`Float32`/`Float64` and (macOS 15.5+) `ComplexF32`/`ComplexF64` values.

### Matrix operations

| Function | Description |
|----------|-------------|
| [`AASparseMatrix`](@ref AppleAccelerate.AASparseMatrix) | Construct from a `SparseMatrixCSC` or from COO triplets `(I, J, V, m, n)` |
| `A * x` | Sparse matrix-vector or matrix-matrix multiply |
| `alpha * A * x` | Scaled sparse multiply |
| [`muladd!`](@ref AppleAccelerate.muladd!) | Multiply-add: `y += A * x` or `y += alpha * A * x` |
| `transpose(A)` | Transpose (sets flag, no copy) |
| `adjoint(A)` / `A'` | Conjugate transpose (complex; equals `transpose` for real) |

### Query functions

| Function | Description |
|----------|-------------|
| `size(A)` | Matrix dimensions |
| `eltype(A)` | Element type |
| `issymmetric(A)` | Check if symmetric |
| `istriu(A)` | Check if upper triangular |
| `istril(A)` | Check if lower triangular |
| `A[i, j]` | Element access |

## Complex-valued matrices

The **entire** sparse surface accepts `ComplexF32`/`ComplexF64` in addition to
`Float32`/`Float64` (complex requires **macOS 15.5+**): the `SparseMatrixCSC` and
COO constructors, direct factorizations (Cholesky/LDLᵀ/LU/QR), `solve`/`solve!`/
`ldiv!`, `refactor!`, `SparseMultiply`/`muladd!`, transpose/adjoint, the iterative
solvers (`:cg`/`:gmres`/`:lsmr`) with preconditioners, sub-factor extraction, the
preallocated-workspace solve, partial-LU update, and introspection.

Complex matrices differ from real ones in two ways:

  - **Cholesky and CG use the Hermitian path.** For a complex matrix, `cholesky(A)`
    and `solve(A, b; method = :cg)` require `A` to be **Hermitian** positive-definite
    (`A == A'`), and dispatch to libSparse's Hermitian factorization — not the
    (complex-)symmetric one. A `SparseMatrixCSC{Complex}` that satisfies `ishermitian`
    is auto-tagged Hermitian when wrapped.
  - **`adjoint` (`A'`) and `transpose` are distinct.** `adjoint` conjugates as well as
    transposes; both are attribute-flag views that share the underlying CSC data.

```@example sparse_complex
using AppleAccelerate, SparseArrays, LinearAlgebra
import AppleAccelerate: AASparseMatrix, cholesky, solve

M = sprandn(ComplexF64, 60, 60, 0.05)
H = AASparseMatrix(SparseMatrixCSC(M * M' + 60I))   # Hermitian positive-definite
b = randn(ComplexF64, 60)
x = solve(cholesky(H), b)                            # Hermitian Cholesky
@assert x ≈ Matrix(H) \ b
nothing # hide
```

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

## Iterative solvers (CG / GMRES / LSMR)

Krylov iterative solvers are available through `solve` with a `method` keyword,
dispatching on `AASparseMatrix` (or a `SparseMatrixCSC` directly). Choose `:cg`
for symmetric positive-definite systems, `:gmres` for square non-symmetric or
indefinite systems, and `:lsmr` for rectangular or singular least-squares systems.

```@example sparse
import AppleAccelerate: AASparseMatrix, solve

M = sprandn(200, 200, 0.02)
A = AASparseMatrix(M * M' + 200I)    # SPD
b = randn(200)
x = solve(A, b; method = :cg, rtol = 1e-10)
nothing # hide
```

`atol`/`rtol` set the convergence tolerances (0 selects the library default),
`maxiter` caps iterations (0 → 100). GMRES accepts `variant`
(`:dqgmres`/`:gmres`/`:fgmres`) and `nvec`; LSMR accepts `lambda` (Tikhonov
damping) and `nvec`.

### Preconditioners

Pass `preconditioner = :diagonal` or `:diagscaling` to `solve`, or build a reusable
[`AAPreconditioner`](@ref AppleAccelerate.AAPreconditioner) handle:

```@example sparse
import AppleAccelerate: AASparseMatrix, AAPreconditioner, solve

A = AASparseMatrix(let M = sprandn(200, 200, 0.02); M * M' + 200I end)
P = AAPreconditioner(A; kind = :diagonal)
x = solve(A, randn(200); method = :cg, preconditioner = P, rtol = 1e-10)
nothing # hide
```

## Preallocated / thread-safe solve workspace

For repeated or concurrent solves that reuse a factorization, a `solve!` overload
takes a caller-owned scratch buffer sized with
[`solve_workspace_size`](@ref AppleAccelerate.solve_workspace_size), avoiding the
per-call internal allocation. Each concurrent thread must use its own workspace and
its own solution buffer.

```@example sparse
import AppleAccelerate: AAFactorization, factor!, solve!, solve_workspace_size

f = AAFactorization(sprandn(100, 100, 0.1) + 20I)
factor!(f)
ws = Vector{UInt8}(undef, solve_workspace_size(f, 1))
b = randn(100); x = similar(b)
solve!(f, b, x, ws)
nothing # hide
```

## Sub-factor extraction (Q, R, L, D, P)

Individual factors of a factorization can be extracted with
[`subfactor`](@ref AppleAccelerate.subfactor) and applied with `*` (multiply) or
`\` (solve). For example, from a QR factorization, `R` is upper-triangular and `Q`
is orthogonal:

```@example sparse
import AppleAccelerate: AAFactorization, factor!, subfactor,
                        SparseFactorizationQR, SparseSubfactorR, SparseSubfactorQ

A = sprandn(40, 15, 0.3)
f = AAFactorization(A)
factor!(f, SparseFactorizationQR)
R = subfactor(f, SparseSubfactorR)   # 15×15 upper-triangular
Q = subfactor(f, SparseSubfactorQ)   # 40×15 orthogonal
y = randn(15)
@assert R * (R \ y) ≈ y
nothing # hide
```

## Partial LU update

For a **pivotless** LU factorization, [`update_partial_lu!`](@ref AppleAccelerate.update_partial_lu!)
applies a partial refactorization when only a few entries change — recomputing only
the L/U values a from-scratch LU would alter (requires macOS 15.5+). Distinct from
`refactor!`, which recomputes the entire numeric factorization.

## Introspection

[`numeric_options`](@ref AppleAccelerate.numeric_options) and
[`symbolic_options`](@ref AppleAccelerate.symbolic_options) read back the options
libSparse recorded for a completed factorization.

### Iterative / advanced solve functions

| Function | Description |
|----------|-------------|
| `solve(A::AASparseMatrix, b; method, …)` | Iterative CG/GMRES/LSMR solve |
| [`AAPreconditioner`](@ref AppleAccelerate.AAPreconditioner) | Diagonal / diagonal-scaling preconditioner |
| [`solve_workspace_size`](@ref AppleAccelerate.solve_workspace_size) | Bytes needed for a preallocated-workspace solve |
| `solve!(f, b, x, ws)` | Solve reusing a caller-owned workspace buffer |
| [`subfactor`](@ref AppleAccelerate.subfactor) | Extract a Q/R/L/D/P sub-factor to apply with `*` / `\` |
| [`update_partial_lu!`](@ref AppleAccelerate.update_partial_lu!) | Partial LU refactorization for a low-rank change |
| [`numeric_options`](@ref AppleAccelerate.numeric_options) | Read back numeric-factor options |
| [`symbolic_options`](@ref AppleAccelerate.symbolic_options) | Read back symbolic-factor options |

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.AAFactorization
AppleAccelerate.muladd!
AppleAccelerate.factor!
AppleAccelerate.solve
AppleAccelerate.solve!
AppleAccelerate.refactor!
AppleAccelerate.AAPreconditioner
AppleAccelerate.solve_workspace_size
AppleAccelerate.subfactor
AppleAccelerate.update_partial_lu!
AppleAccelerate.numeric_options
AppleAccelerate.symbolic_options
```
