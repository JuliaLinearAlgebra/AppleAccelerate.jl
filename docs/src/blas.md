# Dense Linear Algebra (BLAS / LAPACK)

AppleAccelerate forwards [BLAS](https://developer.apple.com/documentation/accelerate/blas) and [LAPACK](https://developer.apple.com/documentation/accelerate/solving-systems-of-linear-equations-with-lapack) calls to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate) via Julia's [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) (LBT) mechanism. This happens automatically when the package is loaded.

## How it works

On `__init__`, AppleAccelerate loads both LP64 and ILP64 BLAS/LAPACK interfaces from Accelerate. OpenBLAS remains as a fallback for operations not provided by Accelerate (e.g., `gemmt`).

Since Accelerate provides a full [BLAS](https://developer.apple.com/documentation/accelerate/blas) and [LAPACK](https://developer.apple.com/documentation/accelerate/solving-systems-of-linear-equations-with-lapack) implementation, all standard Julia [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) operations are accelerated transparently. This includes:

**Factorizations:** `lu`, `qr`, `cholesky`, `svd`, `eigen`, `schur`, `ldlt`, `hessenberg`, and their in-place `!` variants

**Solvers:** `\`, `ldiv!`, `rdiv!`

**Matrix operations:** `mul!`, `*`, `det`, `tr`, `inv`, `pinv`, `rank`, `norm`, `cond`, `opnorm`

**BLAS routines:** All Level 1 (vector), Level 2 (matrix-vector), and Level 3 (matrix-matrix) operations via `LinearAlgebra.BLAS`

For the complete list of available operations, see the [Julia LinearAlgebra documentation](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/).

## Loading

```@example
using AppleAccelerate
```

| Function | Description |
|----------|-------------|
| [`load_accelerate`](@ref AppleAccelerate.load_accelerate) | Load Accelerate BLAS/LAPACK via LBT |

## Threading

On macOS 26+, you can control BLAS threading:

| Function | Description |
|----------|-------------|
| [`set_num_threads`](@ref AppleAccelerate.set_num_threads) | Set the number of Accelerate BLAS threads |
| [`get_num_threads`](@ref AppleAccelerate.get_num_threads) | Get the number of Accelerate BLAS threads |

## Utilities

| Function | Description |
|----------|-------------|
| [`get_macos_version`](@ref AppleAccelerate.get_macos_version) | Return the current macOS version |

```@docs
AppleAccelerate.load_accelerate
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
AppleAccelerate.get_macos_version
AppleAccelerate._read_macos_version
```
