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

On macOS 15+, you can control BLAS threading:

```@docs
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
```

## Utilities

```@docs
AppleAccelerate.get_macos_version
```
