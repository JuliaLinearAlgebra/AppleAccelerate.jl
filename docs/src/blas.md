# Dense Linear Algebra (BLAS / LAPACK)

AppleAccelerate forwards [BLAS](https://developer.apple.com/documentation/accelerate/blas) and [LAPACK](https://developer.apple.com/documentation/accelerate/solving-systems-of-linear-equations-with-lapack) calls to Apple's [Accelerate framework](https://developer.apple.com/documentation/accelerate) via Julia's [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) (LBT) mechanism. This happens automatically when the package is loaded, unless you [turn automatic forwarding off](@ref blas-opt-out).

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
| [`auto_forward_blas`](@ref AppleAccelerate.auto_forward_blas) | Whether forwarding happens automatically on load |
| [`set_auto_forward!`](@ref AppleAccelerate.set_auto_forward!) | Persistently enable/disable automatic forwarding |

## [Using the package without changing BLAS](@id blas-opt-out)

BLAS/LAPACK forwarding is global: it changes what every `LinearAlgebra` call in the session
runs on. If you only want the other subsystems (vDSP, vImage, libSparse, BNNS, …), turn
automatic forwarding off and the session's BLAS is left untouched:

```julia
using AppleAccelerate
AppleAccelerate.set_auto_forward!(false)   # writes LocalPreferences.toml; restart Julia
```

or, for a single process, set the environment variable before the package is loaded:

```sh
APPLEACCELERATE_AUTO_FORWARD=0 julia
```

The environment variable takes precedence over the preference. With automatic forwarding
off, forwarding becomes an explicit opt-in:

```julia
using AppleAccelerate, LinearAlgebra
X = AppleAccelerate.fft(randn(ComplexF64, 1024))   # vDSP; BLAS is still OpenBLAS
AppleAccelerate.load_accelerate()                  # now forward BLAS/LAPACK to Accelerate
```

`AppleAccelerate.set_auto_forward!(true)` restores the default.

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
AppleAccelerate.auto_forward_blas
AppleAccelerate.set_auto_forward!
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
AppleAccelerate.get_macos_version
AppleAccelerate._read_macos_version
```
