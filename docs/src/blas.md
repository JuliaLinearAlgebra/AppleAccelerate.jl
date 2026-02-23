# BLAS & LAPACK

AppleAccelerate forwards BLAS and LAPACK calls to Apple's Accelerate framework via Julia's [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) (LBT) mechanism. This happens automatically when the package is loaded.

## Loading

```julia
using AppleAccelerate
```

On `__init__`, AppleAccelerate loads both LP64 and ILP64 BLAS/LAPACK interfaces from Accelerate. OpenBLAS remains as a fallback for operations not provided by Accelerate (e.g., `gemmt`).

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
