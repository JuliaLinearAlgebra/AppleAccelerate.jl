# [Architecture & Package Extensions](@id extensions)

This page explains how AppleAccelerate is structured internally and how its
optional integrations with other packages are wired together.

## Two-layer design

AppleAccelerate is built in two layers:

- **Raw ABI layer — `AppleAccelerate.LibAccelerate`.** This submodule
  (`src/lib/LibAccelerate.jl`) is **auto-generated** with
  [Clang.jl](https://github.com/JuliaInterop/Clang.jl) directly from Apple's
  Accelerate C headers. It is a near 1:1 mirror of the C API — roughly a thousand
  `@ccall` wrappers plus the matching structs and enums — and is committed but never
  hand-edited (regenerate it with `julia --project=gen gen/generate.jl`). Because
  Accelerate's struct- and enum-heavy subframeworks (Quadrature, Sparse, BNNS) are
  painful and error-prone to bind by hand, generating this layer keeps every field
  offset and enum value in sync with the SDK automatically.

- **Idiomatic layer — the hand-written `src/*.jl` files.** These provide the
  ergonomic Julia API (`AppleAccelerate.exp`, `integrate`, `bnns_matmul`,
  `AASparseMatrix`, …). They own the niceties — Julia functions, keyword options,
  tidy return values, broadcasting — and call into `LibAccelerate` for the actual
  ABI work, never naming a C field offset or enum integer themselves.

Most users only ever touch the idiomatic layer. When you need a symbol that has no
idiomatic wrapper yet, the raw layer is available directly:

```julia
using AppleAccelerate
AppleAccelerate.LibAccelerate.some_unwrapped_symbol(args...)
```

## Package extensions

Two optional integrations are implemented as Julia
[package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)).
An extension is code that ships *inside* AppleAccelerate but only loads when a
particular **trigger package** is also loaded in the same session. This lets
AppleAccelerate integrate with those packages without forcing every user to install
them.

The mechanism, concretely:

- The extension code lives in `ext/` (`AppleAccelerateAbstractFFTsExt.jl`,
  `AppleAccelerateNNlibExt.jl`).
- The trigger packages are declared under `[weakdeps]` in `Project.toml`, and the
  `[extensions]` table maps each extension module to its trigger.
- The extension activates **automatically** the moment both packages are loaded —
  you do not call anything special:

```julia
using AppleAccelerate, NNlib   # NNlib extension activates automatically here
```

You never need `Base.get_extension(...)`; just `using` the trigger package. If the
trigger package is not installed, AppleAccelerate still loads fine — the extension
simply stays dormant, and it adds no dependency for users who do not need it.

### Which `using` activates what

| Integration | Load | Enables |
|-------------|------|---------|
| FFT via the standard Julia interface | `using AppleAccelerate, AbstractFFTs` | `plan_fft`/`plan_ifft`/`plan_bfft`, in-place plans, `mul!`, `inv`, `\` backed by vDSP — see [FFT & Transforms](@ref) |
| Neural-network ops | `using AppleAccelerate, NNlib` | `NNlib.batched_mul!` (`Float32` 3-D) and BNNS-backed activations — see [Neural Network Primitives (BNNS)](@ref) |

## Why `LinearAlgebra` and `SparseArrays` are *not* extensions

The extension mechanism applies **only** to AbstractFFTs and NNlib.
[`LinearAlgebra`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) and
[`SparseArrays`](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) are
deliberately kept as regular (hard) `[deps]`, not weak dependencies, for three
reasons:

1. They are standard libraries that are effectively always present in the Julia
   sysimage, so making them weak would save nothing.
2. BLAS/LAPACK forwarding needs `LinearAlgebra` at **load time** (`__init__`
   installs Accelerate into libblastrampoline), so it cannot be deferred behind an
   extension trigger.
3. Keeping `SparseArrays` a hard dependency lets `AAFactorization` subtype
   `LinearAlgebra.Factorization` and lets `AASparseMatrix` interoperate with
   `SparseMatrixCSC` directly — behavior that would be awkward to gate behind a
   weak dependency.
