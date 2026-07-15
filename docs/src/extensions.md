# [Architecture](@id architecture)

This page explains how AppleAccelerate is structured internally, and why it
integrates with the rest of the ecosystem through its own namespace rather than
by hooking other packages' interfaces.

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

## No package extensions

AppleAccelerate ships **no package extensions**, and defines no methods on other
packages' functions. Everything it provides lives under the `AppleAccelerate.`
namespace, and the package exports nothing.

This is a deliberate choice. Earlier versions did hook
[AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl) and
[NNlib.jl](https://github.com/FluxML/NNlib.jl) through extensions, which meant
that merely loading AppleAccelerate — even for an unrelated feature like BLAS
forwarding — silently changed what `fft`, `plan_fft`, or `batched_mul!` did
elsewhere in the session. Because Accelerate's kernels only cover a subset of the
input space (vDSP's FFT is power-of-2 and 1-D/2-D only), those methods won
dispatch and then failed on inputs the general backend would have handled fine.
See [issue #139](https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl/issues/139).

The consequence for you is simple and predictable:

```julia
using AppleAccelerate, FFTW

fft(x)                  # always FFTW — loading AppleAccelerate changes nothing
AppleAccelerate.fft(x)  # explicitly vDSP, for power-of-2 1-D/2-D input
```

You choose Accelerate per call site, by name. Nothing is intercepted, and
benchmarking one against the other is a matter of changing the prefix.

## Why `LinearAlgebra` and `SparseArrays` are hard dependencies

[`LinearAlgebra`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) and
[`SparseArrays`](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) are
regular `[deps]`, not weak dependencies:

1. They are standard libraries that are effectively always present in the Julia
   sysimage, so making them weak would save nothing.
2. BLAS/LAPACK forwarding needs `LinearAlgebra` at **load time** (`__init__`
   installs Accelerate into libblastrampoline), so it cannot be deferred.
3. Keeping `SparseArrays` a hard dependency lets `AAFactorization` subtype
   `LinearAlgebra.Factorization` and lets `AASparseMatrix` interoperate with
   `SparseMatrixCSC` directly.

BLAS/LAPACK forwarding is the one place AppleAccelerate does change global
behavior — that is the documented purpose of loading it, and it goes through
libblastrampoline, a mechanism designed for exactly this.
