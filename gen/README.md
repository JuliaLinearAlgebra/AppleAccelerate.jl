# `gen/` — Clang.jl bindings generator

This directory regenerates the raw ABI layer in
[`../src/lib/LibAccelerate.jl`](../src/lib/LibAccelerate.jl) directly from Apple's
Accelerate C headers, using [Clang.jl](https://github.com/JuliaInterop/Clang.jl).

## Why

Accelerate's struct- and enum-heavy subframeworks (Quadrature, Sparse, BNNS) are
painful and error-prone to bind by hand — every C struct field offset and enum
value has to be transcribed and kept in sync with the SDK. Clang.jl reads the
headers and emits those definitions automatically, so the hand-written code in
`src/` only has to provide ergonomics, never ABI details.

The result is a conventional two-layer wrapper:

- **`src/lib/LibAccelerate.jl`** — generated, committed, *not* hand-edited. A 1:1
  mirror of the C API (structs, enums, `@ccall`s).
- **`src/quadrature.jl`, `src/array.jl`, …** — hand-written idiomatic Julia that
  calls into `LibAccelerate`.

## Usage

```sh
julia --project=gen gen/generate.jl
```

This reads headers from the active macOS SDK (resolved via `xcrun --show-sdk-path`,
so it needs Xcode or the Command Line Tools) and overwrites `src/lib/LibAccelerate.jl`.
The committed output is self-contained — **end users need only the runtime
framework that ships with every Mac**, not the headers or Clang.jl.

## Reproducibility

Generation is deterministic **per machine**: re-running with the same SDK *and* the
same macOS runtime produces byte-identical output. It is not deterministic across
machines, because:

- The dead-symbol strip pass `dlopen`s the **live** framework at
  `/System/Library/Frameworks/Accelerate.framework/Accelerate` (not the SDK), so the
  set of stripped wrappers depends on which symbols the running macOS version exports.
- Clang.jl names anonymous structs with a global counter (`var"##Ctag#NNN"`). These
  names renumber **wholesale** when SDK headers add or remove any anonymous type, so
  even a small SDK update can produce a large, mechanical diff in the committed output.

There is currently no CI job that regenerates and diffs the output; drift against a
new SDK is caught by re-running the generator manually. To make silent drift harder,
`generate.jl` pins Clang.jl (via `Project.toml` `[compat]`) and every
post-processing pass **errors** if an expected
transformation finds zero matches instead of silently no-opping.

## Scope

Generates the in-scope headers listed in [`generate.jl`](./generate.jl): vDSP,
vForce, vBasicOps, vfp, vectorOps, vBigNum, `Sparse/Solve.h` (the C solver API),
BNNS (+ graph), Quadrature, and the array-based vImage headers (`vImage_Types.h`,
`Alpha.h`, `BasicImageTypes.h`, `Conversion.h`, `Convolution.h`, `Geometry.h`,
`Histogram.h`, `Morphology.h`, `Transform.h`) — ~1400 functions (~540 of them vImage),
140+ structs, 70+ enums. vImage is a separate sub-framework binary, but the Accelerate
umbrella re-exports its symbols, so its wrappers share `libacc` with everything else. Extend
coverage by appending headers to that list.

**Intentionally excluded:**
- BLAS/LAPACK (`cblas*.h`, `lapack*.h`, …) — forwarded via libblastrampoline, not
  ccall. They get transitively pulled in via `Sparse/Types.h`, so `generate.jl`
  post-processes the output to strip the `cblas_*`/`catlas_*`/`clapack_*` wrappers.
- `LinearAlgebra/` — C++ generics, not C-mappable.
- `Sparse/BLAS.h` dense×sparse multiply — C++ name-mangled (hand-wrapped elsewhere).
- `vImage/vImage_Utilities.h`, `vImage/vImage_CVUtilities.h` — CoreGraphics / CoreVideo
  interop (`CGImage`, `CVPixelBuffer`); no array-based surface, and they pull in the whole
  CG/CV header graph. The vImage operation headers are listed individually (not via the
  `vImage.h` umbrella) to keep these out.

## Auditing idiomatic coverage

```sh
julia --project=. gen/coverage_audit.jl          # vDSP_* (default)
julia --project=. gen/coverage_audit.jl BNNS     # any other name prefix
```

Lists the generated functions that nothing in the idiomatic layer calls. It inspects
the **lowered IR** of every method in the package, so it is exact where a text search
is not: most wrappers assemble the C name at macro-expansion time
(`Symbol(string("vDSP_vfix", intname, suff))`), which makes fully wrapped families look
"missing" to `grep`. The surviving vDSP names and the reason each is left to the raw
layer are recorded in `AppleAccelerate.VDSP_COVERAGE`. For opaque-handle families
(Sparse, BNNS) judge the result at the capability level — those wrappers reach public
umbrella symbols, so private per-type implementations legitimately show as unreferenced.

## Known limitations

- **Bitfield structs are not safe to pass by value.** Clang.jl emits any struct that
  contains a C bitfield (e.g. everything embedding `SparseAttributes_t`) as an
  `NTuple{N,UInt8}` blob with alignment 1, while the real struct is 8-aligned. The *size*
  is right, so passing such a struct by pointer is fine, but a by-value `@ccall` argument
  can be laid out wrongly. `src/sparse.jl` therefore keeps field-typed mirrors for those
  structs and asserts their size/offsets against the generated ones in its "libSparse ABI
  parity" tests. The same goes for over-aligned opaque blobs (vImage's `aligned(16)`
  Y'CbCr conversion infos).
- **`__asm__` labels are applied by a post-pass, not by Clang.jl.** bnns_graph.h links its
  entry points to `_v2` symbols while the un-suffixed symbols stay exported with the old
  argument lists; `apply_asm_labels!` in `generate.jl` retargets the `@ccall`s (and runs
  *before* the dead-symbol strip, so functions that only exist as `_v2` are kept).

- A few vDSP functions whose signatures use `arm_neon`/`simd` vector types are dropped by
  libclang under the default GCC artifact include path (we capture ~93% of vDSP, all the
  pointer/length array ops). These SIMD-typed overloads are not part of the idiomatic
  surface anyway.
- Clang.jl omits the SIMD typedef `vUInt32` but leaves references to it in the vBigNum
  union `.v` accessors. `prologue.jl` supplies `NTuple{4,VecElement{UInt32}}`, preserving
  its 16-byte size and alignment; `test/lib_tests.jl` checks the layout and accessors.
- `Sparse/Solve.h` includes the whole `<Accelerate/Accelerate.h>` umbrella when it
  resolves, which reaches CoreGraphics → CoreFoundation → libdispatch and aborts Clang.jl
  (`no definition for dispatch_queue_t's underlying type`). `generate.jl` pre-defines the
  umbrella's include guard and force-includes `cblas.h` instead, so scope is decided only
  by the `headers` list.
- Struct alignment attributes are not carried over. The two opaque vImage Y'CbCr info
  blobs are `aligned(16)` in C but come out as `NTuple{128,UInt8}`; `src/vimage.jl` keeps
  its own 16-byte-aligned definition for those and passes them as untyped pointers.
- Exported *data* symbols (e.g. `kvImage_YpCbCrToARGBMatrix_ITU_R_601_4`) are not emitted;
  `src/vimage.jl` reads them with `dlsym`.

## Files

| File | Purpose |
|------|---------|
| `generate.jl` | Entry point: resolves SDK paths, runs Clang.jl, strips out-of-scope BLAS |
| `coverage_audit.jl` | Exact list of raw functions the idiomatic layer never calls (IR walk, not grep) |
| `generator.toml` | Clang.jl options (module name, library, enum style, …) |
| `prologue.jl` | Spliced into the generated module — `libacc`, `vUInt32`, and BNNSGraph opaque handles |
| `shims/bnns_graph_shim.h` | Neutralizes availability attributes that break Clang.jl |
| `Project.toml` | Generator environment; `[compat]` pins Clang.jl |
