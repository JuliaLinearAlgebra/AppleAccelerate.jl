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

Because a byte-for-byte diff is meaningless across machines, drift is checked on the
*stable surface* instead. [`check_drift.jl`](./check_drift.jl) regenerates into a throwaway
copy of the tree (the committed file is never touched) and compares function `@ccall`
signatures, named struct fields, enum members/values and constants by name, and the
anonymous `##Ctag#` types as a multiset with the counter erased:

```sh
julia gen/check_drift.jl                      # regenerate + compare; exit 1 on drift
julia gen/check_drift.jl --generated <file>   # compare against an existing output
```

A name present on only one side is reported as informational (a newer SDK adds API; a
different macOS strips a different set of unexported symbols). A definition present on
both sides that differs — or a committed anonymous type with no match — is drift and
fails. The [`GenDrift`](../.github/workflows/GenDrift.yml) workflow runs this weekly on
`macOS-latest` (and on PRs touching `gen/` or the generated file); a failing scheduled
run opens, or comments on, a "LibAccelerate generator drift detected" issue. The check
also fails (exit 2) when `generate.jl` itself no longer completes against the SDK.

As a second line of defence, `generate.jl` pins Clang.jl (via `Project.toml` `[compat]`)
and every post-processing pass **errors** if an expected
transformation finds zero matches instead of silently no-opping.

## Scope

Generates the in-scope headers listed in [`generate.jl`](./generate.jl): vDSP,
vForce, vBasicOps, vfp, vectorOps, vBigNum, `Sparse/Solve.h` (the C solver API),
BNNS (+ graph), and Quadrature — ~900 functions, 110+ structs, 60 enums. Extend
coverage by appending headers to that list.

**Intentionally excluded:**
- BLAS/LAPACK (`cblas*.h`, `lapack*.h`, …) — forwarded via libblastrampoline, not
  ccall. They get transitively pulled in via `Sparse/Types.h`, so `generate.jl`
  post-processes the output to strip the `cblas_*`/`catlas_*`/`clapack_*` wrappers.
- `LinearAlgebra/` — C++ generics, not C-mappable.
- `Sparse/BLAS.h` dense×sparse multiply — C++ name-mangled (hand-wrapped elsewhere).
- vImage — a separate Accelerate sub-framework (not part of vecLib), not currently
  in scope.

## Known limitations

- A few vDSP functions whose signatures use `arm_neon`/`simd` vector types or
  CoreFoundation-gated availability are dropped by libclang under the default GCC
  artifact include path (we capture ~93% of vDSP, all the pointer/length array
  ops). These SIMD-typed overloads are not part of the idiomatic surface anyway.

## Files

| File | Purpose |
|------|---------|
| `check_drift.jl` | Regenerates into a temp tree and compares the stable surface with the committed file |
| `generate.jl` | Entry point: resolves SDK paths, runs Clang.jl, strips out-of-scope BLAS |
| `generator.toml` | Clang.jl options (module name, library, enum style, …) |
| `prologue.jl` | Spliced into the generated module — `libacc` + BNNSGraph opaque handles |
| `shims/bnns_graph_shim.h` | Neutralizes availability attributes that break Clang.jl |
| `Project.toml` | Generator environment; `[compat]` pins Clang.jl |
