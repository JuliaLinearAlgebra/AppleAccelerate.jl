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

Generation is deterministic: re-running with the same SDK produces byte-identical
output, so CI can regenerate and diff to catch drift when a new SDK adds symbols.

## Scope

Currently generates only the **Quadrature** subframework, as a proof of concept.
To extend coverage, append headers to the `headers` list in
[`generate.jl`](./generate.jl) (e.g. `BNNS/bnns.h`, `Sparse/Solve.h`,
`vForce.h`). The BLAS/LAPACK path is intentionally excluded — it is forwarded
through libblastrampoline, not ccall.

## Files

| File | Purpose |
|------|---------|
| `generate.jl` | Entry point: resolves SDK paths, runs Clang.jl |
| `generator.toml` | Clang.jl options (module name, library, enum style, …) |
| `prologue.jl` | Spliced into the generated module — defines `libacc` |
| `Project.toml` | Pins Clang.jl for the generator environment |
