---
name: accelerate-wrapper
description: >
  General method for building out AppleAccelerate.jl's coverage of the Apple
  Accelerate framework — for capabilities that do NOT exist yet as well as fixing
  ones that do. Given any Accelerate C symbol (vDSP, vForce, libSparse, BNNS,
  Quadrature, or a subframework not yet in scope), it takes you from scanning the
  header surface → generating the FFI binding → deriving a correct idiomatic Julian
  API from the C signature → cross-validated tests → docs → benchmarks, as one
  coherent pipeline. Use for adding new wrapped functions or whole new API families,
  auditing coverage, or fixing correctness / memory-safety bugs in any wrapper.
---

# Building out Accelerate coverage — a general recipe

This is a **method**, not a catalog. It works for a function whose family already has
wrappers *and* for the first wrapper of a brand-new subframework. The existing modules
are **worked examples of a few archetypes** — identify which archetype your target
matches, then follow that pattern. Work one coherent family at a time through all six
phases before moving on; a family scanned → bound → wrapped → tested → documented →
benchmarked is one shippable unit. Don't batch-generate half-finished wrappers.

## The two layers

1. **Raw layer — `src/lib/LibAccelerate.jl`** (Clang.jl-generated, never hand-edited).
   Correct ABI types already (`vDSP_Stride=Int64`, `vDSP_Length=UInt64`,
   `Ptr{DSPSplitComplex}`, struct layouts, enums). This is your source of truth for
   every signature.
2. **Idiomatic layer — the hand-written `src/*.jl`.** Calls `LibAccelerate.foo(...)`,
   owns Julia types, allocating + mutating (`!`) variants, docstrings. New work lands
   here.

---

## Phase 1 — Scan the surface

Authoritative contract is the **on-disk SDK header** (no network; often richer than
the web page — e.g. exact buffer-size formulas). Apple's web docs
(`developer.apple.com/documentation/accelerate/<fn>`) supplement it.

```bash
SDK=$(xcrun --show-sdk-path)
ACC="$SDK/System/Library/Frameworks/Accelerate.framework"
ls "$ACC/Frameworks"                       # the subframeworks: today only vecLib + vImage
VECLIB="$ACC/Frameworks/vecLib.framework/Headers"
ls "$VECLIB" "$VECLIB/Sparse" "$VECLIB/BNNS" "$VECLIB/Quadrature"
```

In-scope headers today (`gen/generate.jl`): `vDSP.h`, `vForce.h`, `vBasicOps.h`,
`vfp.h`, `vectorOps.h`, `vBigNum.h`, `Sparse/Solve.h`, `BNNS/bnns.h`,
`Quadrature/Quadrature.h`. Accelerate has exactly **two** subframeworks on disk —
`vecLib` (all of the above) and **`vImage`** (image processing: convolution,
geometry/resize/warp, format conversion, histogram, morphology, alpha compositing,
Core Video interop), which is entirely out of scope. **To bring a new capability area
in scope** (vImage, or a future subframework), add its umbrella header to the
`headers` list in `gen/generate.jl` and regenerate (Phase 2).

**Coverage-gap scan.** Don't use a literal `comm -23` of raw-vs-referenced symbols —
it over-reports catastrophically (wrappers build names dynamically, e.g.
`Symbol(string("vDSP_", op, suff))` and from symbol lists like `(:vadd,:vsub,:vmul)`,
so most *wrapped* functions look "missing"). Reconcile at the **base-operation** level:
collect every identifier token from the idiomatic source, then for each raw function
strip its family prefix and type suffix and test the stem against that token set.

```bash
grep -oE 'function [A-Za-z_][A-Za-z0-9_]+' src/lib/LibAccelerate.jl | sed 's/^function //' | sort -u > /tmp/available.txt
cat src/*.jl > /tmp/idiom.txt          # idiomatic layer only — NOT src/lib/
```
```python
import re
avail  = [l.strip() for l in open('/tmp/available.txt') if l.strip()]
tokens = set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', open('/tmp/idiom.txt').read()))
def stems(n):
    c = {n}
    for p in ('vDSP_','vv','BNNS','bnns_','Sparse','_Sparse','sparse_','quadrature_'):
        if n.startswith(p): c.add(n[len(p):])
    for x in list(c):
        for s in ('D','f','_Float','_Double','_Complex_Float','_Complex_Double'):
            if x.endswith(s): c.add(x[:-len(s)])
    return {x for x in c if x}
for n in avail:                        # survivors = genuine candidates to investigate
    if not any(s in tokens for s in stems(n)): print(n)
```

This scored vForce at 100% and cut false positives ~10× versus `comm`. Survivors are
still only *candidates* — verify each by reading the module. Judge **opaque-handle
families at the capability level, not the symbol level**: Sparse and BNNS call *public*
mangled/umbrella symbols, so the private `_Sparse*` implementations and per-layer
`BNNSFilterCreate*` names show as "unwrapped" even when the capability is fully covered
(or, conversely, is exposed only through dead binding code). Confirm by grepping the
base op and reading how the wrapper reaches it.

**Decide whether it deserves an idiomatic wrapper.** Wrap what has a natural Julian
surface and what users reach for — not everything that exists. Register-width SIMD
(`vfp`, `vectorOps`, `vBasicOps`, `vBigNum`) is deliberately left to the raw layer
(rationale in `vmath.jl` `VMATH_COVERAGE`); document any such intentional exclusion
so the next person doesn't re-investigate it.

## Phase 2 — Generate / extend the FFI bindings

The symbol is usually already generated. Regenerate only when it's genuinely absent
or a header changed, or when you added a header in Phase 1:

```bash
julia --project=gen gen/generate.jl     # reproducible; fails loudly on drift (#186)
```

`gen/generator.toml` = module/library/enum settings + an `output_ignorelist` for
macros Clang.jl mistranslates. `gen/generate.jl` = header list + a dead-symbol strip
(dlsym against the live framework). **Always review the raw-layer diff** — a changed
struct layout silently corrupts every call that uses it.

## Phase 3 — Derive the idiomatic API from the archetype

First, **classify your target's API shape** and copy the matching module's pattern.
This is what lets you wrap something with no local precedent:

| Archetype | Shape | Template | Key concerns |
|-----------|-------|----------|--------------|
| Stateless kernel | pointers in/out + length | `array.jl`, `complexarray.jl` | length guards, operand order, strides |
| Setup + execute | create opaque plan once, run many times | FFT/DCT/DFT in `dsp.jl` | cache/finalize the setup, scaling/packing |
| Opaque-handle resource | create → use → destroy | `Biquad`/`BiquadMulti` (dsp), `AAFactorization` (sparse) | `mutable struct` + `finalizer`, keep backing arrays rooted, cleanup fn |
| Struct-descriptor + apply | build C descriptor structs, call apply | `bnns.jl` | exact struct field order/type from raw layer, `Ref`, GC-safe descriptor |
| Caller-owned workspace | pass storage/scratch buffers you size | `SparseConvertFromCoord` (sparse) | size buffers per the header formula, keep them alive, no double-free |
| Callback / trampoline | C calls back into Julia | `integrate` (quadrature) | top-level `@cfunction` (no closures on aarch64<1.12), carry errors out via a boxed ref, `GC.@preserve` the box |

Conventions across all: `for (T, suff) in ((Float32, ""), (Float64, "D"))` + `@eval`
(vDSP suffix `""`/`"D"`; vForce `vv*` suffix `"f"`/`""`); provide allocating `foo`
and mutating `foo!` with the output first; docstring every public name.

### Read the C signature for hazards (general safety checklist)

Every wrapper must clear all of these. They are *classes of hazard* — the parenthetical
functions are only examples. Nearly every historical bug is one of these (the git log
is full of "Harden FFI/GC safety", "Guard against OOB writes", "Add length guards"):

1. **Root every owning object across the ccall** with `GC.@preserve` — arrays behind a
   bare `pointer(...)`, `Ptr` arrays passed as `Ptr{Ptr{T}}`, arrays behind a
   split-complex view, caller-owned workspace, the box behind a `@cfunction` user
   pointer. Build pointer arrays *inside* the preserve block.
2. **Validate every operand length** the call reads/writes before calling — a C
   function sizes all operands from one count arg; a short operand is a silent OOB on
   a raw pointer. Throw `DimensionMismatch`.
3. **Confirm operand order and sign** from the raw signature + header, never from
   intuition or a copied comment (some vDSP ops reverse subtract/divide operands,
   some don't; comments claiming a swap are often wrong — the test is truth).
4. **Match the length/stride argument *kind*, and support strided inputs.** vDSP
   strides/lengths are Int64/UInt64 (Julia `Int` auto-converts); vForce `vv*` take
   the length as a **pointer** (`Ref{Cint}(n)`). Read the raw-layer arg types; don't
   assume by-value. Where the C call takes a stride, accept `StridedVector`/
   `StridedMatrix` (not just `Vector`) and pass `stride(X, dim)` — never a hard-coded
   `1` — so column views, reshapes, and slices work copy-free (see #191); the length
   arg stays the *logical* element count. Guard the ops that can't take an arbitrary
   stride: split-complex packing, matrix/FFT paths, and index-returning reductions
   (vDSP returns a stride-scaled *memory offset*, so require a positive stride via
   `_check_positive_stride` and recover the index with `div(offset, stride)+1`).
   Still `GC.@preserve` across the ccall when passing a view's `pointer` — it roots
   the parent buffer.
5. **Respect layout.** vDSP matrices are **row-major**, Julia is **col-major** — many
   matrix ops need operands and dims passed swapped (see `zmmul`, the confirmed-wrapped
   template). An op that combines a multiply with a subtraction on the multiplicand
   (`D=(A−B)*C`) doesn't compose with the transpose swap without materializing an
   intermediate — do the extra copy or leave it unwrapped rather than ship it wrong.
6. **Verify struct layout field-by-field** against the raw layer for any struct you
   construct (descriptors, options, factor structs) — order/type/alignment must match
   exactly, or it silently corrupts.
7. **Manage resource lifetime.** Opaque handles need a `mutable struct` + `finalizer`
   calling the C destroy/cleanup; caller-owned storage must outlive every use; never
   free memory Accelerate didn't allocate. Cached read-only setups (FFT) may be shared;
   others may not — check the header's threading note.
8. **Know the numeric contract.** Scaling and packing conventions differ per family
   (e.g. vDSP forward real FFT is 2× and packs DC/Nyquist together); encode a
   `is_supported_*` predicate where the C call silently fails on unsupported inputs.

## Phase 4 — Test: cross-validate against an independent package

Assert against a reference, for **both** Float32 and Float64 (Float32 with
`rtol=sqrt(eps(T))`):

| Area | Reference |
|------|-----------|
| array / complexarray | `Base`, `Statistics`, `LinearAlgebra` |
| dsp / FFT | `FFTW` |
| sparse | `SparseArrays.sparse`, dense `\` |
| quadrature | analytic integrals |
| BNNS / NN | plain-Julia reference, `NNlib` |
| a new area | whatever established Julia package computes the same thing |

Also test guards (`@test_throws DimensionMismatch`/`ArgumentError`) and round-trip
identities. Choose **discriminating** inputs — non-square, non-symmetric, unsorted,
partial-length — so a layout/order/lifetime bug fails loudly (a sorted or square-only
test hides real bugs; use FFTW carefully — `MEASURE` overwrites its input during
planning, so never compare against a copy you handed to the plan).

## Phase 5 — Document

Add the function to the page's cross-ref table (`[`foo`](@ref AppleAccelerate.foo)`),
the `@docs` block, and ideally a runnable `@example`. New capability area → new
`docs/src/<area>.md` page (register it in `docs/make.jl`). Gotcha: same-label
`@example` blocks share one module across a page — don't bind a name later blocks
rely on (`I = [...]` shadows `LinearAlgebra.I`). Keep README headline claims in sync
with `benchmarks.md`.

**ALWAYS build the docs locally before pushing — `julia --project=docs docs/make.jl` — and
confirm it exits 0.** The Documentation build is a **separate CI gate from `Pkg.test()`**:
a green test suite says nothing about the docs, and a docs failure blocks the PR. The most
common failure (it has bitten multiple PRs) is Documenter's `cross_references` check:
**every name you `@ref` must also appear in an `@docs` block**, or the build aborts with
`Cannot resolve @ref for ...`. That includes refs in a page's cross-ref *table*, refs inside
a **docstring's own body** (e.g. a "See also [`foo!`](@ref)" line — `foo!` must be `@docs`'d
too), and functions that share a docstring via `@doc (@doc x) y` (both `x` and `y` need
their own `@docs` entry to be linkable). `@example`/`jldoctest` blocks are also *executed*
during the build, so a runtime error or a wrong expected output fails it. When you add a new
page, wire it into `docs/make.jl` AND make sure every `@ref` it introduces resolves — a
subagent that "didn't build docs" is the usual cause of a red Documentation check on an
otherwise-green PR.

## Phase 6 — Benchmark

Add cases to the matching `test/bench/bench_*.jl` (or a new one, registered in
`run_benchmarks.jl`); they print Markdown tables. Each suite runs in its **own fresh
process** so the reference library (OpenBLAS/SuiteSparse/FFTW) is measured *before*
AppleAccelerate loads. Single-threaded, min-of-N. Transcribe into
`docs/src/benchmarks.md` preserving each table's rows/rounding. Re-measure any
surprising ratio in isolation and confirm numerical correctness before publishing it.

---

## Verify — exact commands (macOS)

`julia` is at `~/.juliaup/bin` (put it on PATH). `timeout` does **not** exist on
macOS — use `run_in_background`. Benchmarks are timing-sensitive; don't run other
heavy Julia processes alongside them.

```bash
export PATH="$HOME/.juliaup/bin:$PATH"
julia --project=. -e 'using Pkg; Pkg.test()'                       # full suite (authoritative)
julia --project=test/bench test/bench/run_benchmarks.jl [fft|sparse|dense|array]
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'  # once
julia --project=docs docs/make.jl                                  # runs @example + doctests
julia --project=gen  gen/generate.jl                               # only when the raw layer must change
```

## Reference

Durable, empirically-verified ABI facts (packing layouts, buffer-size formulas,
operand-order findings, the wrapped-function inventory) live in the project
auto-memory `MEMORY.md`. Consult it before re-deriving ABI details; add to it
whenever you verify something non-obvious so the next capability build reuses it.
**But verify its inventory claims against the source before trusting them** — the
"what's wrapped" list drifts (a note claiming `zmma`/`zmms` were added was found to be
false; grep of the module settled it). ABI/packing facts are durable; coverage claims
are not. A grep of the actual module beats a stale memory line every time.
