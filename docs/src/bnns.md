# Neural Network Primitives (BNNS)

AppleAccelerate wraps the **current, non-deprecated** slice of Apple's
[BNNS (Basic Neural Network Subroutines)](https://developer.apple.com/documentation/accelerate/bnns)
library — 61 of the ~136 `BNNS*` C entry points. The bulk of the remainder are
APIs Apple **deprecated in macOS 15 / iOS 18**: the classic filter/layer
construction API and the deprecated classic + DirectApply tensor kernels
(`BNNSMatMul`, `BNNSTile`, `BNNSGather`/`BNNSScatter`, the clip / norm family,
`BNNSOptimizerStep`, …). Those are intentionally **not wrapped** — target the
**BNNS Graph API** instead, which this page covers end to end: compile a Core ML
model, inspect its arguments, and [run inference on Julia arrays](@ref bnns-graph).
Numerically verified helpers (transpose, copy, reductions, top-k, random
generation, nearest neighbors, graph execution) are cross-checked against
plain-Julia references; the remaining thin wrappers expose the rest of the
current surface with exact FFI signatures for callers who need them.

## Element types

`Float32` works everywhere. `Float16` and the integer / `Bool` types are accepted
wherever the underlying kernel was **verified at runtime** to compute correctly —
BNNS returns a success status with wrong values for some unsupported types, so
each wrapper restricts its element types by dispatch:

| Function | Element types |
|----------|---------------|
| [`bnns_transpose`](@ref AppleAccelerate.bnns_transpose), same-type [`bnns_copy!`](@ref AppleAccelerate.bnns_copy!) | `Float16`, `Float32`, `Int8`–`Int64`, `UInt8`–`UInt64`, `Bool` |
| converting [`bnns_copy!`](@ref AppleAccelerate.bnns_copy!) | `Float16 ↔ Float32`, `Int8`/`Int16`/`Int32`/`UInt8` `→ Float32`, `Float32 → Int32` (anything else throws) |
| [`bnns_reduce`](@ref AppleAccelerate.bnns_reduce) | `Float16`, `Float32`; `Int32` for the integer-exact reductions |
| [`bnns_topk`](@ref AppleAccelerate.bnns_topk) | `Float16`, `Float32`, `Int8`, `Int16`, `Int32`, `UInt8`, `UInt16` |
| [`bnns_in_topk`](@ref AppleAccelerate.bnns_in_topk) | `Float16`, `Float32` |
| uniform / normal / categorical random fills | `Float16`, `Float32` |
| [`bnns_random_fill_uniform_int!`](@ref AppleAccelerate.bnns_random_fill_uniform_int!) | `Int8`–`Int64`, `UInt8`–`UInt64` |
| graph inputs / outputs | whatever the model declares: `Float16`, `Float32`, integers, `Bool` |

```@setup bnns
using AppleAccelerate
```

!!! note "Namespace"
    These functions are not exported. Access them via the `AppleAccelerate.`
    prefix (e.g. `AppleAccelerate.bnns_reduce`).

!!! warning "Deprecated APIs are excluded"
    Apple deprecated the classic BNNS filter/layer API and much of the classic
    tensor/DirectApply surface (macOS 15 / iOS 18) in favour of the newer
    **BNNS Graph API** ([`BNNSGraph`](@ref AppleAccelerate.BNNSGraph)). This
    package does **not** wrap any of those deprecated entry points; use the Graph
    API for that functionality. The "What's left to the raw layer" section below
    lists the full excluded set.

## Descriptors

[`BNNSArray`](@ref AppleAccelerate.BNNSArray) builds a GC-safe
`BNNSNDArrayDescriptor` view of a dense, contiguous Julia array. Internally the
N-D op wrappers map a column-major Julia `Array` onto a `BNNSDataLayout{N}DLastMajor`
descriptor with explicit strides, so **BNNS axis `k` corresponds to Julia
dimension `k+1`** (axis 0 is the contiguous/fastest axis).

```@docs
AppleAccelerate.BNNSArray
```

## Tensor manipulation

Stateless tensor ops that remain current, cross-validated against `permutedims`
and plain copies.

| Function | Meaning |
|----------|---------|
| [`bnns_transpose`](@ref AppleAccelerate.bnns_transpose) | swap two axes |
| [`bnns_copy!`](@ref AppleAccelerate.bnns_copy!) | copy, optionally converting the element type |

```@example bnns
M = Float32[1 2 3; 4 5 6]
@assert AppleAccelerate.bnns_transpose(M, 1, 2) == permutedims(M, (2, 1))
@assert AppleAccelerate.bnns_copy!(zeros(Float32, 2, 3), M) == M
H = AppleAccelerate.bnns_copy!(zeros(Float16, 2, 3), M)     # Float32 -> Float16
@assert H == Float16.(M)
nothing # hide
```

```@docs
AppleAccelerate.bnns_transpose
AppleAccelerate.bnns_copy!
```

## Reductions

```@docs
AppleAccelerate.bnns_reduce
```

## DirectApply kernels

Fused kernels that run without an explicit filter handle.

```@docs
AppleAccelerate.bnns_topk
AppleAccelerate.bnns_in_topk
```

## Utility queries

```@docs
AppleAccelerate.bnns_layout_rank
AppleAccelerate.bnns_data_size
AppleAccelerate.bnns_tensor_allocation_size
```

## Random number generation

[`BNNSRandomGenerator`](@ref AppleAccelerate.BNNSRandomGenerator) is an AES-CTR
generator with an optional seed; the fill functions populate arrays in place and
the state can be snapshot and restored for reproducibility.

```@docs
AppleAccelerate.BNNSRandomGenerator
AppleAccelerate.bnns_random_fill_uniform!
AppleAccelerate.bnns_random_fill_uniform_int!
AppleAccelerate.bnns_random_fill_normal!
AppleAccelerate.bnns_random_fill_categorical!
AppleAccelerate.bnns_random_state
AppleAccelerate.bnns_random_state!
```

## Nearest neighbors

```@docs
AppleAccelerate.BNNSNearestNeighbors
AppleAccelerate.bnns_knn_load!
AppleAccelerate.bnns_knn_query
```

## [BNNS Graph API](@id bnns-graph)

The modern, non-deprecated pipeline (macOS 15+): compile a Core ML model into a
[`BNNSGraph`](@ref AppleAccelerate.BNNSGraph), make an executable
[`BNNSGraphContext`](@ref AppleAccelerate.BNNSGraphContext), look at what it
expects with [`bnns_graph_arguments`](@ref AppleAccelerate.bnns_graph_arguments),
and run it on Julia arrays with
[`bnns_graph_run`](@ref AppleAccelerate.bnns_graph_run) /
[`bnns_graph_run!`](@ref AppleAccelerate.bnns_graph_run!).

The input is a **compiled Core ML model** — the `.mlmodelc` directory that Xcode or
`xcrun coremlcompiler compile model.mlpackage out/` produces (ML Program models
only). There is no in-memory graph builder in this API. A `.mlmodelc` is a
directory holding a textual MIL program, `model.mil`, plus a weights blob; the
example below writes a tiny one by hand so that it is self-contained.

```@example bnns
modeldir = joinpath(mktempdir(), "dense.mlmodelc"); mkpath(modeldir)
write(joinpath(modeldir, "model.mil"), """
program(1.3)
[buildInfo = dict<string, string>({{"coremlc-component-MIL", "handwritten"}})]
{
    func main<ios16>(tensor<fp32, [2, 3]> x) {
            tensor<fp32, [3]> b = const()[name = string("b"), val = tensor<fp32, [3]>([1.0, -2.0, 0.5])];
            tensor<fp32, [2, 3]> s = add(x = x, y = b)[name = string("s")];
            tensor<fp32, [2, 3]> z = relu(x = s)[name = string("z")];
        } -> (z);
}
""")

graph = AppleAccelerate.BNNSGraph(modeldir)
ctx   = AppleAccelerate.BNNSGraphContext(graph)
AppleAccelerate.bnns_graph_arguments(ctx)
```

### Memory layout

MIL tensors are row-major, Julia arrays are column-major, so a model tensor of
shape `[2, 3]` has the memory layout of a Julia array of size `(3, 2)`. The graph
functions pass Julia's memory to BNNS untouched, which means **every argument is a
Julia array whose `size` is the reverse of the model's shape** — that is the
`size` field reported above. For an image model, `[N, C, H, W]` is a Julia
`(W, H, C, N)` array. BNNS graphs ignore custom strides, so this is the only
zero-copy mapping.

```@example bnns
X = Float32[1 2 3; -4 0 6]                       # in the model's [2, 3] index order
out = AppleAccelerate.bnns_graph_run(ctx, "x" => permutedims(X))   # (3, 2): reversed dims
@assert permutedims(out["z"]) == max.(X .+ Float32[1 -2 0.5], 0)

# mil_order=true does that permutedims for you, on the way in and on the way out
out = AppleAccelerate.bnns_graph_run(ctx, "x" => X; mil_order = true)
@assert out["z"] == max.(X .+ Float32[1 -2 0.5], 0)
nothing # hide
```

For repeated inference preallocate the outputs and a page-aligned workspace, and
nothing is allocated per call on the BNNS side:

```@example bnns
Z  = zeros(Float32, 3, 2)
ws = AppleAccelerate.bnns_graph_workspace(ctx)
AppleAccelerate.bnns_graph_run!(ctx, "z" => Z, "x" => permutedims(X); workspace = ws)
@assert permutedims(Z) == out["z"]
nothing # hide
```

A half-precision model takes and returns `Float16` arrays; convert with
[`bnns_copy!`](@ref AppleAccelerate.bnns_copy!) (or plain `Float16.(x)`). Models
with a dynamic batch dimension are bound with
[`bnns_graph_context_set_batch_size!`](@ref AppleAccelerate.bnns_graph_context_set_batch_size!),
more general dynamic shapes with
[`bnns_graph_context_set_dynamic_shapes!`](@ref AppleAccelerate.bnns_graph_context_set_dynamic_shapes!).
A context carries mutable state and must be used by one thread at a time; make
one context per task.

```@docs
AppleAccelerate.BNNSGraphCompileOptions
AppleAccelerate.BNNSGraph
AppleAccelerate.BNNSGraphContext
AppleAccelerate.BNNSGraphArgument
AppleAccelerate.bnns_graph_arguments
AppleAccelerate.bnns_graph_run
AppleAccelerate.bnns_graph_run!
AppleAccelerate.bnns_graph_workspace
AppleAccelerate.bnns_graph_context_workspace_size
AppleAccelerate.bnns_graph_context_set_batch_size!
AppleAccelerate.bnns_graph_context_set_dynamic_shapes!
AppleAccelerate.bnns_graph_argument_names
AppleAccelerate.bnns_graph_execute!
```

The compile-options accessors
(`bnns_compile_options_set_single_thread!`,
`…_set_optimization!`, `…_set_output_path!`, and their getters) and the remaining
graph introspection helpers (`bnns_graph_input_count`, `bnns_graph_input_names`,
`bnns_graph_argument_intents`, `bnns_graph_argument_position`, …) round out the
family.

!!! note "Versioned symbols"
    `bnns_graph.h` redirects `BNNSGraphCompileFromFile`, `BNNSGraphContextExecute`
    and eight other functions to `_v2` symbols with an `__asm__` label. The
    un-suffixed symbols that libBNNS still exports have a different, pre-release
    argument list, so the generated `LibAccelerate.BNNSGraphCompileFromFile` &c.
    must **not** be called directly (they crash); the wrappers on this page bind
    the `_v2` symbols.

## What's left to the raw layer

Everything not wrapped above is reachable through the raw
`AppleAccelerate.LibAccelerate` layer. It falls into three groups:

- **Deprecated classic tensor / DirectApply kernels** (macOS 15 / iOS 18) —
  `BNNSMatMul`, the activation-filter path, `BNNSTile`/`BNNSTileBackward`,
  `BNNSCompareTensor`, `BNNSBandPart`, `BNNSGather`/`BNNSScatter` (and their ND
  forms), `BNNSShuffle`, the clip family (`BNNSClipByValue`/`…ByNorm`/
  `…ByGlobalNorm`), `BNNSComputeNorm`, `BNNSOptimizerStep`, and the
  `BNNSDirectApply{ActivationBatch,BroadcastMatMul,Quantizer}` kernels.
  Superseded by the BNNS Graph API; intentionally not given an idiomatic wrapper.
- **The deprecated classic filter/layer API** — the `BNNSFilterCreate*` /
  `BNNSFilterCreateLayer*` constructors and their `*FilterApply*` /
  `BNNSFusedFilterApply*` execute paths (including the two-input / fused / loss /
  normalization / pooling / permute batch variants).
- **Exotic, training-only entry points** that need training caches or opaque
  multi-kilobyte parameter blocks that cannot be validated generically: multi-head
  attention (`BNNSApplyMultiheadAttention` and its backward), the LSTM
  training-cache path (`BNNSComputeLSTMTrainingCacheCapacity`,
  `BNNSDirectApplyLSTMBatchTrainingCaching` / `…Backward`),
  `BNNSComputeNormBackward`, image crop/resize (`BNNSCropResize` /
  `BNNSCropResizeBackward`), and the fully-connected sparsification helpers
  (`BNNSNDArrayFullyConnectedSparsifySparse{COO,CSR}`).
