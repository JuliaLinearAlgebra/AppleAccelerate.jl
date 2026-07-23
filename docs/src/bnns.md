# Neural Network Primitives (BNNS)

AppleAccelerate wraps the **current, non-deprecated** slice of Apple's
[BNNS (Basic Neural Network Subroutines)](https://developer.apple.com/documentation/accelerate/bnns)
library — 61 of the ~136 `BNNS*` C entry points. The bulk of the remainder are
APIs Apple **deprecated in macOS 15 / iOS 18**: the classic filter/layer
construction API and the deprecated classic + DirectApply tensor kernels
(`BNNSMatMul`, `BNNSTile`, `BNNSGather`/`BNNSScatter`, the clip / norm family,
`BNNSOptimizerStep`, …). Those are intentionally **not wrapped** — target the
**BNNS Graph API** instead. The wrappers are **`Float32`-centric**, matching
BNNS's native precision for inference. Numerically verified helpers (transpose,
copy, reductions, top-k, random generation, nearest neighbors) are cross-checked
here against plain-Julia references; the remaining thin wrappers expose the rest
of the current surface with exact FFI signatures for callers who need them.

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
| [`bnns_copy!`](@ref AppleAccelerate.bnns_copy!) | copy with BNNS layout rules |

```@example bnns
M = Float32[1 2 3; 4 5 6]
@assert AppleAccelerate.bnns_transpose(M, 1, 2) == permutedims(M, (2, 1))
@assert AppleAccelerate.bnns_copy!(zeros(Float32, 2, 3), M) == M
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

## BNNS Graph API

The modern, non-deprecated pipeline: build compile options, compile a serialized
graph package into a [`BNNSGraph`](@ref AppleAccelerate.BNNSGraph), make an
executable [`BNNSGraphContext`](@ref AppleAccelerate.BNNSGraphContext), then
introspect and execute. Compiling requires an on-disk graph package (there is no
in-memory graph builder in this API).

```@docs
AppleAccelerate.BNNSGraphCompileOptions
AppleAccelerate.BNNSGraph
AppleAccelerate.BNNSGraphContext
AppleAccelerate.bnns_graph_execute!
AppleAccelerate.bnns_graph_context_workspace_size
```

The compile-options accessors
(`bnns_compile_options_set_single_thread!`,
`…_set_optimization!`, `…_set_output_path!`, and their getters) and the graph
introspection helpers (`bnns_graph_input_count`, `bnns_graph_argument_names`,
`bnns_graph_argument_intents`, …) round out the family.

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
