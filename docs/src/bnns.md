# Neural Network Primitives (BNNS)

AppleAccelerate wraps a broad slice of Apple's
[BNNS (Basic Neural Network Subroutines)](https://developer.apple.com/documentation/accelerate/bnns)
library — 83 of the ~136 `BNNS*` C entry points. Most of what is left unwrapped is
the **deprecated** classic filter/layer construction API (superseded by the Graph
API, see the deprecation note below); of the modern, non-deprecated surface roughly
90% is covered. The wrappers are **`Float32`-centric**, matching BNNS's native
precision for inference. Numerically verified helpers (matrix multiply, activations,
tensor ops, reductions, clipping, random generation, nearest neighbors, the
optimizer step) are cross-checked here against plain-Julia references; the remaining
thin wrappers expose the rest of the surface with exact FFI signatures for callers
who need them.

```@setup bnns
using AppleAccelerate
```

!!! note "Namespace"
    These functions are not exported. Access them via the `AppleAccelerate.`
    prefix (e.g. `AppleAccelerate.bnns_matmul`).

!!! warning "Deprecation"
    Apple deprecated the classic BNNS filter/layer API (macOS 15 / iOS 18) in
    favour of the newer **BNNS Graph API** ([`BNNSGraph`](@ref AppleAccelerate.BNNSGraph)).
    The layer wrappers below still work but emit C-level deprecation warnings;
    new code should target the Graph API.

## Descriptors

[`BNNSArray`](@ref AppleAccelerate.BNNSArray) builds a GC-safe
`BNNSNDArrayDescriptor` view of a dense, contiguous Julia array. Internally the
N-D op wrappers map a column-major Julia `Array` onto a `BNNSDataLayout{N}DLastMajor`
descriptor with explicit strides, so **BNNS axis `k` corresponds to Julia
dimension `k+1`** (axis 0 is the contiguous/fastest axis).

```@docs
AppleAccelerate.BNNSArray
```

## Matrix multiply

[`bnns_matmul`](@ref AppleAccelerate.bnns_matmul) computes `alpha * (A * B)` via
`BNNSMatMul`. `A` is `m×k`, `B` is `k×n`, and the result is `m×n`. The in-place
variant [`bnns_matmul!`](@ref AppleAccelerate.bnns_matmul!) writes into a
preallocated `C`.

```@example bnns
A = rand(Float32, 4, 3)
B = rand(Float32, 3, 5)

C = AppleAccelerate.bnns_matmul(A, B)          # == A * B
C2 = AppleAccelerate.bnns_matmul(A, B; alpha = 2.0f0)  # == 2 .* (A * B)

@assert C  ≈ A * B
@assert C2 ≈ 2 .* (A * B)
nothing # hide
```

```@docs
AppleAccelerate.bnns_matmul
AppleAccelerate.bnns_matmul!
```

## Activations

[`bnns_activation`](@ref AppleAccelerate.bnns_activation) applies a pointwise
activation to a `Float32` array, returning a new array;
[`bnns_activation!`](@ref AppleAccelerate.bnns_activation!) writes into a
preallocated output. Supported functions are `:identity`, `:relu`, `:sigmoid`,
`:tanh`, and `:abs`.

```@example bnns
x = Float32[-2, -1, 0, 1, 2]
y = AppleAccelerate.bnns_activation(:sigmoid, x)
@assert y ≈ 1f0 ./ (1f0 .+ exp.(-x))
nothing # hide
```

```@docs
AppleAccelerate.bnns_activation
AppleAccelerate.bnns_activation!
```

## Tensor manipulation

Stateless tensor ops, cross-validated against `repeat`, `permutedims`, `tril`/
`triu`, indexing and broadcast comparisons.

| Function | Meaning |
|----------|---------|
| [`bnns_tile`](@ref AppleAccelerate.bnns_tile) / [`bnns_tile_backward`](@ref AppleAccelerate.bnns_tile_backward) | `repeat` and its adjoint |
| [`bnns_transpose`](@ref AppleAccelerate.bnns_transpose) | swap two axes |
| [`bnns_compare`](@ref AppleAccelerate.bnns_compare) | element-wise relational op → `Bool` |
| [`bnns_band_part`](@ref AppleAccelerate.bnns_band_part) | keep a diagonal band (`tril`/`triu`) |
| [`bnns_gather`](@ref AppleAccelerate.bnns_gather) / [`bnns_scatter`](@ref AppleAccelerate.bnns_scatter) | gather/scatter along an axis |
| [`bnns_gather_nd`](@ref AppleAccelerate.bnns_gather_nd) / [`bnns_scatter_nd`](@ref AppleAccelerate.bnns_scatter_nd) | N-D coordinate gather/scatter |
| [`bnns_shuffle`](@ref AppleAccelerate.bnns_shuffle) | pixel / depth shuffle (NCHW) |
| [`bnns_copy!`](@ref AppleAccelerate.bnns_copy!) | copy with BNNS layout rules |

```@example bnns
M = Float32[1 2 3; 4 5 6]
@assert AppleAccelerate.bnns_tile(M, (2, 1)) == repeat(M, 2, 1)
@assert AppleAccelerate.bnns_band_part(M[:, 1:2], -1, 0) == [1f0 0; 4 5]  # tril
nothing # hide
```

```@docs
AppleAccelerate.bnns_tile
AppleAccelerate.bnns_tile_backward
AppleAccelerate.bnns_transpose
AppleAccelerate.bnns_compare
AppleAccelerate.bnns_band_part
AppleAccelerate.bnns_gather
AppleAccelerate.bnns_scatter
AppleAccelerate.bnns_gather_nd
AppleAccelerate.bnns_scatter_nd
AppleAccelerate.bnns_shuffle
AppleAccelerate.bnns_copy!
```

## Reductions, norms and clipping

```@docs
AppleAccelerate.bnns_reduce
AppleAccelerate.bnns_compute_norm
AppleAccelerate.bnns_clip_by_value
AppleAccelerate.bnns_clip_by_norm
AppleAccelerate.bnns_clip_by_global_norm
```

## DirectApply kernels

Fused kernels that run without an explicit filter handle.

```@docs
AppleAccelerate.bnns_broadcast_matmul
AppleAccelerate.bnns_topk
AppleAccelerate.bnns_in_topk
AppleAccelerate.bnns_activation_batch!
AppleAccelerate.bnns_quantizer!
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

## Optimizer

```@docs
AppleAccelerate.bnns_optimizer_step_sgd!
```

## What's left to the raw layer

Two groups of entry points are left to the raw `AppleAccelerate.LibAccelerate`
layer:

- **The deprecated classic filter/layer API** — the `BNNSFilterCreate*` /
  `BNNSFilterCreateLayer*` constructors and their `*FilterApply*` /
  `BNNSFusedFilterApply*` execute paths (including the two-input / fused / loss /
  normalization / pooling / permute batch variants). This is the largest unwrapped
  block; it is superseded by the BNNS Graph API (above) and is intentionally
  not given an idiomatic wrapper.
- **Exotic, training-only entry points** that need training caches or opaque
  multi-kilobyte parameter blocks that cannot be validated generically: multi-head
  attention (`BNNSApplyMultiheadAttention` and its backward), the LSTM
  training-cache path (`BNNSComputeLSTMTrainingCacheCapacity`,
  `BNNSDirectApplyLSTMBatchTrainingCaching` / `…Backward`),
  `BNNSComputeNormBackward`, image crop/resize (`BNNSCropResize` /
  `BNNSCropResizeBackward`), and the fully-connected sparsification helpers
  (`BNNSNDArrayFullyConnectedSparsifySparse{COO,CSR}`).
```
