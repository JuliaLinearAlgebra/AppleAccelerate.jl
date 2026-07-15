# Neural Network Primitives (BNNS)

AppleAccelerate wraps a small, numerically verified subset of Apple's
[BNNS (Basic Neural Network Subroutines)](https://developer.apple.com/documentation/accelerate/bnns)
library — the primitives that are broadly useful and that can be checked against a
plain-Julia reference. These wrappers are **`Float32`-centric**, matching BNNS's
native precision for inference.

```@setup bnns
using AppleAccelerate
```

!!! note "Namespace"
    Like the rest of the package, these functions are not exported. Access them via
    the `AppleAccelerate.` prefix (e.g. `AppleAccelerate.bnns_matmul`).

## Descriptors

[`BNNSArray`](@ref AppleAccelerate.BNNSArray) builds a GC-safe
`BNNSNDArrayDescriptor` view of a dense, contiguous Julia array. It keeps the
backing array alive and reports a column-major-friendly BNNS layout (1D vectors as
`BNNSDataLayoutVector`, 2D matrices as `BNNSDataLayoutColumnMajorMatrix`). Most
users do not need to construct one directly — the matmul and activation helpers
build descriptors internally — but it is the building block if you reach into the
raw layer.

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

## Activations

[`bnns_activation`](@ref AppleAccelerate.bnns_activation) applies a pointwise
activation function to every element of a `Float32` array using the BNNS
activation-layer filter API, returning a new array;
[`bnns_activation!`](@ref AppleAccelerate.bnns_activation!) writes into a
preallocated output of the same shape. Supported functions are `:identity`,
`:relu`, `:sigmoid`, `:tanh`, and `:abs`.

```@example bnns
x = Float32[-2, -1, 0, 1, 2]

AppleAccelerate.bnns_activation(:relu, x)     # max.(x, 0)
AppleAccelerate.bnns_activation(:abs, x)      # abs.(x)
y = AppleAccelerate.bnns_activation(:sigmoid, x)
@assert y ≈ 1f0 ./ (1f0 .+ exp.(-x))
nothing # hide
```

```@docs
AppleAccelerate.BNNSArray
AppleAccelerate.bnns_matmul
AppleAccelerate.bnns_matmul!
AppleAccelerate.bnns_activation
AppleAccelerate.bnns_activation!
```

## Using BNNS with NNlib

AppleAccelerate defines no methods on NNlib's functions — loading it never changes
what `batched_mul!` or an activation function does. To use BNNS for a batched
`Float32` matmul, call [`bnns_matmul!`](@ref AppleAccelerate.bnns_matmul!) per
slice explicitly:

```@example bnns
A = rand(Float32, 4, 3, 8)   # 8 batched 4x3 matrices
B = rand(Float32, 3, 5, 8)   # 8 batched 3x5 matrices
C = Array{Float32}(undef, 4, 5, 8)

for i in axes(C, 3)
    AppleAccelerate.bnns_matmul!(view(C, :, :, i), view(A, :, :, i), view(B, :, :, i))
end
@assert C[:, :, 1] ≈ A[:, :, 1] * B[:, :, 1]
nothing # hide
```

For activations, [`bnns_activation`](@ref AppleAccelerate.bnns_activation) takes a
`Symbol` (`:relu`, `:sigmoid`, `:tanh`, `:abs`, `:identity`) rather than an NNlib
function object, so no NNlib dependency is involved:

```@example bnns
x = randn(Float32, 6)
@assert AppleAccelerate.bnns_activation(:relu, x) ≈ max.(x, 0f0)
nothing # hide
```

## What's left to the raw layer

BNNS is large (~139 C functions, 110+ structs). Deliberately *not* given idiomatic
wrappers — use the raw `AppleAccelerate.LibAccelerate` layer
directly if you need them (see [Architecture](@ref architecture)) — are the BNNS graph-compiler API (`bnns_graph.h`),
convolution / pooling / normalization / LSTM / attention layers, quantization, and
training. These require substantial stateful plumbing or are hard to verify
generically, so they fall outside this idiomatic core.
