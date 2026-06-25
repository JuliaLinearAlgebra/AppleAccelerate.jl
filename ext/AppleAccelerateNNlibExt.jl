module AppleAccelerateNNlibExt

@static if Sys.isapple()

using NNlib
using AppleAccelerate
using AppleAccelerate: bnns_matmul!, bnns_activation!

# This extension provides Accelerate (BNNS)-backed methods for a small, verified
# subset of NNlib operations:
#
#   * `NNlib.batched_mul!` for dense `Float32` batches — each `m×k × k×n` slice is
#     evaluated with `BNNSMatMul` (no transpose, β == 0 only).
#   * `bnns_act` — a convenience that applies a BNNS activation function over an
#     `Array{Float32}`, matching the corresponding NNlib scalar activation.
#
# We only intercept the cases BNNS handles cleanly and that we test numerically;
# everything else falls through to NNlib's generic implementation. Convolution,
# pooling, etc. are left to NNlib / the raw `LibAccelerate` layer.

# --- batched_mul -------------------------------------------------------------

# 5-argument batched_mul! is NNlib's core mutating entry point: C = α·A⊠B + β·C.
# BNNSMatMul computes C = α·A·B with no accumulation, so we only take the β == 0
# case with contiguous (non-transposed) Float32 arrays; anything else defers to
# NNlib via `@invoke`.
function NNlib.batched_mul!(C::Array{Float32,3}, A::Array{Float32,3}, B::Array{Float32,3},
                            α::Number, β::Number)
    if iszero(β) && size(A, 3) == size(B, 3) == size(C, 3)
        @inbounds for i in axes(C, 3)
            bnns_matmul!(view(C, :, :, i), view(A, :, :, i), view(B, :, :, i);
                         alpha = Float32(α))
        end
        return C
    end
    # Fall back to NNlib's generic method for the unsupported cases.
    return @invoke NNlib.batched_mul!(C::AbstractArray{Float32,3},
                                      A::AbstractArray{<:Any,3},
                                      B::AbstractArray{<:Any,3},
                                      α::Number, β::Number)
end

# --- activations -------------------------------------------------------------

# Map NNlib activation functions to BNNS activation symbols understood by
# `bnns_activation!`.
const _NNLIB_ACT = IdDict{Any,Symbol}(
    NNlib.relu      => :relu,
    NNlib.sigmoid   => :sigmoid,
    NNlib.σ         => :sigmoid,
    tanh            => :tanh,
    NNlib.tanh_fast => :tanh,
    abs             => :abs,
    identity        => :identity,
)

"""
    AppleAccelerate.bnns_act(f, X::Array{Float32}) -> Array{Float32}

Apply the NNlib activation function `f` (one of `relu`, `sigmoid`/`σ`, `tanh`,
`tanh_fast`, `abs`, `identity`) to `X` using BNNS. Available once both
`AppleAccelerate` and `NNlib` are loaded. Falls back to `f.(X)` for unsupported
activations.
"""
function AppleAccelerate.bnns_act(f, X::Array{Float32})
    sym = get(_NNLIB_ACT, f, nothing)
    sym === nothing && return f.(X)
    Y = similar(X)
    return bnns_activation!(sym, Y, X)
end

end # @static if Sys.isapple()

end # module
