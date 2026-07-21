# BNNS (Basic Neural Network Subroutines) idiomatic wrappers.
#
# This module wraps a small, coherent, *numerically verified* subset of Apple's
# BNNS API exposed by the raw `LibAccelerate` layer. BNNS is large (~139 C
# functions, 110+ structs, including a full graph-compiler API and many exotic
# layer types). Rather than mechanically wrapping everything, we expose the
# primitives that are broadly useful and that we can prove correct against a
# plain-Julia reference:
#
#   * `BNNSArray`            — a GC-safe `BNNSNDArrayDescriptor` built from a Julia array.
#   * `bnns_matmul` / `!`    — matrix multiply via `BNNSMatMul`.
#   * `bnns_activation` / `!`— pointwise activation via the activation-layer filter API.
#
# Intentionally left to the raw `LibAccelerate` layer (use it directly if you
# need these): the BNNS graph API (`bnns_graph.h`), convolution / pooling /
# normalization / LSTM / attention layers, quantization, and training. These
# either require substantial stateful plumbing or are hard to verify generically;
# they are out of scope for this idiomatic core.
#
# NOTE: The BNNS filter/matmul/activation and descriptor API wrapped here is
# DEPRECATED by Apple as of macOS 15.0 / iOS 18.0 in favor of the newer BNNS
# Graph API (`bnns_graph.h`, left to the raw `LibAccelerate` layer). These
# wrappers still function but emit C-level deprecation warnings and may be
# removed in a future macOS; new code should target the BNNS Graph API.
#
# Apple docs: https://developer.apple.com/documentation/accelerate/bnns

using .LibAccelerate:
    BNNSNDArrayDescriptor, BNNSDataType, BNNSDataLayout,
    BNNSDataTypeFloat32, BNNSDataTypeInt32,
    BNNSDataLayoutVector, BNNSDataLayoutColumnMajorMatrix,
    BNNSActivation, BNNSActivationFunction,
    BNNSActivationFunctionIdentity, BNNSActivationFunctionRectifiedLinear,
    BNNSActivationFunctionSigmoid, BNNSActivationFunctionTanh,
    BNNSActivationFunctionAbs,
    BNNSLayerParametersActivation, BNNSFilterParameters,
    BNNSFilterCreateLayerActivation, BNNSFilterApply, BNNSFilterDestroy,
    BNNSMatMul, BNNSMatMulWorkspaceSize

# --- Data-type mapping -------------------------------------------------------

# Map a Julia element type to the corresponding BNNS data-type enum value.
# The struct fields are typed `BNNSDataType` (a `UInt32` alias), so the enum
# value is converted to its underlying integer code.
bnns_data_type(::Type{Float32}) = BNNSDataType(BNNSDataTypeFloat32)
bnns_data_type(::Type{Int32})   = BNNSDataType(BNNSDataTypeInt32)
bnns_data_type(::Type{T}) where {T} =
    throw(ArgumentError("BNNS: unsupported element type $T (supported: Float32, Int32)"))

const BNNSFloat = Float32  # element type supported by these wrappers

# --- Descriptor construction -------------------------------------------------

"""
    BNNSArray(A::AbstractArray)

A GC-safe `BNNSNDArrayDescriptor` view of a Julia array `A`, suitable for passing
to BNNS routines via `Ref`. The wrapper keeps a reference to the backing array so
it is not collected while the descriptor is alive; pass the underlying descriptor
with `Base.cconvert`/`Ref` only inside a `GC.@preserve` block guarding `A`.

BNNS descriptors are *layout aware*. Julia stores arrays in column-major order, so
this constructor reports the array using a column-major-friendly BNNS layout:

  * 1D `Vector`  -> `BNNSDataLayoutVector`
  * 2D `Matrix`  -> `BNNSDataLayoutColumnMajorMatrix` (BNNS `size = (rows, cols)`).

Only `Float32` (and `Int32`) dense, contiguous arrays are supported here; other
element types or strided/transposed arrays should use the raw `LibAccelerate`
layer directly.

See also [`bnns_matmul`](@ref), [`bnns_activation`](@ref).

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
struct BNNSArray{T,N}
    desc::BNNSNDArrayDescriptor
    parent::Array{T,N}    # keep alive; also asserts contiguity
end

function BNNSArray(A::Array{T}) where {T}
    dt = bnns_data_type(T)
    n = ndims(A)
    sz = ntuple(i -> i <= n ? Csize_t(size(A, i)) : Csize_t(0), 8)
    # Contiguous column-major strides; BNNS uses 0 to mean "default/contiguous"
    # but we set them explicitly for clarity.
    st = ntuple(8) do i
        if i > n
            Csize_t(0)
        elseif i == 1
            Csize_t(1)
        else
            Csize_t(prod(size(A)[1:i-1]))
        end
    end
    layout = if n == 1
        BNNSDataLayoutVector
    elseif n == 2
        BNNSDataLayoutColumnMajorMatrix
    else
        throw(ArgumentError("BNNSArray supports 1D and 2D arrays; got $(n)D"))
    end
    desc = BNNSNDArrayDescriptor(
        UInt32(0),                 # flags
        BNNSDataLayout(layout),    # layout
        sz,                        # size
        st,                        # stride
        Ptr{Cvoid}(pointer(A)),    # data
        dt,                        # data_type
        C_NULL,                    # table_data
        dt,                        # table_data_type (unused)
        1.0f0,                     # data_scale
        0.0f0,                     # data_bias
    )
    return BNNSArray{T,ndims(A)}(desc, A)
end

# --- Matrix multiply ---------------------------------------------------------

"""
    bnns_matmul(A::Matrix{Float32}, B::Matrix{Float32}; alpha = 1.0f0) -> Matrix{Float32}

Compute `alpha * (A * B)` using BNNS (`BNNSMatMul`). `A` is `m×k`, `B` is `k×n`,
and the result is `m×n`. Equivalent to `alpha .* (A * B)`.

See also [`bnns_matmul!`](@ref).

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
function bnns_matmul(A::AbstractMatrix{BNNSFloat}, B::AbstractMatrix{BNNSFloat};
                     alpha::Real = 1.0f0)
    m = size(A, 1)
    n = size(B, 2)
    C = Matrix{BNNSFloat}(undef, m, n)
    return bnns_matmul!(C, A, B; alpha)
end

"""
    bnns_matmul!(C, A, B; alpha = 1.0f0) -> C

In-place matrix multiply: `C .= alpha .* (A * B)`. All arguments are `Float32`
matrices with conforming dimensions. Returns `C`.

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
function bnns_matmul!(C::AbstractMatrix{BNNSFloat}, A::AbstractMatrix{BNNSFloat},
                      B::AbstractMatrix{BNNSFloat}; alpha::Real = 1.0f0)
    m, k = size(A)
    k2, n = size(B)
    k == k2 || throw(DimensionMismatch("A is $(size(A)), B is $(size(B))"))
    size(C) == (m, n) || throw(DimensionMismatch("C is $(size(C)), expected ($m, $n)"))

    Ac = A isa Array ? A : Array(A)
    Bc = B isa Array ? B : Array(B)
    # BNNSArray requires a dense `Array`; if C is a view (e.g. a batched slice),
    # compute into a temporary and copy the result back at the end.
    Cc = C isa Array{BNNSFloat,2} ? C : Matrix{BNNSFloat}(undef, m, n)

    da = BNNSArray(Ac)
    db = BNNSArray(Bc)
    dc = BNNSArray(Cc)
    fp = BNNSFilterParameters(UInt32(0), Csize_t(0), C_NULL, C_NULL)

    GC.@preserve Ac Bc Cc da db dc begin
        ra = Ref(da.desc); rb = Ref(db.desc); rc = Ref(dc.desc); rfp = Ref(fp)
        # BNNSMatMul computes out = alpha * inputA * inputB (no transpose).
        ws = BNNSMatMulWorkspaceSize(false, false, Float32(alpha), ra, rb, rc, rfp)
        workspace = ws > 0 ? Vector{UInt8}(undef, Int(ws)) : UInt8[]
        GC.@preserve workspace begin
            wptr = isempty(workspace) ? C_NULL : pointer(workspace)
            status = BNNSMatMul(false, false, Float32(alpha), ra, rb, rc, wptr, rfp)
            status == 0 || error("BNNSMatMul failed with status $status")
        end
    end
    Cc === C || copyto!(C, Cc)
    return C
end

# --- Pointwise activation ----------------------------------------------------

"""
    bnns_activation(f::Symbol, X::AbstractArray{Float32}; alpha=0.0f0, beta=0.0f0) -> AbstractArray{Float32}

Apply a pointwise activation function `f` to every element of `X` using the BNNS
activation-layer filter API, returning a new array. Supported `f`:

  * `:identity`
  * `:relu`     — rectified linear
  * `:sigmoid`
  * `:tanh`
  * `:abs`

`alpha`/`beta` are passed to BNNS for activations that take parameters (e.g. leaky
variants); they are ignored by the activations listed above.

See also [`bnns_activation!`](@ref).

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
function bnns_activation(f::Symbol, X::AbstractArray{BNNSFloat};
                         alpha::Real = 0.0f0, beta::Real = 0.0f0)
    Y = similar(Array(X))
    return bnns_activation!(f, Y, X; alpha, beta)
end

const _BNNS_ACT = Dict(
    :identity => BNNSActivationFunctionIdentity,
    :relu     => BNNSActivationFunctionRectifiedLinear,
    :sigmoid  => BNNSActivationFunctionSigmoid,
    :tanh     => BNNSActivationFunctionTanh,
    :abs      => BNNSActivationFunctionAbs,
)

"""
    bnns_activation!(f::Symbol, Y::AbstractArray{Float32}, X::AbstractArray{Float32}; alpha=0.0f0, beta=0.0f0) -> Y

In-place pointwise activation: `Y .= f.(X)`. See [`bnns_activation`](@ref) for the
list of supported `f`. `X` and `Y` must have the same shape. Returns `Y`.

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
function bnns_activation!(f::Symbol, Y::AbstractArray{BNNSFloat}, X::AbstractArray{BNNSFloat};
                          alpha::Real = 0.0f0, beta::Real = 0.0f0)
    size(Y) == size(X) || throw(DimensionMismatch("Y is $(size(Y)), X is $(size(X))"))
    haskey(_BNNS_ACT, f) ||
        throw(ArgumentError("BNNS: unsupported activation $(repr(f)); supported: $(sort(collect(keys(_BNNS_ACT))))"))

    # Activation is fully pointwise; treat the data as a flat vector so the
    # descriptor layout is trivial and identical for any array shape.
    Xc = X isa Vector ? X : vec(Array(X))
    Yc = (Y isa Vector) ? Y : Vector{BNNSFloat}(undef, length(Y))

    di = BNNSArray(Xc)
    do_ = BNNSArray(Yc)
    act = BNNSActivation(BNNSActivationFunction(_BNNS_ACT[f]),
                         Float32(alpha), Float32(beta),
                         Int32(0), Int32(0), Int32(0),
                         C_NULL, C_NULL, C_NULL)
    fp = BNNSFilterParameters(UInt32(0), Csize_t(0), C_NULL, C_NULL)

    GC.@preserve Xc Yc di do_ begin
        lp = BNNSLayerParametersActivation(di.desc, do_.desc, act, UInt32(0))
        filter = GC.@preserve lp fp begin
            BNNSFilterCreateLayerActivation(Ref(lp), Ref(fp))
        end
        filter == C_NULL && error("BNNSFilterCreateLayerActivation returned NULL")
        try
            status = BNNSFilterApply(filter, pointer(Xc), pointer(Yc))
            status == 0 || error("BNNSFilterApply failed with status $status")
        finally
            BNNSFilterDestroy(filter)
        end
    end

    if !(Y isa Vector)
        copyto!(Y, reshape(Yc, size(Y)))
    end
    return Y
end
