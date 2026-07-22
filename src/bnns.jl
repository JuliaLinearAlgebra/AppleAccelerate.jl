# BNNS (Basic Neural Network Subroutines) idiomatic wrappers.
#
# This module wraps a broad slice (~90%) of Apple's BNNS API exposed by the raw
# `LibAccelerate` layer — tensor ops, reductions/clipping, random generation,
# nearest neighbors, the whole BNNS Graph compiler pipeline, the classic layer
# filters, and the optimizer step. The numerically important helpers are
# cross-validated against plain-Julia references; the deprecated layer-create /
# gradient wrappers are thin passthroughs over the raw parameter structs.
#
# A handful of exotic, training-only entry points are left to the raw layer:
# multi-head attention (forward+backward), the LSTM training-cache path, the
# two-input/fused/loss backward variants, and the fully-connected sparsification
# helpers. They need training caches or opaque multi-kilobyte parameter blocks
# that cannot be validated generically.
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

# Alias the raw layer so the many BNNS enum values, structs and ccall wrappers can
# be reached without an unwieldy explicit import list. All new wrappers below go
# through `LA.<name>`.
import .LibAccelerate as LA

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

# =============================================================================
# General N-D descriptor plumbing shared by the tensor-op / random / graph
# wrappers below. These cover a much wider slice of the BNNS surface than the
# `BNNSArray` core above.
#
# Layout convention (verified empirically, Jul 2026): BNNS honours the explicit
# per-axis `stride` array, so a column-major Julia `Array` maps 1:1 onto a
# `BNNSDataLayout{N}DLastMajor` descriptor with `size[i] = size(A,i+1)` and
# `stride[i] = stride(A,i+1)` (axis 0 is the contiguous/fastest axis, exactly
# like Julia's first dimension). BNNS axis `k` therefore corresponds to Julia
# dimension `k+1`. `{N}DFirstMajor` is the row-major transpose and is *not* used
# here. Confirmed against `repeat`/`BNNSTile`: only LastMajor returns status 0.
# =============================================================================

# Map a Julia element type to the BNNS data-type enum. Broader than
# `bnns_data_type` (which only accepts the Float32/Int32 core) so index and
# boolean tensors work.
_bnns_dt(::Type{Float32}) = BNNSDataType(LA.BNNSDataTypeFloat32)
_bnns_dt(::Type{Int32})   = BNNSDataType(LA.BNNSDataTypeInt32)
_bnns_dt(::Type{Int64})   = BNNSDataType(LA.BNNSDataTypeInt64)
_bnns_dt(::Type{Int16})   = BNNSDataType(LA.BNNSDataTypeInt16)
_bnns_dt(::Type{Int8})    = BNNSDataType(LA.BNNSDataTypeInt8)
_bnns_dt(::Type{UInt32})  = BNNSDataType(LA.BNNSDataTypeUInt32)
_bnns_dt(::Type{UInt8})   = BNNSDataType(LA.BNNSDataTypeUInt8)
_bnns_dt(::Type{Bool})    = BNNSDataType(LA.BNNSDataTypeBoolean)
_bnns_dt(::Type{T}) where {T} =
    throw(ArgumentError("BNNS: unsupported element type $T"))

const _ND_LAST_MAJOR = (
    LA.BNNSDataLayout1DLastMajor, LA.BNNSDataLayout2DLastMajor,
    LA.BNNSDataLayout3DLastMajor, LA.BNNSDataLayout4DLastMajor,
    LA.BNNSDataLayout5DLastMajor, LA.BNNSDataLayout6DLastMajor,
    LA.BNNSDataLayout7DLastMajor, LA.BNNSDataLayout8DLastMajor,
)

_bnns_layout(n::Integer) =
    1 <= n <= 8 ? BNNSDataLayout(_ND_LAST_MAJOR[n]) :
        throw(ArgumentError("BNNS: rank $n out of range 1:8"))

# Build a `BNNSNDArrayDescriptor` viewing a dense Julia `Array`. The caller is
# responsible for `GC.@preserve`-ing `A` across every ccall that reads the
# returned descriptor's `data` pointer.
function _desc(A::Array{T}) where {T}
    n = ndims(A)
    sz = ntuple(i -> i <= n ? Csize_t(size(A, i)) : Csize_t(0), 8)
    st = ntuple(i -> i <= n ? Csize_t(stride(A, i)) : Csize_t(0), 8)
    return BNNSNDArrayDescriptor(
        LA.BNNSNDArrayFlags(0), _bnns_layout(n), sz, st,
        Ptr{Cvoid}(pointer(A)), _bnns_dt(T), C_NULL, _bnns_dt(T), 1.0f0, 0.0f0)
end

_bnns_check(status, name) =
    status == 0 || error("$name failed with status $status")

# =============================================================================
# Tensor manipulation ops (stateless kernels)
# =============================================================================

"""
    bnns_tile(A::Array, reps) -> Array

Tile `A` by integer repeats `reps` (a tuple/vector, one entry per dimension),
equivalent to `repeat(A, reps...)`, via `BNNSTile`. `Float32`/integer arrays are
supported.

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
function bnns_tile(A::Array{T}, reps) where {T}
    r = ntuple(i -> i <= ndims(A) ? Int(reps[i]) : 1, ndims(A))
    O = Array{T}(undef, ntuple(i -> size(A, i) * r[i], ndims(A)))
    di = _desc(A); do_ = _desc(O)
    GC.@preserve A O begin
        _bnns_check(LA.BNNSTile(Ref(di), Ref(do_), C_NULL), "BNNSTile")
    end
    return O
end

"""
    bnns_tile_backward(out_delta::Array, insize) -> Array

Reduce a tiled gradient `out_delta` back to the untiled shape `insize` by summing
the contributions from each tile, via `BNNSTileBackward`. This is the adjoint of
[`bnns_tile`](@ref).
"""
function bnns_tile_backward(out_delta::Array{T}, insize) where {T}
    isz = ntuple(i -> Int(insize[i]), ndims(out_delta))
    ind = Array{T}(undef, isz)
    di = _desc(ind); do_ = _desc(out_delta)
    GC.@preserve ind out_delta begin
        _bnns_check(LA.BNNSTileBackward(Ref(di), Ref(do_), C_NULL), "BNNSTileBackward")
    end
    return ind
end

"""
    bnns_transpose(A::Array, dim0, dim1) -> Array

Swap Julia dimensions `dim0` and `dim1` of `A` (1-based) via `BNNSTranspose`,
equivalent to a `permutedims` that exchanges those two axes.
"""
function bnns_transpose(A::Array{T}, dim0::Integer, dim1::Integer) where {T}
    p = collect(1:ndims(A)); p[dim0], p[dim1] = p[dim1], p[dim0]
    O = Array{T}(undef, ntuple(i -> size(A, p[i]), ndims(A)))
    di = _desc(A); do_ = _desc(O)
    GC.@preserve A O begin
        _bnns_check(LA.BNNSTranspose(Ref(do_), Ref(di),
                    Csize_t(dim0 - 1), Csize_t(dim1 - 1), C_NULL), "BNNSTranspose")
    end
    return O
end

const _BNNS_RELOP = Dict(
    :(==) => LA.BNNSRelationalOperatorEqual, :eq => LA.BNNSRelationalOperatorEqual,
    :< => LA.BNNSRelationalOperatorLess, :lt => LA.BNNSRelationalOperatorLess,
    :<= => LA.BNNSRelationalOperatorLessEqual, :le => LA.BNNSRelationalOperatorLessEqual,
    :> => LA.BNNSRelationalOperatorGreater, :gt => LA.BNNSRelationalOperatorGreater,
    :>= => LA.BNNSRelationalOperatorGreaterEqual, :ge => LA.BNNSRelationalOperatorGreaterEqual,
    :(!=) => LA.BNNSRelationalOperatorNotEqual, :ne => LA.BNNSRelationalOperatorNotEqual,
)

"""
    bnns_compare(op::Symbol, A::Array{Float32}, B::Array{Float32}) -> Array{Bool}

Element-wise relational comparison `A op B` via `BNNSCompareTensor`, returning a
`Bool` array. `op` is one of `:eq (==)`, `:lt (<)`, `:le (<=)`, `:gt (>)`,
`:ge (>=)`, `:ne (!=)`.
"""
function bnns_compare(op::Symbol, A::Array{Float32}, B::Array{Float32})
    haskey(_BNNS_RELOP, op) ||
        throw(ArgumentError("BNNS: unsupported comparison $(repr(op))"))
    size(A) == size(B) || throw(DimensionMismatch("A is $(size(A)), B is $(size(B))"))
    O = Array{Bool}(undef, size(A))
    da = _desc(A); db = _desc(B); do_ = _desc(O)
    GC.@preserve A B O begin
        _bnns_check(LA.BNNSCompareTensor(Ref(da), Ref(db),
                    LA.BNNSRelationalOperator(_BNNS_RELOP[op]), Ref(do_)),
                    "BNNSCompareTensor")
    end
    return O
end

"""
    bnns_band_part(A::Matrix, num_lower, num_upper) -> Matrix

Keep only the central band of `A` via `BNNSBandPart`: entries more than
`num_lower` below or `num_upper` above the diagonal (in Julia's column-major
sense) are zeroed. A negative count keeps the entire lower/upper triangle, so
`bnns_band_part(A, -1, 0)` is `tril(A)` and `bnns_band_part(A, 0, -1)` is
`triu(A)`.

!!! note
    BNNS counts bands in row-major orientation, so the underlying
    `num_lower`/`num_upper` are swapped internally to give the Julia-native
    `tril`/`triu` meaning documented here.
"""
function bnns_band_part(A::Array{T}, num_lower::Integer, num_upper::Integer) where {T}
    O = similar(A)
    di = _desc(A); do_ = _desc(O)
    # BNNS band orientation is row-major (transposed vs Julia col-major): pass the
    # counts swapped so callers get the intuitive tril/triu semantics.
    GC.@preserve A O begin
        _bnns_check(LA.BNNSBandPart(Cint(num_upper), Cint(num_lower),
                    Ref(di), Ref(do_), C_NULL), "BNNSBandPart")
    end
    return O
end

"""
    bnns_gather(dim::Integer, input::Array, indices::Array{Int32}) -> Array

Gather slices of `input` along Julia dimension `dim` using 0-based `indices`
(as BNNS/C expects) via `BNNSGather`, i.e. `selectdim(input, dim, indices.+1)`
stacked along `dim`. `indices` is a 1-D `Int32` vector.
"""
function bnns_gather(dim::Integer, input::Array{T}, indices::Array{Int32}) where {T}
    osz = collect(size(input)); osz[dim] = length(indices)
    O = Array{T}(undef, osz...)
    di = _desc(input); dind = _desc(indices); do_ = _desc(O)
    GC.@preserve input indices O begin
        _bnns_check(LA.BNNSGather(Csize_t(dim - 1), Ref(di), Ref(dind), Ref(do_), C_NULL),
                    "BNNSGather")
    end
    return O
end

const _BNNS_REDUCE = Dict(
    :max => LA.BNNSReduceFunctionMax, :min => LA.BNNSReduceFunctionMin,
    :sum => LA.BNNSReduceFunctionSum, :mean => LA.BNNSReduceFunctionMean,
    :sumsquare => LA.BNNSReduceFunctionSumSquare, :l1 => LA.BNNSReduceFunctionL1Norm,
    :l2 => LA.BNNSReduceFunctionL2Norm, :product => LA.BNNSReduceFunctionProduct,
    :logsumexp => LA.BNNSReduceFunctionLogSumExp,
)

"""
    bnns_scatter(dim::Integer, op::Symbol, input::Array, indices::Array{Int32}, output::Array) -> output

Scatter `input` into a copy of `output` along Julia dimension `dim` at 0-based
`indices`, combining collisions with reduce `op` (`:sum`, `:max`, `:min`, ...)
via `BNNSScatter`. `output` is updated in place and returned.
"""
function bnns_scatter(dim::Integer, op::Symbol, input::Array{T},
                      indices::Array{Int32}, output::Array{T}) where {T}
    haskey(_BNNS_REDUCE, op) || throw(ArgumentError("BNNS: unsupported reduce $(repr(op))"))
    di = _desc(input); dind = _desc(indices); do_ = _desc(output)
    GC.@preserve input indices output begin
        _bnns_check(LA.BNNSScatter(Csize_t(dim - 1), LA.BNNSReduceFunction(_BNNS_REDUCE[op]),
                    Ref(di), Ref(dind), Ref(do_), C_NULL), "BNNSScatter")
    end
    return output
end

"""
    bnns_gather_nd(input::Array, indices::Array{Int32}, output::Array) -> output

N-dimensional gather (`BNNSGatherND`): `indices` holds coordinate vectors along
its first (fastest) axis and `output` collects `input[coord...]`. Thin wrapper —
supply pre-shaped descriptors' backing arrays; `output` is filled and returned.
"""
function bnns_gather_nd(input::Array{T}, indices::Array{Int32}, output::Array{T}) where {T}
    di = _desc(input); dind = _desc(indices); do_ = _desc(output)
    GC.@preserve input indices output begin
        _bnns_check(LA.BNNSGatherND(Ref(di), Ref(dind), Ref(do_), C_NULL), "BNNSGatherND")
    end
    return output
end

"""
    bnns_scatter_nd(op::Symbol, input::Array, indices::Array{Int32}, output::Array) -> output

N-dimensional scatter (`BNNSScatterND`), the adjoint of [`bnns_gather_nd`](@ref),
combining collisions with reduce `op`. `output` is updated in place and returned.
"""
function bnns_scatter_nd(op::Symbol, input::Array{T}, indices::Array{Int32},
                         output::Array{T}) where {T}
    haskey(_BNNS_REDUCE, op) || throw(ArgumentError("BNNS: unsupported reduce $(repr(op))"))
    di = _desc(input); dind = _desc(indices); do_ = _desc(output)
    GC.@preserve input indices output begin
        _bnns_check(LA.BNNSScatterND(LA.BNNSReduceFunction(_BNNS_REDUCE[op]),
                    Ref(di), Ref(dind), Ref(do_), C_NULL), "BNNSScatterND")
    end
    return output
end

const _BNNS_SHUFFLE = Dict(
    :pixel_shuffle => LA.BNNSShuffleTypePixelShuffleNCHW,
    :pixel_unshuffle => LA.BNNSShuffleTypePixelUnshuffleNCHW,
    :depth_to_space => LA.BNNSShuffleTypeDepthToSpaceNCHW,
    :space_to_depth => LA.BNNSShuffleTypeSpaceToDepthNCHW,
)

"""
    bnns_shuffle(type::Symbol, input::Array, output::Array) -> output

Pixel / depth shuffle (`BNNSShuffle`). `type` is one of `:pixel_shuffle`,
`:pixel_unshuffle`, `:depth_to_space`, `:space_to_depth` (all NCHW). `output`
must be pre-sized for the transform and is returned.
"""
function bnns_shuffle(type::Symbol, input::Array{T}, output::Array{T}) where {T}
    haskey(_BNNS_SHUFFLE, type) || throw(ArgumentError("BNNS: unsupported shuffle $(repr(type))"))
    di = _desc(input); do_ = _desc(output)
    GC.@preserve input output begin
        _bnns_check(LA.BNNSShuffle(LA.BNNSShuffleType(_BNNS_SHUFFLE[type]),
                    Ref(di), Ref(do_), C_NULL), "BNNSShuffle")
    end
    return output
end

"""
    bnns_copy!(dest::Array, src::Array) -> dest

Copy (with BNNS's broadcasting / layout conversion rules) `src` into `dest` via
`BNNSCopy`. For equal shapes this is a plain element copy.
"""
function bnns_copy!(dest::Array{T}, src::Array{T}) where {T}
    ds = _desc(src); dd = _desc(dest)
    GC.@preserve src dest begin
        _bnns_check(LA.BNNSCopy(Ref(dd), Ref(ds), C_NULL), "BNNSCopy")
    end
    return dest
end

# =============================================================================
# Clipping / norm (stateless kernels)
# =============================================================================

"""
    bnns_clip_by_value(src::Array{Float32}, lo, hi) -> Array{Float32}

Element-wise clamp to `[lo, hi]` via `BNNSClipByValue` (`≡ clamp.(src, lo, hi)`).
"""
function bnns_clip_by_value(src::Array{Float32}, lo::Real, hi::Real)
    dest = similar(src)
    ds = _desc(src); dd = _desc(dest)
    GC.@preserve src dest begin
        _bnns_check(LA.BNNSClipByValue(Ref(dd), Ref(ds), Float32(lo), Float32(hi)),
                    "BNNSClipByValue")
    end
    return dest
end

"""
    bnns_clip_by_norm(src::Array{Float32}, max_norm; axis_flags=0) -> Array{Float32}

Scale `src` down so its L2 norm does not exceed `max_norm` via `BNNSClipByNorm`
(`axis_flags` selects the reduction axes as a bitmask; `0` = whole tensor).
"""
function bnns_clip_by_norm(src::Array{Float32}, max_norm::Real; axis_flags::Integer = 0)
    dest = similar(src)
    ds = _desc(src); dd = _desc(dest)
    GC.@preserve src dest begin
        _bnns_check(LA.BNNSClipByNorm(Ref(dd), Ref(ds), Float32(max_norm), UInt32(axis_flags)),
                    "BNNSClipByNorm")
    end
    return dest
end

"""
    bnns_clip_by_global_norm(srcs::Vector{<:Array{Float32}}, max_norm; use_norm=0) -> Vector{Array{Float32}}

Rescale a collection of tensors by a single global-norm factor via
`BNNSClipByGlobalNorm`. When `use_norm <= 0` the global norm is computed from
`srcs`; otherwise the supplied `use_norm` is used. Returns freshly-allocated
clipped copies.
"""
function bnns_clip_by_global_norm(srcs::Vector{<:Array{Float32}}, max_norm::Real;
                                  use_norm::Real = 0)
    n = length(srcs)
    dests = [similar(s) for s in srcs]
    sdescs = [Ref(_desc(s)) for s in srcs]
    ddescs = [Ref(_desc(d)) for d in dests]
    GC.@preserve srcs dests sdescs ddescs begin
        sptrs = [Base.unsafe_convert(Ptr{BNNSNDArrayDescriptor}, r) for r in sdescs]
        dptrs = [Base.unsafe_convert(Ptr{BNNSNDArrayDescriptor}, r) for r in ddescs]
        GC.@preserve sptrs dptrs begin
            _bnns_check(LA.BNNSClipByGlobalNorm(pointer(dptrs), pointer(sptrs),
                        Csize_t(n), Float32(max_norm), Float32(use_norm)),
                        "BNNSClipByGlobalNorm")
        end
    end
    return dests
end

"""
    bnns_compute_norm(src::Array{Float32}; axis_flags=0) -> Array{Float32}

L2 norm reduction over the axes selected by the `axis_flags` bitmask via
`BNNSComputeNorm` (`axis_flags = 0` reduces the whole tensor to a scalar). The
reduced axes collapse to length 1 in the returned array.
"""
function bnns_compute_norm(src::Array{Float32}; axis_flags::Integer = 0)
    if axis_flags == 0
        dest = zeros(Float32, ntuple(_ -> 1, ndims(src)))
    else
        osz = ntuple(i -> (axis_flags & (1 << (i - 1))) != 0 ? 1 : size(src, i), ndims(src))
        dest = zeros(Float32, osz)
    end
    ds = _desc(src); dd = _desc(dest)
    GC.@preserve src dest begin
        _bnns_check(LA.BNNSComputeNorm(Ref(dd), Ref(ds),
                    LA.BNNSNormType(LA.BNNSL2Norm), UInt32(axis_flags)), "BNNSComputeNorm")
    end
    return dest
end

# =============================================================================
# Utility / introspection queries (stateless)
# =============================================================================

"""
    bnns_layout_rank(layout::BNNSDataLayout) -> Int

Rank (number of dimensions) encoded by a `BNNSDataLayout` constant, via
`BNNSDataLayoutGetRank`.
"""
bnns_layout_rank(layout) = Int(LA.BNNSDataLayoutGetRank(BNNSDataLayout(layout)))

"""
    bnns_data_size(A::Array) -> Int

Number of bytes of tensor data described by `A` (`BNNSNDArrayGetDataSize`).
"""
function bnns_data_size(A::Array)
    d = _desc(A)
    GC.@preserve A begin
        return Int(LA.BNNSNDArrayGetDataSize(Ref(d)))
    end
end

"""
    bnns_tensor_allocation_size(A::Array) -> Int

Bytes required to allocate a `BNNSTensor` describing `A`
(`BNNSTensorGetAllocationSize`). Uses the modern `BNNSTensor` struct (rank +
shape/stride), distinct from the legacy `BNNSNDArrayDescriptor`.
"""
function bnns_tensor_allocation_size(A::Array{T}) where {T}
    n = ndims(A)
    shape = ntuple(i -> i <= n ? Cssize_t(size(A, i)) : Cssize_t(0), 8)
    st = ntuple(i -> i <= n ? Cssize_t(stride(A, i)) : Cssize_t(0), 8)
    t = LA.BNNSTensor(_bnns_dt(T), UInt8(n), shape, st, Ptr{Cvoid}(pointer(A)),
                      Csize_t(sizeof(A)), C_NULL)
    GC.@preserve A begin
        return Int(LA.BNNSTensorGetAllocationSize(Ref(t)))
    end
end

# =============================================================================
# DirectApply family (stateless kernels operating on descriptors)
# =============================================================================

"""
    bnns_broadcast_matmul(A, B; transA=false, transB=false, alpha=1f0) -> Array{Float32}

Batched / broadcasting matrix multiply `alpha * A * B` via
`BNNSDirectApplyBroadcastMatMul`. For 2-D inputs this equals
`alpha .* (A * B)`; higher-rank inputs multiply the trailing two axes and
broadcast the batch axes.
"""
function bnns_broadcast_matmul(A::Array{Float32}, B::Array{Float32};
                               transA::Bool = false, transB::Bool = false,
                               alpha::Real = 1.0f0)
    m = transA ? size(A, 2) : size(A, 1)
    n = transB ? size(B, 1) : size(B, 2)
    batch = ntuple(i -> max(size(A, i + 2), size(B, i + 2)), max(ndims(A), ndims(B)) - 2)
    O = Array{Float32}(undef, m, n, batch...)
    da = _desc(A); db = _desc(B); do_ = _desc(O)
    GC.@preserve A B O begin
        LA.BNNSDirectApplyBroadcastMatMul(transA, transB, Float32(alpha),
            Ref(da), Ref(db), Ref(do_), C_NULL)
    end
    return O
end

"""
    bnns_topk(input::Array{Float32}, K; dim=1) -> (values, indices)

Top-`K` values and their 0-based indices along Julia dimension `dim` via
`BNNSDirectApplyTopK`. `values` is `Float32`, `indices` is `Int32`. Comparable to
`sort`-based `partialsortperm` per slice.
"""
function bnns_topk(input::Array{Float32}, K::Integer; dim::Integer = 1)
    osz = collect(size(input)); osz[dim] = K
    vals = Array{Float32}(undef, osz...)
    inds = Array{Int32}(undef, osz...)
    di = _desc(input); dv = _desc(vals); dind = _desc(inds)
    GC.@preserve input vals inds begin
        _bnns_check(LA.BNNSDirectApplyTopK(Csize_t(K), Csize_t(dim - 1), Csize_t(1),
            Ref(di), Csize_t(0), Ref(dv), Csize_t(0), Ref(dind), Csize_t(0), C_NULL),
            "BNNSDirectApplyTopK")
    end
    return vals, inds
end

"""
    bnns_reduce(func::Symbol, input::Array{Float32}; dim=1) -> Array{Float32}

Reduce `input` along Julia dimension `dim` with `func`
(`:sum`, `:mean`, `:max`, `:min`, `:sumsquare`, `:l1`, `:l2`, `:product`,
`:logsumexp`) via `BNNSDirectApplyReduction`. The reduced axis collapses to
length 1.
"""
function bnns_reduce(func::Symbol, input::Array{Float32}; dim::Integer = 1)
    haskey(_BNNS_REDUCE, func) || throw(ArgumentError("BNNS: unsupported reduce $(repr(func))"))
    osz = collect(size(input)); osz[dim] = 1
    O = zeros(Float32, osz...)
    di = _desc(input); do_ = _desc(O)
    lp = Ref(LA.BNNSLayerParametersReduction(di, do_, _empty_desc(),
            LA.BNNSReduceFunction(_BNNS_REDUCE[func]), 0.0f0))
    GC.@preserve input O lp begin
        _bnns_check(LA.BNNSDirectApplyReduction(
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersReduction}, lp), C_NULL),
            "BNNSDirectApplyReduction")
    end
    return O
end

# An all-zero descriptor for optional slots (no backing data).
_empty_desc() = BNNSNDArrayDescriptor(LA.BNNSNDArrayFlags(0), _bnns_layout(1),
    ntuple(_ -> Csize_t(0), 8), ntuple(_ -> Csize_t(0), 8), C_NULL,
    _bnns_dt(Float32), C_NULL, _bnns_dt(Float32), 1.0f0, 0.0f0)

# =============================================================================
# Random number generation (opaque-handle resource)
# =============================================================================

"""
    BNNSRandomGenerator([seed]) -> BNNSRandomGenerator

A BNNS random number generator handle (AES-CTR method). Construct with an
optional 64-bit `seed` for reproducibility (`BNNSCreateRandomGeneratorWithSeed`,
or `BNNSCreateRandomGenerator` when omitted). The handle is destroyed
automatically by a finalizer (`BNNSDestroyRandomGenerator`).

Use with [`bnns_random_fill_uniform!`](@ref),
[`bnns_random_fill_normal!`](@ref), [`bnns_random_fill_uniform_int!`](@ref),
[`bnns_random_fill_categorical!`](@ref) and the
[`bnns_random_state`](@ref)/[`bnns_random_state!`](@ref) round-trip.
"""
mutable struct BNNSRandomGenerator
    handle::Ptr{Cvoid}
    function BNNSRandomGenerator(seed::Union{Integer,Nothing} = nothing)
        h = seed === nothing ?
            LA.BNNSCreateRandomGenerator(LA.BNNSRandomGeneratorMethod(LA.BNNSRandomGeneratorMethodAES_CTR), C_NULL) :
            LA.BNNSCreateRandomGeneratorWithSeed(LA.BNNSRandomGeneratorMethod(LA.BNNSRandomGeneratorMethodAES_CTR), UInt64(seed), C_NULL)
        h == C_NULL && error("BNNSCreateRandomGenerator returned NULL")
        g = new(h)
        finalizer(g) do x
            x.handle == C_NULL || LA.BNNSDestroyRandomGenerator(x.handle)
            x.handle = C_NULL
        end
        return g
    end
end

"""
    bnns_random_fill_uniform!(g::BNNSRandomGenerator, A::Array{Float32}, lo=0f0, hi=1f0) -> A

Fill `A` with i.i.d. uniform samples on `[lo, hi)` (`BNNSRandomFillUniformFloat`).
"""
function bnns_random_fill_uniform!(g::BNNSRandomGenerator, A::Array{Float32},
                                   lo::Real = 0.0f0, hi::Real = 1.0f0)
    d = _desc(A)
    GC.@preserve A begin
        _bnns_check(LA.BNNSRandomFillUniformFloat(g.handle, Ref(d), Float32(lo), Float32(hi)),
                    "BNNSRandomFillUniformFloat")
    end
    return A
end

"""
    bnns_random_fill_uniform_int!(g::BNNSRandomGenerator, A::Array{Int32}, lo, hi) -> A

Fill integer array `A` with i.i.d. uniform samples on the half-open range
`[lo, hi)` (`BNNSRandomFillUniformInt`).
"""
function bnns_random_fill_uniform_int!(g::BNNSRandomGenerator, A::Array{Int32},
                                       lo::Integer, hi::Integer)
    d = _desc(A)
    GC.@preserve A begin
        _bnns_check(LA.BNNSRandomFillUniformInt(g.handle, Ref(d), Int64(lo), Int64(hi)),
                    "BNNSRandomFillUniformInt")
    end
    return A
end

"""
    bnns_random_fill_normal!(g::BNNSRandomGenerator, A::Array{Float32}, mean=0f0, stddev=1f0) -> A

Fill `A` with i.i.d. Gaussian samples (`BNNSRandomFillNormalFloat`).
"""
function bnns_random_fill_normal!(g::BNNSRandomGenerator, A::Array{Float32},
                                  mean::Real = 0.0f0, stddev::Real = 1.0f0)
    d = _desc(A)
    GC.@preserve A begin
        _bnns_check(LA.BNNSRandomFillNormalFloat(g.handle, Ref(d), Float32(mean), Float32(stddev)),
                    "BNNSRandomFillNormalFloat")
    end
    return A
end

"""
    bnns_random_fill_categorical!(g::BNNSRandomGenerator, out::Array{Float32}, probs::Vector{Float32}; log_probs=false) -> out

Draw categorical samples (0-based category indices, stored as `Float32`) into
`out` using per-category weights `probs` (`BNNSRandomFillCategoricalFloat`). Pass
`log_probs=true` if `probs` holds log probabilities.
"""
function bnns_random_fill_categorical!(g::BNNSRandomGenerator, out::Array{Float32},
                                       probs::Array{Float32}; log_probs::Bool = false)
    dout = _desc(out); dp = _desc(probs)
    GC.@preserve out probs begin
        _bnns_check(LA.BNNSRandomFillCategoricalFloat(g.handle, Ref(dout), Ref(dp), log_probs),
                    "BNNSRandomFillCategoricalFloat")
    end
    return out
end

"""
    bnns_random_state(g::BNNSRandomGenerator) -> Vector{UInt8}

Snapshot the generator's internal state (`BNNSRandomGeneratorStateSize` +
`BNNSRandomGeneratorGetState`). Restore it with [`bnns_random_state!`](@ref).
"""
function bnns_random_state(g::BNNSRandomGenerator)
    sz = Int(LA.BNNSRandomGeneratorStateSize(g.handle))
    state = Vector{UInt8}(undef, sz)
    GC.@preserve state begin
        _bnns_check(LA.BNNSRandomGeneratorGetState(g.handle, Csize_t(sz), pointer(state)),
                    "BNNSRandomGeneratorGetState")
    end
    return state
end

"""
    bnns_random_state!(g::BNNSRandomGenerator, state::Vector{UInt8}) -> g

Restore a generator state captured by [`bnns_random_state`](@ref)
(`BNNSRandomGeneratorSetState`).
"""
function bnns_random_state!(g::BNNSRandomGenerator, state::Vector{UInt8})
    GC.@preserve state begin
        _bnns_check(LA.BNNSRandomGeneratorSetState(g.handle, Csize_t(length(state)), pointer(state)),
                    "BNNSRandomGeneratorSetState")
    end
    return g
end

# =============================================================================
# k-Nearest-Neighbours (opaque-handle resource)
# =============================================================================

"""
    BNNSNearestNeighbors(max_samples, n_features, n_neighbors; T=Float32) -> BNNSNearestNeighbors

A brute-force k-nearest-neighbours index (`BNNSCreateNearestNeighbors`) holding up
to `max_samples` reference points of dimension `n_features`, answering
`n_neighbors`-NN queries. Destroyed automatically (`BNNSDestroyNearestNeighbors`).

Add reference points with [`bnns_knn_load!`](@ref) and query with
[`bnns_knn_query`](@ref).
"""
mutable struct BNNSNearestNeighbors
    handle::Ptr{Cvoid}
    n_neighbors::Int
    function BNNSNearestNeighbors(max_samples::Integer, n_features::Integer,
                                  n_neighbors::Integer; T::Type = Float32)
        h = LA.BNNSCreateNearestNeighbors(Cuint(max_samples), Cuint(n_features),
                Cuint(n_neighbors), _bnns_dt(T), C_NULL)
        h == C_NULL && error("BNNSCreateNearestNeighbors returned NULL")
        knn = new(h, Int(n_neighbors))
        finalizer(knn) do x
            x.handle == C_NULL || LA.BNNSDestroyNearestNeighbors(x.handle)
            x.handle = C_NULL
        end
        return knn
    end
end

"""
    bnns_knn_load!(knn::BNNSNearestNeighbors, data::Matrix{Float32}) -> Int

Append reference samples to the index (`BNNSNearestNeighborsLoad`). `data` is
`n_features × n_new_samples` (each column is one sample, matching BNNS's
feature-major layout). Returns the number of samples loaded.
"""
function bnns_knn_load!(knn::BNNSNearestNeighbors, data::Matrix{Float32})
    n_new = size(data, 2)
    GC.@preserve data begin
        _bnns_check(LA.BNNSNearestNeighborsLoad(knn.handle, Cuint(n_new), pointer(data)),
                    "BNNSNearestNeighborsLoad")
    end
    return n_new
end

"""
    bnns_knn_query(knn::BNNSNearestNeighbors, sample_number) -> (indices, distances)

Return the `n_neighbors` nearest reference points to the (0-based) loaded sample
`sample_number` (`BNNSNearestNeighborsGetInfo`): their 0-based `indices`
(`Vector{Int32}`) and `Float32` `distances`.
"""
function bnns_knn_query(knn::BNNSNearestNeighbors, sample_number::Integer)
    k = knn.n_neighbors
    indices = Vector{Int32}(undef, k)
    distances = Vector{Float32}(undef, k)
    GC.@preserve indices distances begin
        _bnns_check(LA.BNNSNearestNeighborsGetInfo(knn.handle, Cint(sample_number),
                    pointer(indices), pointer(distances)), "BNNSNearestNeighborsGetInfo")
    end
    return indices, distances
end

# =============================================================================
# BNNS Graph API (the modern, non-deprecated compiler pipeline).
#
# Lifecycle: build compile options -> `BNNSGraphCompileFromFile` a serialized
# `.bnns`/MIL package into a `BNNSGraph` -> `BNNSGraphContextMake` an executable
# context -> introspect / `execute`. Compiling requires an on-disk graph package
# (there is no in-memory graph builder in this API), so the wrappers here cover
# the full surface but the graph/context calls need a real package to exercise.
# =============================================================================

# The raw BNNSGraph symbols take `Ptr{Cchar}` (not `Cstring`); build a
# null-terminated byte buffer and pass `pointer` under `GC.@preserve`.
_cstr(s) = push!(Vector{UInt8}(codeunits(String(s))), 0x00)
_cptr(buf::Vector{UInt8}) = Ptr{Cchar}(pointer(buf))

const _OPT_PREF = Dict(:performance => LA.BNNSGraphOptimizationPreferencePerformance,
                       :ir_size => LA.BNNSGraphOptimizationPreferenceIRSize)
const _OPT_PREF_REV = Dict(UInt32(v) => k for (k, v) in _OPT_PREF)

"""
    BNNSGraphCompileOptions(; single_thread=nothing, generate_debug_info=nothing,
                              optimization=nothing, log_mask=nothing,
                              output_path=nothing, output_fd=nothing) -> BNNSGraphCompileOptions

Options controlling `BNNSGraphCompileFromFile`, backed by
`BNNSGraphCompileOptionsMakeDefault` and destroyed by a finalizer
(`BNNSGraphCompileOptionsDestroy`). Any keyword left `nothing` keeps the BNNS
default. `optimization` is `:performance` or `:ir_size`. Individual fields can
also be read/written with the accessor functions below.
"""
mutable struct BNNSGraphCompileOptions
    opts::LA.bnns_graph_compile_options_t
    function BNNSGraphCompileOptions(; single_thread = nothing, generate_debug_info = nothing,
                                     optimization = nothing, log_mask = nothing,
                                     output_path = nothing, output_fd = nothing)
        o = new(LA.BNNSGraphCompileOptionsMakeDefault())
        finalizer(x -> LA.BNNSGraphCompileOptionsDestroy(x.opts), o)
        single_thread === nothing || bnns_compile_options_set_single_thread!(o, single_thread)
        generate_debug_info === nothing || bnns_compile_options_set_debug_info!(o, generate_debug_info)
        optimization === nothing || bnns_compile_options_set_optimization!(o, optimization)
        log_mask === nothing || bnns_compile_options_set_log_mask!(o, log_mask)
        output_path === nothing || bnns_compile_options_set_output_path!(o, output_path)
        output_fd === nothing || bnns_compile_options_set_output_fd!(o, output_fd)
        return o
    end
end

"Set the single-threaded-target flag (`BNNSGraphCompileOptionsSetTargetSingleThread`)."
bnns_compile_options_set_single_thread!(o::BNNSGraphCompileOptions, v::Bool) =
    (LA.BNNSGraphCompileOptionsSetTargetSingleThread(o.opts, v); o)
"Query the single-threaded-target flag (`BNNSGraphCompileOptionsGetTargetSingleThread`)."
bnns_compile_options_get_single_thread(o::BNNSGraphCompileOptions) =
    LA.BNNSGraphCompileOptionsGetTargetSingleThread(o.opts)

"Set the debug-info flag (`BNNSGraphCompileOptionsSetGenerateDebugInfo`)."
bnns_compile_options_set_debug_info!(o::BNNSGraphCompileOptions, v::Bool) =
    (LA.BNNSGraphCompileOptionsSetGenerateDebugInfo(o.opts, v); o)
"Query the debug-info flag (`BNNSGraphCompileOptionsGetGenerateDebugInfo`)."
bnns_compile_options_get_debug_info(o::BNNSGraphCompileOptions) =
    LA.BNNSGraphCompileOptionsGetGenerateDebugInfo(o.opts)

"Set the optimization preference `:performance`/`:ir_size` (`...SetOptimizationPreference`)."
function bnns_compile_options_set_optimization!(o::BNNSGraphCompileOptions, pref::Symbol)
    haskey(_OPT_PREF, pref) || throw(ArgumentError("BNNS: optimization must be :performance or :ir_size"))
    LA.BNNSGraphCompileOptionsSetOptimizationPreference(o.opts, LA.BNNSGraphOptimizationPreference(_OPT_PREF[pref]))
    return o
end
"Query the optimization preference (`...GetOptimizationPreference`)."
bnns_compile_options_get_optimization(o::BNNSGraphCompileOptions) =
    _OPT_PREF_REV[UInt32(LA.BNNSGraphCompileOptionsGetOptimizationPreference(o.opts))]

"Set the message-log level bitmask (`BNNSGraphCompileOptionsSetMessageLogMask`)."
bnns_compile_options_set_log_mask!(o::BNNSGraphCompileOptions, mask::Integer) =
    (LA.BNNSGraphCompileOptionsSetMessageLogMask(o.opts, UInt32(mask)); o)

"Install a C message-log callback (`BNNSGraphCompileOptionsSetMessageLogCallback`)."
bnns_compile_options_set_log_callback!(o::BNNSGraphCompileOptions, cb::Ptr{Cvoid},
                                       data::Ptr = C_NULL) =
    (LA.BNNSGraphCompileOptionsSetMessageLogCallback(o.opts, cb, Ptr{LA.bnns_user_message_data_t}(data)); o)

"Set the compiled-artifact output path (`BNNSGraphCompileOptionsSetOutputPath`)."
function bnns_compile_options_set_output_path!(o::BNNSGraphCompileOptions, path::AbstractString)
    buf = _cstr(path)
    GC.@preserve buf LA.BNNSGraphCompileOptionsSetOutputPath(o.opts, _cptr(buf))
    return o
end
"Query the compiled-artifact output path (`BNNSGraphCompileOptionsGetOutputPath`)."
function bnns_compile_options_get_output_path(o::BNNSGraphCompileOptions)
    p = LA.BNNSGraphCompileOptionsGetOutputPath(o.opts)
    p == C_NULL ? "" : unsafe_string(p)
end

"Set the compiled-artifact output file descriptor (`BNNSGraphCompileOptionsSetOutputFD`)."
bnns_compile_options_set_output_fd!(o::BNNSGraphCompileOptions, fd::Integer) =
    (LA.BNNSGraphCompileOptionsSetOutputFD(o.opts, Cint(fd)); o)
"Query the compiled-artifact output file descriptor (`BNNSGraphCompileOptionsGetOutputFD`)."
bnns_compile_options_get_output_fd(o::BNNSGraphCompileOptions) =
    Int(LA.BNNSGraphCompileOptionsGetOutputFD(o.opts))

# --- Graph object ------------------------------------------------------------

"""
    BNNSGraph(filename; func=nothing, options=BNNSGraphCompileOptions()) -> BNNSGraph

Compile a serialized BNNS graph package at `filename` (optionally selecting a
named `func` inside it) into an executable graph via `BNNSGraphCompileFromFile`.
The returned handle feeds [`BNNSGraphContext`](@ref) and the graph-introspection
helpers.
"""
mutable struct BNNSGraph
    graph::LA.bnns_graph_t
    function BNNSGraph(filename::AbstractString; func = nothing,
                       options::BNNSGraphCompileOptions = BNNSGraphCompileOptions())
        fn = _cstr(filename)
        fckeep, fcptr = _fnarg(func)
        g = GC.@preserve fn fckeep options begin
            LA.BNNSGraphCompileFromFile(_cptr(fn), fcptr, options.opts)
        end
        return new(g)
    end
end

# Return `(keep, ptr)` where `keep` roots the buffer that `ptr` points into. For
# a `nothing` function name pass a NULL `Ptr{Cchar}`.
function _fnarg(func)
    func === nothing && return (nothing, Ptr{Cchar}(C_NULL))
    buf = _cstr(func)
    return (buf, _cptr(buf))
end

"Number of graph functions (`BNNSGraphGetFunctionCount`)."
bnns_graph_function_count(g::BNNSGraph) = Int(LA.BNNSGraphGetFunctionCount(g.graph))

"Number of inputs to `func` (`BNNSGraphGetInputCount`)."
function bnns_graph_input_count(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func)
    GC.@preserve keep Int(LA.BNNSGraphGetInputCount(g.graph, fptr))
end
"Number of outputs of `func` (`BNNSGraphGetOutputCount`)."
function bnns_graph_output_count(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func)
    GC.@preserve keep Int(LA.BNNSGraphGetOutputCount(g.graph, fptr))
end
"Number of arguments of `func` (`BNNSGraphGetArgumentCount`)."
function bnns_graph_argument_count(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func)
    GC.@preserve keep Int(LA.BNNSGraphGetArgumentCount(g.graph, fptr))
end

_names_from(cptrs) = String[p == C_NULL ? "" : unsafe_string(p) for p in cptrs]

"Names of the graph's functions (`BNNSGraphGetFunctionNames`)."
function bnns_graph_function_names(g::BNNSGraph)
    n = bnns_graph_function_count(g)
    buf = fill(Ptr{Cchar}(C_NULL), n)
    GC.@preserve buf _bnns_check(LA.BNNSGraphGetFunctionNames(g.graph, Csize_t(n), pointer(buf)), "BNNSGraphGetFunctionNames")
    return _names_from(buf)
end
"Input names of `func` (`BNNSGraphGetInputNames`)."
function bnns_graph_input_names(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func); n = bnns_graph_input_count(g, func)
    buf = fill(Ptr{Cchar}(C_NULL), n)
    GC.@preserve keep buf _bnns_check(LA.BNNSGraphGetInputNames(g.graph, fptr, Csize_t(n), pointer(buf)), "BNNSGraphGetInputNames")
    return _names_from(buf)
end
"Output names of `func` (`BNNSGraphGetOutputNames`)."
function bnns_graph_output_names(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func); n = bnns_graph_output_count(g, func)
    buf = fill(Ptr{Cchar}(C_NULL), n)
    GC.@preserve keep buf _bnns_check(LA.BNNSGraphGetOutputNames(g.graph, fptr, Csize_t(n), pointer(buf)), "BNNSGraphGetOutputNames")
    return _names_from(buf)
end
"Argument names of `func` (`BNNSGraphGetArgumentNames`)."
function bnns_graph_argument_names(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func); n = bnns_graph_argument_count(g, func)
    buf = fill(Ptr{Cchar}(C_NULL), n)
    GC.@preserve keep buf _bnns_check(LA.BNNSGraphGetArgumentNames(g.graph, fptr, Csize_t(n), pointer(buf)), "BNNSGraphGetArgumentNames")
    return _names_from(buf)
end

const _INTENT_REV = Dict(UInt32(LA.BNNSGraphArgumentIntentIn) => :in,
                         UInt32(LA.BNNSGraphArgumentIntentOut) => :out,
                         UInt32(LA.BNNSGraphArgumentIntentInOut) => :inout)

"Per-argument intents (`:in`/`:out`/`:inout`) of `func` (`BNNSGraphGetArgumentIntents`)."
function bnns_graph_argument_intents(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func); n = bnns_graph_argument_count(g, func)
    buf = zeros(LA.BNNSGraphArgumentIntent, n)
    GC.@preserve keep buf _bnns_check(LA.BNNSGraphGetArgumentIntents(g.graph, fptr, Csize_t(n), pointer(buf)), "BNNSGraphGetArgumentIntents")
    return [get(_INTENT_REV, UInt32(v), v) for v in buf]
end

"0-based position of a named `argument` within `func` (`BNNSGraphGetArgumentPosition`)."
function bnns_graph_argument_position(g::BNNSGraph, argument::AbstractString, func = nothing)
    keep, fptr = _fnarg(func); a = _cstr(argument)
    GC.@preserve keep a Int(LA.BNNSGraphGetArgumentPosition(g.graph, fptr, _cptr(a)))
end

"Per-argument interleave factors of `func` (`BNNSGraphGetArgumentInterleaveFactors`)."
function bnns_graph_argument_interleave_factors(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func); n = bnns_graph_argument_count(g, func)
    ptrbuf = fill(Ptr{UInt16}(C_NULL), n)
    counts = zeros(Csize_t, n)
    GC.@preserve keep ptrbuf counts _bnns_check(
        LA.BNNSGraphGetArgumentInterleaveFactors(g.graph, fptr, Csize_t(n), pointer(ptrbuf), pointer(counts)),
        "BNNSGraphGetArgumentInterleaveFactors")
    return [p == C_NULL ? UInt16[] : unsafe_wrap(Vector{UInt16}, p, Int(c); own = false) |> copy
            for (p, c) in zip(ptrbuf, counts)]
end

"""
    bnns_graph_fill_strides!(g::BNNSGraph, argument, tensor::Ref{BNNSTensor}; func=nothing) -> tensor

Populate the stride fields of `tensor` for the named `argument` of `func` from
the graph's known layout (`BNNSGraphTensorFillStrides`).
"""
function bnns_graph_fill_strides!(g::BNNSGraph, argument::AbstractString,
                                  tensor::Ref{LA.BNNSTensor}; func = nothing)
    keep, fptr = _fnarg(func); a = _cstr(argument)
    GC.@preserve keep a tensor _bnns_check(
        LA.BNNSGraphTensorFillStrides(g.graph, fptr, _cptr(a),
            Base.unsafe_convert(Ptr{LA.BNNSTensor}, tensor)), "BNNSGraphTensorFillStrides")
    return tensor
end

# --- Graph execution context -------------------------------------------------

"""
    BNNSGraphContext(g::BNNSGraph) -> BNNSGraphContext

An executable context for a compiled [`BNNSGraph`](@ref) (`BNNSGraphContextMake`),
destroyed by a finalizer (`BNNSGraphContextDestroy`). Feed it to
[`bnns_graph_execute!`](@ref).
"""
mutable struct BNNSGraphContext
    ctx::LA.bnns_graph_context_t
    function BNNSGraphContext(c::LA.bnns_graph_context_t)
        ctx = new(c)
        finalizer(x -> LA.BNNSGraphContextDestroy(x.ctx), ctx)
        return ctx
    end
end

BNNSGraphContext(g::BNNSGraph) = BNNSGraphContext(LA.BNNSGraphContextMake(g.graph))

"""
    BNNSGraphContext(g, func, initial_states::Vector{BNNSTensor}) -> BNNSGraphContext

A streaming context (`BNNSGraphContextMakeStreaming`) seeded with recurrent
`initial_states`. Same finalizer as the plain constructor.
"""
function BNNSGraphContext(g::BNNSGraph, func, initial_states::Vector{LA.BNNSTensor})
    keep, fptr = _fnarg(func)
    c = GC.@preserve keep initial_states LA.BNNSGraphContextMakeStreaming(
        g.graph, fptr, Csize_t(length(initial_states)), pointer(initial_states))
    return BNNSGraphContext(c)
end

const _ARGTYPE = Dict(:pointer => LA.BNNSGraphArgumentTypePointer, :tensor => LA.BNNSGraphArgumentTypeTensor)

"Set how execute arguments are interpreted, `:pointer` or `:tensor` (`BNNSGraphContextSetArgumentType`)."
function bnns_graph_context_set_argument_type!(c::BNNSGraphContext, t::Symbol)
    haskey(_ARGTYPE, t) || throw(ArgumentError("BNNS: argument type must be :pointer or :tensor"))
    _bnns_check(LA.BNNSGraphContextSetArgumentType(c.ctx, LA.BNNSGraphArgumentType(_ARGTYPE[t])), "BNNSGraphContextSetArgumentType")
    return c
end

"Enable/disable NaN & Inf checking during execute (`BNNSGraphContextEnableNanAndInfChecks`)."
bnns_graph_context_enable_nan_inf_checks!(c::BNNSGraphContext, enable::Bool = true) =
    (LA.BNNSGraphContextEnableNanAndInfChecks(c.ctx, enable); c)

"Set the streaming advance count (`BNNSGraphContextSetStreamingAdvanceCount`)."
bnns_graph_context_set_streaming_advance_count!(c::BNNSGraphContext, n::Integer) =
    (_bnns_check(LA.BNNSGraphContextSetStreamingAdvanceCount(c.ctx, Csize_t(n)), "BNNSGraphContextSetStreamingAdvanceCount"); c)

"Set the execute message-log level bitmask (`BNNSGraphContextSetMessageLogMask`)."
bnns_graph_context_set_log_mask!(c::BNNSGraphContext, mask::Integer) =
    (_bnns_check(LA.BNNSGraphContextSetMessageLogMask(c.ctx, UInt32(mask)), "BNNSGraphContextSetMessageLogMask"); c)

"Install a C execute message-log callback (`BNNSGraphContextSetMessageLogCallback`)."
bnns_graph_context_set_log_callback!(c::BNNSGraphContext, cb::Ptr{Cvoid}, data::Ptr = C_NULL) =
    (_bnns_check(LA.BNNSGraphContextSetMessageLogCallback(c.ctx, cb, Ptr{LA.bnns_user_message_data_t}(data)), "BNNSGraphContextSetMessageLogCallback"); c)

"Workspace size (bytes) required to execute `func` (`BNNSGraphContextGetWorkspaceSize`)."
function bnns_graph_context_workspace_size(c::BNNSGraphContext, func = nothing)
    keep, fptr = _fnarg(func)
    GC.@preserve keep Int(LA.BNNSGraphContextGetWorkspaceSize(c.ctx, fptr))
end

"""
    bnns_graph_context_set_dynamic_shapes!(c, shapes::Vector{bnns_graph_shape_t}; func=nothing) -> c

Bind concrete shapes to a graph compiled with dynamic dimensions
(`BNNSGraphContextSetDynamicShapes`).
"""
function bnns_graph_context_set_dynamic_shapes!(c::BNNSGraphContext,
                                                shapes::Vector{LA.bnns_graph_shape_t}; func = nothing)
    keep, fptr = _fnarg(func)
    GC.@preserve keep shapes _bnns_check(
        LA.BNNSGraphContextSetDynamicShapes(c.ctx, fptr, Csize_t(length(shapes)), pointer(shapes)),
        "BNNSGraphContextSetDynamicShapes")
    return c
end

"""
    bnns_graph_context_get_tensor(c, argument; func=nothing, fill_shapes=true) -> BNNSTensor

Fetch the `BNNSTensor` descriptor for a named `argument`
(`BNNSGraphContextGetTensor`).
"""
function bnns_graph_context_get_tensor(c::BNNSGraphContext, argument::AbstractString;
                                       func = nothing, fill_shapes::Bool = true)
    keep, fptr = _fnarg(func); a = _cstr(argument)
    t = Ref{LA.BNNSTensor}()
    GC.@preserve keep a t _bnns_check(
        LA.BNNSGraphContextGetTensor(c.ctx, fptr, _cptr(a), fill_shapes,
            Base.unsafe_convert(Ptr{LA.BNNSTensor}, t)), "BNNSGraphContextGetTensor")
    return t[]
end

"""
    bnns_graph_execute!(c, arguments::Vector{bnns_graph_argument_t}; func=nothing, workspace=UInt8[]) -> c

Execute `func` with the supplied argument buffers (`BNNSGraphContextExecute`).
Size `workspace` from [`bnns_graph_context_workspace_size`](@ref).
"""
function bnns_graph_execute!(c::BNNSGraphContext, arguments::Vector{LA.bnns_graph_argument_t};
                             func = nothing, workspace::Vector{UInt8} = UInt8[])
    keep, fptr = _fnarg(func)
    GC.@preserve keep arguments workspace begin
        wptr = isempty(workspace) ? Ptr{Cchar}(C_NULL) : Ptr{Cchar}(pointer(workspace))
        _bnns_check(LA.BNNSGraphContextExecute(c.ctx, fptr, Csize_t(length(arguments)),
            pointer(arguments), Csize_t(length(workspace)), wptr), "BNNSGraphContextExecute")
    end
    return c
end

# =============================================================================
# Layer filters (the classic BNNSFilter create/apply API).
#
# Apple DEPRECATED this whole family (macOS 15 / iOS 18) in favour of the BNNS
# Graph API above. We still expose it because it remains the only in-process way
# to build individual layers without a serialized graph package. Every filter is
# an opaque handle wrapped in `BNNSFilter` (create -> apply -> destroy, with a
# finalizer). The high-level [`bnns_fully_connected`](@ref) constructor is fully
# numerically verified; the remaining `bnns_create_layer_*` constructors are thin
# passthroughs over the raw parameter structs for callers that need them.
# =============================================================================

"""
    BNNSFilter

An opaque BNNS layer-filter handle (create → [`bnns_filter_apply!`](@ref) →
destroy). Destroyed automatically by a finalizer (`BNNSFilterDestroy`). Backing
weight/bias/parameter arrays passed at construction are kept alive for the
filter's lifetime.

!!! warning
    The BNNSFilter layer API is deprecated by Apple (macOS 15+) in favour of the
    BNNS Graph API. Prefer [`BNNSGraph`](@ref) for new code.
"""
mutable struct BNNSFilter
    handle::Ptr{Cvoid}
    roots::Vector{Any}
    function BNNSFilter(handle::Ptr{Cvoid}, roots::Vector{Any} = Any[])
        handle == C_NULL && error("BNNS filter creation returned NULL")
        f = new(handle, roots)
        finalizer(f) do x
            x.handle == C_NULL || LA.BNNSFilterDestroy(x.handle)
            x.handle = C_NULL
        end
        return f
    end
end

_default_fp() = LA.BNNSFilterParameters(UInt32(0), Csize_t(0), C_NULL, C_NULL)

const _BNNS_ACT_ALL = Dict(
    :identity => LA.BNNSActivationFunctionIdentity,
    :relu => LA.BNNSActivationFunctionRectifiedLinear,
    :sigmoid => LA.BNNSActivationFunctionSigmoid,
    :tanh => LA.BNNSActivationFunctionTanh,
    :abs => LA.BNNSActivationFunctionAbs,
    :gelu => LA.BNNSActivationFunctionGELU,
    :softmax => LA.BNNSActivationFunctionSoftmax,
)

function _activation_struct(f::Symbol; alpha = 0.0f0, beta = 0.0f0)
    haskey(_BNNS_ACT_ALL, f) || throw(ArgumentError("BNNS: unsupported activation $(repr(f))"))
    LA.BNNSActivation(BNNSActivationFunction(_BNNS_ACT_ALL[f]), Float32(alpha), Float32(beta),
                      Int32(0), Int32(0), Int32(0), C_NULL, C_NULL, C_NULL)
end

# Shape-only vector / column-major-matrix descriptors used by the FC layer.
_vec_shape_desc(sz; data = C_NULL) = BNNSNDArrayDescriptor(
    LA.BNNSNDArrayFlags(0), BNNSDataLayout(LA.BNNSDataLayoutVector),
    ntuple(i -> i == 1 ? Csize_t(sz) : Csize_t(0), 8),
    ntuple(i -> i == 1 ? Csize_t(1) : Csize_t(0), 8),
    Ptr{Cvoid}(data), _bnns_dt(Float32), C_NULL, _bnns_dt(Float32), 1.0f0, 0.0f0)

_col_matrix_desc(rows, cols; data = C_NULL) = BNNSNDArrayDescriptor(
    LA.BNNSNDArrayFlags(0), BNNSDataLayout(LA.BNNSDataLayoutColumnMajorMatrix),
    ntuple(i -> i == 1 ? Csize_t(rows) : i == 2 ? Csize_t(cols) : Csize_t(0), 8),
    ntuple(i -> i == 1 ? Csize_t(1) : i == 2 ? Csize_t(rows) : Csize_t(0), 8),
    Ptr{Cvoid}(data), _bnns_dt(Float32), C_NULL, _bnns_dt(Float32), 1.0f0, 0.0f0)

"""
    bnns_fully_connected(W::Matrix{Float32}, bias::Vector{Float32}=Float32[]; activation=:identity) -> BNNSFilter

Create a fully-connected layer computing `activation.(W * x .+ bias)` where `W` is
`out × in` (`BNNSFilterCreateLayerFullyConnected`). Apply it with
[`bnns_filter_apply!`](@ref) (single vector) or [`bnns_filter_apply_batch!`](@ref)
(an `in × batch` matrix). `activation` is one of `:identity`, `:relu`, `:sigmoid`,
`:tanh`, `:abs`, `:gelu`, `:softmax`. Pass an empty `bias` for no bias.

!!! warning
    Deprecated by Apple (macOS 15+); prefer the BNNS Graph API for new code.
"""
function bnns_fully_connected(W::Matrix{Float32}, bias::Vector{Float32} = Float32[];
                              activation::Symbol = :identity)
    out, inn = size(W)
    isempty(bias) || length(bias) == out ||
        throw(DimensionMismatch("bias length $(length(bias)) != out $out"))
    act = _activation_struct(activation)
    roots = Any[W]
    biasd = if isempty(bias)
        _vec_shape_desc(0)
    else
        push!(roots, bias); _vec_shape_desc(out; data = pointer(bias))
    end
    wd = _col_matrix_desc(out, inn; data = pointer(W))
    lp = Ref(LA.BNNSLayerParametersFullyConnected(_vec_shape_desc(inn), wd,
             _vec_shape_desc(out), biasd, act))
    fp = Ref(_default_fp())
    h = GC.@preserve W bias lp fp begin
        LA.BNNSFilterCreateLayerFullyConnected(
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersFullyConnected}, lp),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_filter_apply!(f::BNNSFilter, out::Array{Float32}, inp::Array{Float32}) -> out

Apply a single-input filter to `in`, writing `out` (`BNNSFilterApply`).
"""
function bnns_filter_apply!(f::BNNSFilter, out::Array{Float32}, inp::Array{Float32})
    GC.@preserve out inp begin
        _bnns_check(LA.BNNSFilterApply(f.handle, pointer(inp), pointer(out)), "BNNSFilterApply")
    end
    return out
end

"""
    bnns_filter_apply_batch!(f::BNNSFilter, out::Matrix{Float32}, inp::Matrix{Float32}) -> out

Apply a filter across a batch (`BNNSFilterApplyBatch`). `in` is `in_size × batch`,
`out` is `out_size × batch`; the per-sample strides are the column lengths.
"""
function bnns_filter_apply_batch!(f::BNNSFilter, out::Matrix{Float32}, inp::Matrix{Float32})
    batch = size(inp, 2)
    size(out, 2) == batch || throw(DimensionMismatch("batch mismatch"))
    GC.@preserve out inp begin
        _bnns_check(LA.BNNSFilterApplyBatch(f.handle, Csize_t(batch),
            pointer(inp), Csize_t(size(inp, 1)), pointer(out), Csize_t(size(out, 1))),
            "BNNSFilterApplyBatch")
    end
    return out
end

"""
    bnns_filter_apply_two_input!(f::BNNSFilter, out, inA, inB) -> out

Apply a two-input filter (e.g. an arithmetic layer) to `inA`, `inB`
(`BNNSFilterApplyTwoInput`).
"""
function bnns_filter_apply_two_input!(f::BNNSFilter, out::Array{Float32},
                                      inA::Array{Float32}, inB::Array{Float32})
    GC.@preserve out inA inB begin
        _bnns_check(LA.BNNSFilterApplyTwoInput(f.handle, pointer(inA), pointer(inB), pointer(out)),
                    "BNNSFilterApplyTwoInput")
    end
    return out
end

"""
    bnns_filter_apply_two_input_batch!(f::BNNSFilter, out, inA, inB) -> out

Batched two-input apply (`BNNSFilterApplyTwoInputBatch`); operands are
`feature × batch` matrices.
"""
function bnns_filter_apply_two_input_batch!(f::BNNSFilter, out::Matrix{Float32},
                                            inA::Matrix{Float32}, inB::Matrix{Float32})
    batch = size(inA, 2)
    GC.@preserve out inA inB begin
        _bnns_check(LA.BNNSFilterApplyTwoInputBatch(f.handle, Csize_t(batch),
            pointer(inA), Csize_t(size(inA, 1)), pointer(inB), Csize_t(size(inB, 1)),
            pointer(out), Csize_t(size(out, 1))), "BNNSFilterApplyTwoInputBatch")
    end
    return out
end

"""
    bnns_pooling_apply_batch!(f, out, in; indices=nothing) -> out

Batched pooling apply (`BNNSPoolingFilterApplyBatch`). Optional `indices`
(`Matrix{Csize_t}`) receives argmax positions for max-pooling.
"""
function bnns_pooling_apply_batch!(f::BNNSFilter, out::Matrix{Float32}, inp::Matrix{Float32};
                                   indices::Union{Nothing,Matrix{Csize_t}} = nothing)
    batch = size(inp, 2)
    GC.@preserve out inp indices begin
        iptr = indices === nothing ? Ptr{Csize_t}(C_NULL) : pointer(indices)
        istride = indices === nothing ? Csize_t(0) : Csize_t(size(indices, 1))
        _bnns_check(LA.BNNSPoolingFilterApplyBatch(f.handle, Csize_t(batch),
            pointer(inp), Csize_t(size(inp, 1)), pointer(out), Csize_t(size(out, 1)),
            iptr, istride), "BNNSPoolingFilterApplyBatch")
    end
    return out
end

"""
    bnns_pooling_apply_batch_ex!(f, out, in; indices=nothing, index_type=Int32) -> out

Extended batched pooling apply with a typed index buffer
(`BNNSPoolingFilterApplyBatchEx`).
"""
function bnns_pooling_apply_batch_ex!(f::BNNSFilter, out::Matrix{Float32}, inp::Matrix{Float32};
                                      indices::Union{Nothing,Matrix} = nothing,
                                      index_type::Type = Int32)
    batch = size(inp, 2)
    GC.@preserve out inp indices begin
        iptr = indices === nothing ? Ptr{Cvoid}(C_NULL) : Ptr{Cvoid}(pointer(indices))
        istride = indices === nothing ? Csize_t(0) : Csize_t(size(indices, 1))
        _bnns_check(LA.BNNSPoolingFilterApplyBatchEx(f.handle, Csize_t(batch),
            pointer(inp), Csize_t(size(inp, 1)), pointer(out), Csize_t(size(out, 1)),
            _bnns_dt(index_type), iptr, istride), "BNNSPoolingFilterApplyBatchEx")
    end
    return out
end

"""
    bnns_normalization_apply_batch!(f, out, in; training=false) -> out

Batched normalization apply (`BNNSNormalizationFilterApplyBatch`).
"""
function bnns_normalization_apply_batch!(f::BNNSFilter, out::Matrix{Float32},
                                         inp::Matrix{Float32}; training::Bool = false)
    batch = size(inp, 2)
    GC.@preserve out inp begin
        _bnns_check(LA.BNNSNormalizationFilterApplyBatch(f.handle, Csize_t(batch),
            pointer(inp), Csize_t(size(inp, 1)), pointer(out), Csize_t(size(out, 1)), training),
            "BNNSNormalizationFilterApplyBatch")
    end
    return out
end

"""
    bnns_arithmetic_apply_batch!(f, out, ins::Vector{Matrix{Float32}}) -> out

Batched arithmetic apply over `number_of_inputs` operands
(`BNNSArithmeticFilterApplyBatch`).
"""
function bnns_arithmetic_apply_batch!(f::BNNSFilter, out::Matrix{Float32},
                                      ins::Vector{Matrix{Float32}})
    batch = size(out, 2); nin = length(ins)
    inptrs = [Ptr{Cvoid}(pointer(x)) for x in ins]
    strides = [Csize_t(size(x, 1)) for x in ins]
    GC.@preserve out ins inptrs strides begin
        _bnns_check(LA.BNNSArithmeticFilterApplyBatch(f.handle, Csize_t(batch), Csize_t(nin),
            pointer(inptrs), pointer(strides), pointer(out), Csize_t(size(out, 1))),
            "BNNSArithmeticFilterApplyBatch")
    end
    return out
end

"""
    bnns_fused_apply_batch!(f, out, in; training=false) -> out

Batched fused-layer apply (`BNNSFusedFilterApplyBatch`).
"""
function bnns_fused_apply_batch!(f::BNNSFilter, out::Matrix{Float32}, inp::Matrix{Float32};
                                 training::Bool = false)
    batch = size(inp, 2)
    GC.@preserve out inp begin
        _bnns_check(LA.BNNSFusedFilterApplyBatch(f.handle, Csize_t(batch),
            pointer(inp), Csize_t(size(inp, 1)), pointer(out), Csize_t(size(out, 1)), training),
            "BNNSFusedFilterApplyBatch")
    end
    return out
end

"""
    bnns_fused_apply_multi_input_batch!(f, out, ins::Vector{Matrix{Float32}}; training=false) -> out

Batched multi-input fused-layer apply (`BNNSFusedFilterApplyMultiInputBatch`).
"""
function bnns_fused_apply_multi_input_batch!(f::BNNSFilter, out::Matrix{Float32},
                                             ins::Vector{Matrix{Float32}}; training::Bool = false)
    batch = size(out, 2); nin = length(ins)
    inptrs = [Ptr{Cvoid}(pointer(x)) for x in ins]
    strides = [Csize_t(size(x, 1)) for x in ins]
    GC.@preserve out ins inptrs strides begin
        _bnns_check(LA.BNNSFusedFilterApplyMultiInputBatch(f.handle, Csize_t(batch), Csize_t(nin),
            pointer(inptrs), pointer(strides), pointer(out), Csize_t(size(out, 1)), training),
            "BNNSFusedFilterApplyMultiInputBatch")
    end
    return out
end

# --- Thin create wrappers for the remaining layer types ----------------------
# Each takes the raw `LibAccelerate` parameter struct (whose embedded descriptors
# point at the caller's arrays) and returns a managed `BNNSFilter`. Pass the
# backing arrays via `roots` so they outlive the filter. These are deprecated by
# Apple in favour of the Graph API; they are provided for completeness.

for (jl, cfun, PT) in (
    (:bnns_create_layer_convolution, :BNNSFilterCreateLayerConvolution, :BNNSLayerParametersConvolution),
    (:bnns_create_layer_transposed_convolution, :BNNSFilterCreateLayerTransposedConvolution, :BNNSLayerParametersConvolution),
    (:bnns_create_layer_fully_connected, :BNNSFilterCreateLayerFullyConnected, :BNNSLayerParametersFullyConnected),
    (:bnns_create_layer_pooling, :BNNSFilterCreateLayerPooling, :BNNSLayerParametersPooling),
    (:bnns_create_layer_activation, :BNNSFilterCreateLayerActivation, :BNNSLayerParametersActivation),
    (:bnns_create_layer_arithmetic, :BNNSFilterCreateLayerArithmetic, :BNNSLayerParametersArithmetic),
    (:bnns_create_layer_permute, :BNNSFilterCreateLayerPermute, :BNNSLayerParametersPermute),
    (:bnns_create_layer_dropout, :BNNSFilterCreateLayerDropout, :BNNSLayerParametersDropout),
    (:bnns_create_layer_padding, :BNNSFilterCreateLayerPadding, :BNNSLayerParametersPadding),
    (:bnns_create_layer_broadcast_matmul, :BNNSFilterCreateLayerBroadcastMatMul, :BNNSLayerParametersBroadcastMatMul),
    (:bnns_create_layer_tensor_contraction, :BNNSFilterCreateLayerTensorContraction, :BNNSLayerParametersTensorContraction),
    (:bnns_create_layer_gram, :BNNSFilterCreateLayerGram, :BNNSLayerParametersGram),
    (:bnns_create_layer_resize, :BNNSFilterCreateLayerResize, :BNNSLayerParametersResize),
    (:bnns_create_layer_multihead_attention, :BNNSFilterCreateLayerMultiheadAttention, :BNNSLayerParametersMultiheadAttention),
    (:bnns_create_layer_reduction, :BNNSFilterCreateLayerReduction, :BNNSLayerParametersReduction),
    (:bnns_create_layer_embedding, :BNNSFilterCreateLayerEmbedding, :BNNSLayerParametersEmbedding),
)
    @eval begin
        """
            $($(string(jl)))(params::LibAccelerate.$($(string(PT))); roots=Any[]) -> BNNSFilter

        Thin wrapper over `$($(string(cfun)))`. Build `params` (embedding
        descriptors that point at your backing arrays) and pass those arrays as
        `roots` to keep them alive. Deprecated by Apple in favour of the Graph API.
        """
        function $jl(params::LA.$PT; roots::Vector{Any} = Any[])
            pref = Ref(params); fp = Ref(_default_fp())
            h = GC.@preserve pref fp roots begin
                LA.$cfun(Base.unsafe_convert(Ptr{LA.$PT}, pref),
                         Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
            end
            return BNNSFilter(h, roots)
        end
    end
end

const _NORM_TYPE = Dict(:batch => LA.BNNSBatchNorm, :instance => LA.BNNSInstanceNorm,
                        :layer => LA.BNNSLayerNorm, :group => LA.BNNSGroupNorm)

"""
    bnns_create_layer_normalization(norm_type::Symbol, params::LibAccelerate.BNNSLayerParametersNormalization; roots=Any[]) -> BNNSFilter

Thin wrapper over `BNNSFilterCreateLayerNormalization`. `norm_type` is `:batch`,
`:instance`, `:layer` or `:group`. Deprecated by Apple in favour of the Graph API.
"""
function bnns_create_layer_normalization(norm_type::Symbol,
                                         params::LA.BNNSLayerParametersNormalization;
                                         roots::Vector{Any} = Any[])
    haskey(_NORM_TYPE, norm_type) || throw(ArgumentError("BNNS: bad norm type $(repr(norm_type))"))
    pref = Ref(params); fp = Ref(_default_fp())
    h = GC.@preserve pref fp roots begin
        LA.BNNSFilterCreateLayerNormalization(LA.BNNSFilterType(_NORM_TYPE[norm_type]),
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersNormalization}, pref),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_create_layer_loss(params; roots=Any[]) -> BNNSFilter

Thin wrapper over `BNNSFilterCreateLayerLoss`. `params` is one of the loss
parameter structs (e.g. `BNNSLayerParametersLossSoftmaxCrossEntropy`), passed by
`Ref` as an opaque `void*`. Deprecated by Apple in favour of the Graph API.
"""
function bnns_create_layer_loss(params; roots::Vector{Any} = Any[])
    pref = Ref(params); fp = Ref(_default_fp())
    h = GC.@preserve pref fp roots begin
        LA.BNNSFilterCreateLayerLoss(Ptr{Cvoid}(Base.unsafe_convert(Ptr{typeof(params)}, pref)),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_create_fused_layer(filter_types::Vector, layer_param_ptrs::Vector{Ptr{Cvoid}}; roots=Any[]) -> BNNSFilter

Thin wrapper over `BNNSFilterCreateFusedLayer`, fusing `length(filter_types)`
sub-filters. Deprecated by Apple in favour of the Graph API.
"""
function bnns_create_fused_layer(filter_types::Vector{<:Integer},
                                 layer_param_ptrs::Vector{Ptr{Cvoid}};
                                 roots::Vector{Any} = Any[])
    n = length(filter_types)
    ftypes = LA.BNNSFilterType.(filter_types)
    fp = Ref(_default_fp())
    h = GC.@preserve ftypes layer_param_ptrs fp roots begin
        LA.BNNSFilterCreateFusedLayer(Csize_t(n), pointer(ftypes),
            pointer(layer_param_ptrs), Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

# --- Legacy (pre-NDArray) create wrappers ------------------------------------

"""
    bnns_create_convolution_layer(in_desc, out_desc, params::LibAccelerate.BNNSConvolutionLayerParameters; roots=Any[]) -> BNNSFilter

Thin wrapper over the legacy `BNNSFilterCreateConvolutionLayer` (image-stack
descriptors). Deprecated by Apple in favour of the Graph API.
"""
function bnns_create_convolution_layer(in_desc::LA.BNNSImageStackDescriptor,
                                       out_desc::LA.BNNSImageStackDescriptor,
                                       params::LA.BNNSConvolutionLayerParameters;
                                       roots::Vector{Any} = Any[])
    ir = Ref(in_desc); orf = Ref(out_desc); pr = Ref(params); fp = Ref(_default_fp())
    h = GC.@preserve ir orf pr fp roots begin
        LA.BNNSFilterCreateConvolutionLayer(
            Base.unsafe_convert(Ptr{LA.BNNSImageStackDescriptor}, ir),
            Base.unsafe_convert(Ptr{LA.BNNSImageStackDescriptor}, orf),
            Base.unsafe_convert(Ptr{LA.BNNSConvolutionLayerParameters}, pr),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_create_fully_connected_layer(in_desc, out_desc, params::LibAccelerate.BNNSFullyConnectedLayerParameters; roots=Any[]) -> BNNSFilter

Thin wrapper over the legacy `BNNSFilterCreateFullyConnectedLayer` (vector
descriptors). Deprecated by Apple in favour of the Graph API.
"""
function bnns_create_fully_connected_layer(in_desc::LA.BNNSVectorDescriptor,
                                           out_desc::LA.BNNSVectorDescriptor,
                                           params::LA.BNNSFullyConnectedLayerParameters;
                                           roots::Vector{Any} = Any[])
    ir = Ref(in_desc); orf = Ref(out_desc); pr = Ref(params); fp = Ref(_default_fp())
    h = GC.@preserve ir orf pr fp roots begin
        LA.BNNSFilterCreateFullyConnectedLayer(
            Base.unsafe_convert(Ptr{LA.BNNSVectorDescriptor}, ir),
            Base.unsafe_convert(Ptr{LA.BNNSVectorDescriptor}, orf),
            Base.unsafe_convert(Ptr{LA.BNNSFullyConnectedLayerParameters}, pr),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_create_pooling_layer(in_desc, out_desc, params::LibAccelerate.BNNSPoolingLayerParameters; roots=Any[]) -> BNNSFilter

Thin wrapper over the legacy `BNNSFilterCreatePoolingLayer`. Deprecated by Apple
in favour of the Graph API.
"""
function bnns_create_pooling_layer(in_desc::LA.BNNSImageStackDescriptor,
                                   out_desc::LA.BNNSImageStackDescriptor,
                                   params::LA.BNNSPoolingLayerParameters;
                                   roots::Vector{Any} = Any[])
    ir = Ref(in_desc); orf = Ref(out_desc); pr = Ref(params); fp = Ref(_default_fp())
    h = GC.@preserve ir orf pr fp roots begin
        LA.BNNSFilterCreatePoolingLayer(
            Base.unsafe_convert(Ptr{LA.BNNSImageStackDescriptor}, ir),
            Base.unsafe_convert(Ptr{LA.BNNSImageStackDescriptor}, orf),
            Base.unsafe_convert(Ptr{LA.BNNSPoolingLayerParameters}, pr),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_create_vector_activation_layer(in_desc, out_desc, activation::Symbol; alpha=0f0, beta=0f0, roots=Any[]) -> BNNSFilter

Thin wrapper over the legacy `BNNSFilterCreateVectorActivationLayer`. Deprecated
by Apple in favour of the Graph API.
"""
function bnns_create_vector_activation_layer(in_desc::LA.BNNSVectorDescriptor,
                                             out_desc::LA.BNNSVectorDescriptor,
                                             activation::Symbol; alpha = 0.0f0, beta = 0.0f0,
                                             roots::Vector{Any} = Any[])
    ir = Ref(in_desc); orf = Ref(out_desc)
    ar = Ref(_activation_struct(activation; alpha, beta)); fp = Ref(_default_fp())
    h = GC.@preserve ir orf ar fp roots begin
        LA.BNNSFilterCreateVectorActivationLayer(
            Base.unsafe_convert(Ptr{LA.BNNSVectorDescriptor}, ir),
            Base.unsafe_convert(Ptr{LA.BNNSVectorDescriptor}, orf),
            Base.unsafe_convert(Ptr{LA.BNNSActivation}, ar),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp))
    end
    return BNNSFilter(h, roots)
end

"""
    bnns_get_pointer(f::BNNSFilter, which::Symbol=:alpha) -> BNNSNDArrayDescriptor

Retrieve a descriptor for an internal filter buffer (`BNNSGetPointer`); `which`
is `:alpha` or `:beta`.
"""
function bnns_get_pointer(f::BNNSFilter, which::Symbol = :alpha)
    spec = which === :alpha ? LA.BNNSPointerSpecifierAlpha :
           which === :beta ? LA.BNNSPointerSpecifierBeta :
           throw(ArgumentError("BNNS: which must be :alpha or :beta"))
    return LA.BNNSGetPointer(f.handle, LA.BNNSPointerSpecifier(spec))
end

# =============================================================================
# Remaining DirectApply kernels
# =============================================================================

"""
    bnns_in_topk(input::Array{Float32}, targets::Array{Int32}, K; dim=1) -> Array{Bool}

For each batch column, test whether the `targets` class index is among the top-`K`
scores of `input` along Julia dimension `dim` (`BNNSDirectApplyInTopK`).
"""
function bnns_in_topk(input::Array{Float32}, targets::Array{Int32}, K::Integer; dim::Integer = 1)
    batch = length(targets)
    out = Array{Bool}(undef, size(targets))
    di = _desc(input); dt = _desc(targets); do_ = _desc(out)
    GC.@preserve input targets out begin
        _bnns_check(LA.BNNSDirectApplyInTopK(Csize_t(K), Csize_t(dim - 1), Csize_t(1),
            Ref(di), Csize_t(0), Ref(dt), Csize_t(0), Ref(do_), Csize_t(0), C_NULL),
            "BNNSDirectApplyInTopK")
    end
    return out
end

"""
    bnns_activation_batch!(func::Symbol, out::Array{Float32}, inp::Array{Float32}; alpha=0f0, beta=0f0) -> out

Apply a pointwise activation over a batch directly (`BNNSDirectApplyActivationBatch`),
without an explicit filter handle. `func` is any activation in
[`bnns_fully_connected`](@ref)'s list.
"""
function bnns_activation_batch!(func::Symbol, out::Array{Float32}, inp::Array{Float32};
                                alpha = 0.0f0, beta = 0.0f0)
    di = _desc(inp); do_ = _desc(out)
    act = _activation_struct(func; alpha, beta)
    lp = Ref(LA.BNNSLayerParametersActivation(di, do_, act, UInt32(0)))
    fp = Ref(_default_fp())
    GC.@preserve inp out lp fp begin
        _bnns_check(LA.BNNSDirectApplyActivationBatch(
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersActivation}, lp),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp),
            Csize_t(1), Csize_t(0), Csize_t(0)), "BNNSDirectApplyActivationBatch")
    end
    return out
end

const _QUANT_FN = Dict(:quantize => LA.BNNSQuantizerFunctionQuantize,
                       :dequantize => LA.BNNSQuantizerFunctionDequantize)

"""
    bnns_quantizer!(func::Symbol, out::Array, inp::Array, scale::Array{Float32}, bias::Array{Float32}; axis_mask=0) -> out

Quantize or dequantize `in` into `out` with per-`axis_mask` `scale`/`bias`
(`BNNSDirectApplyQuantizer`). `func` is `:quantize` or `:dequantize`.
"""
function bnns_quantizer!(func::Symbol, out::Array, inp::Array, scale::Array{Float32},
                         bias::Array{Float32}; axis_mask::Integer = 0)
    haskey(_QUANT_FN, func) || throw(ArgumentError("BNNS: func must be :quantize/:dequantize"))
    di = _desc(inp); do_ = _desc(out); ds = _desc(scale); db = _desc(bias)
    lp = Ref(LA.BNNSLayerParametersQuantization(Csize_t(axis_mask),
        LA.BNNSQuantizerFunction(_QUANT_FN[func]), di, do_, ds, db))
    fp = Ref(_default_fp())
    GC.@preserve inp out scale bias lp fp begin
        _bnns_check(LA.BNNSDirectApplyQuantizer(
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersQuantization}, lp),
            Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp),
            Csize_t(1), Csize_t(0), Csize_t(0)), "BNNSDirectApplyQuantizer")
    end
    return out
end

# =============================================================================
# Optimizer step (stateless kernel over parameter/gradient/accumulator lists)
# =============================================================================

"""
    bnns_optimizer_step_sgd!(params, grads, accumulators, fields::BNNSOptimizerSGDMomentumFields) -> params

Apply one SGD-with-momentum update (`BNNSOptimizerStep`,
`BNNSOptimizerFunctionSGDMomentum`) across the parallel lists of parameter,
gradient and momentum-accumulator `Array{Float32}`s, updating each in place.
`fields` is a raw `LibAccelerate.BNNSOptimizerSGDMomentumFields`.
"""
function bnns_optimizer_step_sgd!(params::Vector{<:Array{Float32}},
                                  grads::Vector{<:Array{Float32}},
                                  accumulators::Vector{<:Array{Float32}},
                                  fields::LA.BNNSOptimizerSGDMomentumFields)
    n = length(params)
    (length(grads) == n && length(accumulators) == n) ||
        throw(DimensionMismatch("params/grads/accumulators length mismatch"))
    pd = [Ref(_desc(x)) for x in params]
    gd = [Ref(_desc(x)) for x in grads]
    ad = [Ref(_desc(x)) for x in accumulators]
    fr = Ref(fields); fp = Ref(_default_fp())
    GC.@preserve params grads accumulators pd gd ad fr fp begin
        pp = [Base.unsafe_convert(Ptr{BNNSNDArrayDescriptor}, r) for r in pd]
        gp = [Base.unsafe_convert(Ptr{BNNSNDArrayDescriptor}, r) for r in gd]
        ap = [Base.unsafe_convert(Ptr{BNNSNDArrayDescriptor}, r) for r in ad]
        GC.@preserve pp gp ap begin
            _bnns_check(LA.BNNSOptimizerStep(
                LA.BNNSOptimizerFunction(LA.BNNSOptimizerFunctionSGDMomentum),
                Ptr{Cvoid}(Base.unsafe_convert(Ptr{LA.BNNSOptimizerSGDMomentumFields}, fr)),
                Csize_t(n), pointer(pp), pointer(gp), pointer(ap),
                Base.unsafe_convert(Ptr{LA.BNNSFilterParameters}, fp)), "BNNSOptimizerStep")
        end
    end
    return params
end

# =============================================================================
# Additional forward / gradient / sparsify kernels (thin wrappers).
#
# These round out the coverage of the deprecated filter API. Their FFI signatures
# are exact, but (unlike the verified helpers above) they are not numerically
# cross-validated here — they exist so callers who need the full BNNS surface do
# not have to drop to the raw layer. All are deprecated by Apple in favour of the
# Graph API.
# =============================================================================

_descptr(r) = Base.unsafe_convert(Ptr{BNNSNDArrayDescriptor}, r)

"""
    bnns_crop_resize!(out, input, roi, params::LibAccelerate.BNNSLayerParametersCropResize) -> out

Crop-and-resize regions of `input` selected by `roi` into `out` (`BNNSCropResize`).
"""
function bnns_crop_resize!(out::Array{Float32}, input::Array{Float32}, roi::Array{Float32},
                           params::LA.BNNSLayerParametersCropResize)
    di = Ref(_desc(input)); dr = Ref(_desc(roi)); do_ = Ref(_desc(out)); pr = Ref(params)
    GC.@preserve input roi out pr begin
        _bnns_check(LA.BNNSCropResize(
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersCropResize}, pr),
            _descptr(di), _descptr(dr), _descptr(do_), C_NULL), "BNNSCropResize")
    end
    return out
end

"""
    bnns_crop_resize_backward!(in_delta, roi, out_delta, params) -> in_delta

Backward pass of [`bnns_crop_resize!`](@ref) (`BNNSCropResizeBackward`).
"""
function bnns_crop_resize_backward!(in_delta::Array{Float32}, roi::Array{Float32},
                                    out_delta::Array{Float32},
                                    params::LA.BNNSLayerParametersCropResize)
    did = Ref(_desc(in_delta)); dr = Ref(_desc(roi)); dod = Ref(_desc(out_delta)); pr = Ref(params)
    GC.@preserve in_delta roi out_delta pr begin
        _bnns_check(LA.BNNSCropResizeBackward(
            Base.unsafe_convert(Ptr{LA.BNNSLayerParametersCropResize}, pr),
            _descptr(did), _descptr(dr), _descptr(dod), C_NULL), "BNNSCropResizeBackward")
    end
    return in_delta
end

"""
    bnns_compute_norm_backward!(in_delta, out_delta, input, out; axis_flags=0) -> in_delta

Backward pass of [`bnns_compute_norm`](@ref) (`BNNSComputeNormBackward`, L2 norm).
"""
function bnns_compute_norm_backward!(in_delta::Array{Float32}, out_delta::Array{Float32},
                                     input::Array{Float32}, out::Array{Float32};
                                     axis_flags::Integer = 0)
    did = Ref(_desc(in_delta)); dod = Ref(_desc(out_delta))
    GC.@preserve in_delta out_delta input out did dod begin
        _bnns_check(LA.BNNSComputeNormBackward(pointer(input), _descptr(did),
            pointer(out), _descptr(dod), LA.BNNSNormType(LA.BNNSL2Norm), UInt32(axis_flags)),
            "BNNSComputeNormBackward")
    end
    return in_delta
end

"""
    bnns_loss_apply_batch!(f::BNNSFilter, out, input, labels; weights=Float32[], in_delta=nothing) -> out

Apply a loss filter over a batch (`BNNSLossFilterApplyBatch`), writing the loss to
`out` and, optionally, the input gradient into `in_delta`.
"""
function bnns_loss_apply_batch!(f::BNNSFilter, out::Array{Float32}, input::Matrix{Float32},
                                labels::Matrix{Float32}; weights::Vector{Float32} = Float32[],
                                in_delta::Union{Nothing,Array{Float32}} = nothing)
    batch = size(input, 2)
    idd = in_delta === nothing ? Ref(_empty_desc()) : Ref(_desc(in_delta))
    GC.@preserve out input labels weights in_delta idd begin
        wptr = isempty(weights) ? Ptr{Cvoid}(C_NULL) : Ptr{Cvoid}(pointer(weights))
        _bnns_check(LA.BNNSLossFilterApplyBatch(f.handle, Csize_t(batch),
            pointer(input), Csize_t(size(input, 1)), pointer(labels), Csize_t(size(labels, 1)),
            wptr, Csize_t(length(weights)), pointer(out),
            in_delta === nothing ? Ptr{BNNSNDArrayDescriptor}(C_NULL) : _descptr(idd),
            Csize_t(in_delta === nothing ? 0 : size(in_delta, 1))), "BNNSLossFilterApplyBatch")
    end
    return out
end

# Generic backward-apply wrappers. `in`/`out` are raw data buffers; the `*_delta`
# operands are gradient tensors described by `BNNSNDArrayDescriptor`s built from
# the passed arrays. Pass `nothing` for optional weight/bias gradients.
_optdescptr(r, arr) = arr === nothing ? Ptr{BNNSNDArrayDescriptor}(C_NULL) : _descptr(r)

"""
    bnns_filter_apply_backward_batch!(f, input, in_delta, out, out_delta; weights_delta=nothing, bias_delta=nothing) -> in_delta

Backward pass of a single-input filter (`BNNSFilterApplyBackwardBatch`).
"""
function bnns_filter_apply_backward_batch!(f::BNNSFilter, input::Matrix{Float32},
        in_delta::Matrix{Float32}, out::Matrix{Float32}, out_delta::Matrix{Float32};
        weights_delta::Union{Nothing,Array{Float32}} = nothing,
        bias_delta::Union{Nothing,Array{Float32}} = nothing)
    batch = size(input, 2)
    idd = Ref(_desc(in_delta)); odd = Ref(_desc(out_delta))
    wdd = Ref(weights_delta === nothing ? _empty_desc() : _desc(weights_delta))
    bdd = Ref(bias_delta === nothing ? _empty_desc() : _desc(bias_delta))
    GC.@preserve input in_delta out out_delta weights_delta bias_delta idd odd wdd bdd begin
        _bnns_check(LA.BNNSFilterApplyBackwardBatch(f.handle, Csize_t(batch),
            pointer(input), Csize_t(size(input, 1)), _descptr(idd), Csize_t(size(in_delta, 1)),
            pointer(out), Csize_t(size(out, 1)), _descptr(odd), Csize_t(size(out_delta, 1)),
            _optdescptr(wdd, weights_delta), _optdescptr(bdd, bias_delta)),
            "BNNSFilterApplyBackwardBatch")
    end
    return in_delta
end

"""
    bnns_normalization_apply_backward_batch!(f, in_delta, out, out_delta; beta_delta=nothing, gamma_delta=nothing) -> in_delta

Backward pass of a normalization filter (`BNNSNormalizationFilterApplyBackwardBatch`).
"""
function bnns_normalization_apply_backward_batch!(f::BNNSFilter, in_delta::Matrix{Float32},
        out::Matrix{Float32}, out_delta::Matrix{Float32};
        beta_delta::Union{Nothing,Array{Float32}} = nothing,
        gamma_delta::Union{Nothing,Array{Float32}} = nothing)
    batch = size(out, 2)
    idd = Ref(_desc(in_delta)); odd = Ref(_desc(out_delta))
    bdd = Ref(beta_delta === nothing ? _empty_desc() : _desc(beta_delta))
    gdd = Ref(gamma_delta === nothing ? _empty_desc() : _desc(gamma_delta))
    GC.@preserve in_delta out out_delta beta_delta gamma_delta idd odd bdd gdd begin
        _bnns_check(LA.BNNSNormalizationFilterApplyBackwardBatch(f.handle, Csize_t(batch),
            _descptr(idd), Csize_t(size(in_delta, 1)), pointer(out), Csize_t(size(out, 1)),
            _descptr(odd), Csize_t(size(out_delta, 1)),
            _optdescptr(bdd, beta_delta), _optdescptr(gdd, gamma_delta)),
            "BNNSNormalizationFilterApplyBackwardBatch")
    end
    return in_delta
end

"""
    bnns_permute_apply_backward_batch!(f, in_delta, out_delta) -> in_delta

Backward pass of a permute filter (`BNNSPermuteFilterApplyBackwardBatch`).
"""
function bnns_permute_apply_backward_batch!(f::BNNSFilter, in_delta::Matrix{Float32},
                                            out_delta::Matrix{Float32})
    batch = size(out_delta, 2)
    idd = Ref(_desc(in_delta)); odd = Ref(_desc(out_delta))
    GC.@preserve in_delta out_delta idd odd begin
        _bnns_check(LA.BNNSPermuteFilterApplyBackwardBatch(f.handle, Csize_t(batch),
            _descptr(idd), Csize_t(size(in_delta, 1)), _descptr(odd), Csize_t(size(out_delta, 1))),
            "BNNSPermuteFilterApplyBackwardBatch")
    end
    return in_delta
end

"""
    bnns_pooling_apply_backward_batch!(f, input, in_delta, out, out_delta; bias_delta=nothing, indices=nothing) -> in_delta

Backward pass of a pooling filter (`BNNSPoolingFilterApplyBackwardBatch`).
"""
function bnns_pooling_apply_backward_batch!(f::BNNSFilter, input::Matrix{Float32},
        in_delta::Matrix{Float32}, out::Matrix{Float32}, out_delta::Matrix{Float32};
        bias_delta::Union{Nothing,Array{Float32}} = nothing,
        indices::Union{Nothing,Matrix{Csize_t}} = nothing)
    batch = size(input, 2)
    idd = Ref(_desc(in_delta)); odd = Ref(_desc(out_delta))
    bdd = Ref(bias_delta === nothing ? _empty_desc() : _desc(bias_delta))
    GC.@preserve input in_delta out out_delta bias_delta indices idd odd bdd begin
        iptr = indices === nothing ? Ptr{Csize_t}(C_NULL) : pointer(indices)
        istride = indices === nothing ? Csize_t(0) : Csize_t(size(indices, 1))
        _bnns_check(LA.BNNSPoolingFilterApplyBackwardBatch(f.handle, Csize_t(batch),
            pointer(input), Csize_t(size(input, 1)), _descptr(idd), Csize_t(size(in_delta, 1)),
            pointer(out), Csize_t(size(out, 1)), _descptr(odd), Csize_t(size(out_delta, 1)),
            _optdescptr(bdd, bias_delta), iptr, istride), "BNNSPoolingFilterApplyBackwardBatch")
    end
    return in_delta
end
