# BNNS (Basic Neural Network Subroutines) idiomatic wrappers.
#
# This module wraps the Apple-CURRENT (non-deprecated) slice of BNNS exposed by
# the raw `LibAccelerate` layer: dense array descriptors (`BNNSArray`), the
# stateless tensor ops that remain current (`bnns_transpose`, `bnns_copy!`), the
# DirectApply reduction / top-k kernels (`bnns_reduce`, `bnns_topk`,
# `bnns_in_topk`), random generation, nearest neighbors, layout/size queries, and
# the full BNNS Graph compiler pipeline, including an end-to-end inference path
# (`bnns_graph_arguments`, `bnns_graph_run`, `bnns_graph_run!`). The numerically
# important helpers are cross-validated against plain-Julia references.
#
# Element types: `Float32` everywhere, plus `Float16` and the integer / `Bool`
# types wherever the kernel was verified to honour them at runtime (each
# docstring lists its own set — BNNS returns status 0 with garbage for some
# unsupported types, so the sets are enforced here by dispatch).
#
# NOTE: Every classic and DirectApply BNNS entry point that Apple DEPRECATED in
# macOS 15.0 / iOS 18.0 is intentionally NOT wrapped here — target the BNNS Graph
# API (below) instead. That excluded set includes `BNNSMatMul`, `BNNSTile`,
# `BNNSGather`/`BNNSScatter` (and their ND forms), `BNNSCompareTensor`,
# `BNNSBandPart`, `BNNSShuffle`, the clip / `BNNSComputeNorm` family,
# `BNNSOptimizerStep`, the `BNNSDirectApply{ActivationBatch,BroadcastMatMul,
# Quantizer}` kernels, and the whole `BNNSFilter*` layer create/apply API.
#
# Apple docs: https://developer.apple.com/documentation/accelerate/bnns

using .LibAccelerate:
    BNNSNDArrayDescriptor, BNNSDataType, BNNSDataLayout,
    BNNSDataTypeFloat16, BNNSDataTypeFloat32, BNNSDataTypeInt32,
    BNNSDataLayoutVector, BNNSDataLayoutColumnMajorMatrix

# Alias the raw layer so the many BNNS enum values, structs and ccall wrappers can
# be reached without an unwieldy explicit import list. All new wrappers below go
# through `LA.<name>`.
import .LibAccelerate as LA

# --- Data-type mapping -------------------------------------------------------

# Map a Julia element type to the corresponding BNNS data-type enum value.
# The struct fields are typed `BNNSDataType` (a `UInt32` alias), so the enum
# value is converted to its underlying integer code.
bnns_data_type(::Type{Float16}) = BNNSDataType(BNNSDataTypeFloat16)
bnns_data_type(::Type{Float32}) = BNNSDataType(BNNSDataTypeFloat32)
bnns_data_type(::Type{Int32})   = BNNSDataType(BNNSDataTypeInt32)
bnns_data_type(::Type{T}) where {T} =
    throw(ArgumentError("BNNS: unsupported element type $T (supported: Float16, Float32, Int32)"))

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

Only `Float16`, `Float32` and `Int32` dense, contiguous arrays are supported
here; other element types or strided/transposed arrays should use the raw
`LibAccelerate` layer directly.
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
const _BNNS_DTYPES = (
    Float16 => LA.BNNSDataTypeFloat16, Float32 => LA.BNNSDataTypeFloat32,
    Int8 => LA.BNNSDataTypeInt8, Int16 => LA.BNNSDataTypeInt16,
    Int32 => LA.BNNSDataTypeInt32, Int64 => LA.BNNSDataTypeInt64,
    UInt8 => LA.BNNSDataTypeUInt8, UInt16 => LA.BNNSDataTypeUInt16,
    UInt32 => LA.BNNSDataTypeUInt32, UInt64 => LA.BNNSDataTypeUInt64,
    Bool => LA.BNNSDataTypeBoolean,
)
for (T, dt) in _BNNS_DTYPES
    @eval _bnns_dt(::Type{$T}) = BNNSDataType($dt)
end
_bnns_dt(::Type{T}) where {T} =
    throw(ArgumentError("BNNS: unsupported element type $T"))

# Inverse of `_bnns_dt`: the Julia element type for a BNNS data-type code (as
# reported by graph introspection). Sub-byte / indexed / bfloat16 codes have no
# Julia `Array` element type and throw.
function _julia_type(dt)
    for (T, code) in _BNNS_DTYPES
        UInt32(code) == UInt32(dt) && return T
    end
    throw(ArgumentError("BNNS: data type code $(repr(UInt32(dt))) has no Julia array element type"))
end

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
    bnns_transpose(A::Array, dim0, dim1) -> Array

Swap Julia dimensions `dim0` and `dim1` of `A` (1-based) via `BNNSTranspose`,
equivalent to a `permutedims` that exchanges those two axes. Works for every
element type BNNS can describe: `Float16`, `Float32`, `Int8`–`Int64`,
`UInt8`–`UInt64` and `Bool`.
"""
function bnns_transpose(A::Array{T}, dim0::Integer, dim1::Integer) where {T}
    (1 <= dim0 <= ndims(A) && 1 <= dim1 <= ndims(A)) ||
        throw(ArgumentError("BNNS: dims ($dim0, $dim1) out of range for a $(ndims(A))-D array"))
    p = collect(1:ndims(A)); p[dim0], p[dim1] = p[dim1], p[dim0]
    O = Array{T}(undef, ntuple(i -> size(A, p[i]), ndims(A)))
    di = _desc(A); do_ = _desc(O)
    GC.@preserve A O begin
        _bnns_check(LA.BNNSTranspose(Ref(do_), Ref(di),
                    Csize_t(dim0 - 1), Csize_t(dim1 - 1), C_NULL), "BNNSTranspose")
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

const _BNNS_COPY_CONVERSIONS = Set{Tuple{DataType,DataType}}([
    (Float16, Float32), (Float32, Float16), (Float32, Int32),
    (Int8, Float32), (Int16, Float32), (Int32, Float32), (UInt8, Float32),
])

"""
    bnns_copy!(dest::Array, src::Array) -> dest

Copy `src` into the equally-sized `dest` via `BNNSCopy`. For equal element types
this is a plain element copy.

`dest` and `src` may have **different element types**, in which case BNNS
converts: `Float16 ↔ Float32` (the usual way to move data in and out of a
half-precision graph), `Int8`/`Int16`/`Int32`/`UInt8` `→ Float32`, and
`Float32 → Int32`. Any other pair throws an `ArgumentError` without calling BNNS.
"""
function bnns_copy!(dest::Array{D}, src::Array{S}) where {D,S}
    # Only conversions verified on every supported platform are let through. An
    # unimplemented pair is not a clean error everywhere: Float32 -> Int8 returns a
    # status on Apple silicon but aborts the process inside BNNSCopy on Intel macOS 15.
    (D === S || (S, D) in _BNNS_COPY_CONVERSIONS) || throw(ArgumentError(
        "bnns_copy!: BNNS cannot convert $S to $D; supported conversions are " *
        "Float16 <-> Float32, Int8/Int16/Int32/UInt8 -> Float32 and Float32 -> Int32"))
    # BNNSCopy does NOT broadcast: given a smaller `src` it returns status 0 and
    # leaves the rest of `dest` unwritten, so equal shapes are required here.
    size(dest) == size(src) || throw(DimensionMismatch(
        "bnns_copy!: dest has size $(size(dest)), src has size $(size(src))"))
    ds = _desc(src); dd = _desc(dest)
    GC.@preserve src dest begin
        _bnns_check(LA.BNNSCopy(Ref(dd), Ref(ds), C_NULL), "BNNSCopy")
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

# Element types each kernel was verified to compute correctly (macOS 26). The
# reduction kernel returns status 0 but wrong values for the other integer
# widths, so the restriction is enforced by dispatch rather than by status.
const _BNNSTopKTypes = Union{Float32,Float16,Int8,Int16,Int32,UInt8,UInt16}
const _BNNSReduceTypes = Union{Float32,Float16,Int32}

"""
    bnns_topk(input::Array, K; dim=1) -> (values, indices)

Top-`K` values and their 0-based indices along Julia dimension `dim` via
`BNNSDirectApplyTopK`. `values` has the element type of `input`, `indices` is
`Int32`. Comparable to `sort`-based `partialsortperm` per slice.

Supported element types: `Float32`, `Float16`, `Int8`, `Int16`, `Int32`,
`UInt8`, `UInt16` (BNNS rejects the wider integer types).
"""
function bnns_topk(input::Array{T}, K::Integer; dim::Integer = 1) where {T<:_BNNSTopKTypes}
    1 <= dim <= ndims(input) || throw(ArgumentError("BNNS: dim $dim out of range for a $(ndims(input))-D array"))
    1 <= K <= size(input, dim) || throw(ArgumentError("BNNS: K = $K out of range 1:$(size(input, dim))"))
    osz = collect(size(input)); osz[dim] = K
    vals = Array{T}(undef, osz...)
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
    bnns_reduce(func::Symbol, input::Array; dim=1) -> Array

Reduce `input` along Julia dimension `dim` with `func`
(`:sum`, `:mean`, `:max`, `:min`, `:sumsquare`, `:l1`, `:l2`, `:product`,
`:logsumexp`) via `BNNSDirectApplyReduction`. The reduced axis collapses to
length 1 and the result has the element type of `input`.

Supported element types: `Float32`, `Float16` (computed *in* half precision, so
sums saturate at `floatmax(Float16)` = 65504) and `Int32`. For `Int32` only the
reductions that are exact in integers are offered (`:sum`, `:max`, `:min`,
`:sumsquare`, `:l1`, `:product`); `:mean`, `:l2` and `:logsumexp` throw.
"""
function bnns_reduce(func::Symbol, input::Array{T}; dim::Integer = 1) where {T<:_BNNSReduceTypes}
    haskey(_BNNS_REDUCE, func) || throw(ArgumentError("BNNS: unsupported reduce $(repr(func))"))
    (T === Int32 && func in (:mean, :l2, :logsumexp)) && throw(ArgumentError(
        "BNNS: reduce $(repr(func)) is not exact for Int32 input; convert to Float32 first"))
    1 <= dim <= ndims(input) || throw(ArgumentError("BNNS: dim $dim out of range for a $(ndims(input))-D array"))
    osz = collect(size(input)); osz[dim] = 1
    O = zeros(T, osz...)
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

const _BNNSIntTypes = Union{Int8,Int16,Int32,Int64,UInt8,UInt16,UInt32,UInt64}

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
    bnns_random_fill_uniform!(g::BNNSRandomGenerator, A::Array, lo=0f0, hi=1f0) -> A

Fill `A` (`Float32` or `Float16`) with i.i.d. uniform samples on `[lo, hi)`
(`BNNSRandomFillUniformFloat`). For `Float16` the samples are rounded to half
precision, so a value can round up to exactly `hi`.
"""
function bnns_random_fill_uniform!(g::BNNSRandomGenerator, A::Array{<:Union{Float32,Float16}},
                                   lo::Real = 0.0f0, hi::Real = 1.0f0)
    d = _desc(A)
    GC.@preserve A begin
        _bnns_check(LA.BNNSRandomFillUniformFloat(g.handle, Ref(d), Float32(lo), Float32(hi)),
                    "BNNSRandomFillUniformFloat")
    end
    return A
end

"""
    bnns_random_fill_uniform_int!(g::BNNSRandomGenerator, A::Array{<:Integer}, lo, hi) -> A

Fill integer array `A` with i.i.d. uniform samples on the half-open range
`[lo, hi)` (`BNNSRandomFillUniformInt`). `A` may be `Int8`, `Int16`, `Int32`,
`Int64`, `UInt8`, `UInt16`, `UInt32` or `UInt64`; the range must fit the
element type.
"""
function bnns_random_fill_uniform_int!(g::BNNSRandomGenerator, A::Array{T},
                                       lo::Integer, hi::Integer) where {T<:_BNNSIntTypes}
    lo < hi || throw(ArgumentError("BNNS: need lo < hi, got [$lo, $hi)"))
    (typemin(T) <= lo && hi - 1 <= typemax(T)) ||
        throw(ArgumentError("BNNS: range [$lo, $hi) does not fit element type $T"))
    d = _desc(A)
    GC.@preserve A begin
        _bnns_check(LA.BNNSRandomFillUniformInt(g.handle, Ref(d), Int64(lo), Int64(hi)),
                    "BNNSRandomFillUniformInt")
    end
    return A
end

"""
    bnns_random_fill_normal!(g::BNNSRandomGenerator, A::Array, mean=0f0, stddev=1f0) -> A

Fill `A` (`Float32` or `Float16`) with i.i.d. Gaussian samples
(`BNNSRandomFillNormalFloat`).
"""
function bnns_random_fill_normal!(g::BNNSRandomGenerator, A::Array{<:Union{Float32,Float16}},
                                  mean::Real = 0.0f0, stddev::Real = 1.0f0)
    d = _desc(A)
    GC.@preserve A begin
        _bnns_check(LA.BNNSRandomFillNormalFloat(g.handle, Ref(d), Float32(mean), Float32(stddev)),
                    "BNNSRandomFillNormalFloat")
    end
    return A
end

"""
    bnns_random_fill_categorical!(g::BNNSRandomGenerator, out::Array{T}, probs::Array{T}; log_probs=false) -> out

Draw categorical samples (0-based category indices, stored as floating point)
into `out` using per-category weights `probs` (`BNNSRandomFillCategoricalFloat`).
Pass `log_probs=true` if `probs` holds log probabilities. `T` is `Float32` or
`Float16`; `out` and `probs` must share it (BNNS silently mis-samples mixed
precisions).
"""
function bnns_random_fill_categorical!(g::BNNSRandomGenerator, out::Array{T},
                                       probs::Array{T}; log_probs::Bool = false) where {T<:Union{Float32,Float16}}
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
# Lifecycle: build compile options -> `BNNSGraphCompileFromFile` a compiled Core
# ML model (`.mlmodelc`, i.e. a directory holding a MIL program `model.mil` plus
# its weights) into a `BNNSGraph` -> `BNNSGraphContextMake` an executable context
# -> introspect (`bnns_graph_arguments`) -> run (`bnns_graph_run`/`run!`).
# Compiling requires an on-disk model (there is no in-memory graph builder in this
# API); the tests exercise the whole pipeline with small hand-written MIL programs.
# =============================================================================

# --- `_v2` entry points --------------------------------------------------------
# bnns_graph.h redirects several functions to versioned symbols with
# `__asm__("_<name>_v2")`. The un-suffixed symbols still exported by libBNNS are
# the pre-release ABI (different argument lists): calling them with the header's
# signature segfaults. Clang.jl does not see asm labels, so the generated raw
# layer binds the un-suffixed names; the correct `_v2` symbols are bound here with
# the header's exact signatures. `bnns_graph_shape_t` is likewise re-declared: the
# header's field order is `{ size_t rank; uint64_t *shape; }`.
struct _BNNSGraphShape
    rank::Csize_t
    shape::Ptr{UInt64}
end

_graph_compile_v2(filename, func, opts) =
    @ccall LA.libacc.BNNSGraphCompileFromFile_v2(filename::Ptr{Cchar}, func::Ptr{Cchar},
        opts::LA.bnns_graph_compile_options_t)::LA.bnns_graph_t
_graph_input_names_v2(g, func, n, names) =
    @ccall LA.libacc.BNNSGraphGetInputNames_v2(g::LA.bnns_graph_t, func::Ptr{Cchar},
        n::Csize_t, names::Ptr{Ptr{Cchar}})::Cint
_graph_output_names_v2(g, func, n, names) =
    @ccall LA.libacc.BNNSGraphGetOutputNames_v2(g::LA.bnns_graph_t, func::Ptr{Cchar},
        n::Csize_t, names::Ptr{Ptr{Cchar}})::Cint
_graph_context_destroy_v2(c) =
    @ccall LA.libacc.BNNSGraphContextDestroy_v2(c::LA.bnns_graph_context_t)::Cvoid
_graph_context_workspace_size_v2(c, func) =
    @ccall LA.libacc.BNNSGraphContextGetWorkspaceSize_v2(c::LA.bnns_graph_context_t,
        func::Ptr{Cchar})::Csize_t
_graph_context_set_batch_size_v2(c, func, n) =
    @ccall LA.libacc.BNNSGraphContextSetBatchSize_v2(c::LA.bnns_graph_context_t,
        func::Ptr{Cchar}, n::UInt64)::Cint
_graph_context_set_dynamic_shapes_v2(c, func, n, shapes) =
    @ccall LA.libacc.BNNSGraphContextSetDynamicShapes_v2(c::LA.bnns_graph_context_t,
        func::Ptr{Cchar}, n::Csize_t, shapes::Ptr{_BNNSGraphShape})::Cint
_graph_context_execute_v2(c, func, n, args, wsize, w) =
    @ccall LA.libacc.BNNSGraphContextExecute_v2(c::LA.bnns_graph_context_t, func::Ptr{Cchar},
        n::Csize_t, args::Ptr{LA.bnns_graph_argument_t}, wsize::Csize_t, w::Ptr{Cchar})::Cint

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

Compile the compiled Core ML model (`.mlmodelc` directory, or the `model.mil`
inside it) at `filename` — optionally only the named `func` inside it — into an
executable graph via `BNNSGraphCompileFromFile`. Throws if BNNS cannot compile the
model (unsupported op, malformed program, missing file). The returned handle feeds
[`BNNSGraphContext`](@ref) and the graph-introspection helpers, and the compiled
graph's memory is released by a finalizer once the graph and every context made
from it are unreachable.

`.mlmodelc` is what Xcode / `xcrun coremlcompiler compile` produce from an
`.mlpackage`; only *ML Program* models (MIL), not the older NeuralNetwork format,
are accepted by BNNS.
"""
mutable struct BNNSGraph
    graph::LA.bnns_graph_t
    mapped::Bool                    # compiled into an mmap'd output file -> munmap, not free
    refs::Threads.Atomic{Int}       # 1 for the graph + 1 per live context
    function BNNSGraph(filename::AbstractString; func = nothing,
                       options::BNNSGraphCompileOptions = BNNSGraphCompileOptions())
        ispath(filename) || throw(ArgumentError("BNNSGraph: no such model: $filename"))
        fn = _cstr(filename)
        fckeep, fcptr = _fnarg(func)
        g = GC.@preserve fn fckeep options begin
            _graph_compile_v2(_cptr(fn), fcptr, options.opts)
        end
        g.data == C_NULL && error("BNNSGraphCompileFromFile failed to compile $(repr(filename))" *
                                  (func === nothing ? "" : " (function $(repr(func)))"))
        mapped = !isempty(bnns_compile_options_get_output_path(options)) ||
                 bnns_compile_options_get_output_fd(options) != -1
        obj = new(g, mapped, Threads.Atomic{Int}(1))
        finalizer(_graph_release, obj)
        return obj
    end
end

# Drop one reference; the compiled graph is released when the graph object and
# every context made from it are gone (finalizer order is unspecified, and a
# context reads the graph's memory until it is destroyed).
function _graph_release(g::BNNSGraph)
    if Threads.atomic_sub!(g.refs, 1) == 0 && g.graph.data != C_NULL
        if g.mapped
            @ccall munmap(g.graph.data::Ptr{Cvoid}, g.graph.size::Csize_t)::Cint
        else
            Libc.free(g.graph.data)
        end
        g.graph = LA.bnns_graph_t(C_NULL, 0)
    end
    return nothing
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
    GC.@preserve keep buf _bnns_check(_graph_input_names_v2(g.graph, fptr, Csize_t(n), pointer(buf)), "BNNSGraphGetInputNames")
    return _names_from(buf)
end
"Output names of `func` (`BNNSGraphGetOutputNames`)."
function bnns_graph_output_names(g::BNNSGraph, func = nothing)
    keep, fptr = _fnarg(func); n = bnns_graph_output_count(g, func)
    buf = fill(Ptr{Cchar}(C_NULL), n)
    GC.@preserve keep buf _bnns_check(_graph_output_names_v2(g.graph, fptr, Csize_t(n), pointer(buf)), "BNNSGraphGetOutputNames")
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
destroyed by a finalizer (`BNNSGraphContextDestroy`). It holds the mutable
execution state (dynamic shapes, streaming state), keeps its graph alive, and
must be used by **one thread at a time**; make one context per task for
concurrent inference. Feed it to [`bnns_graph_run`](@ref) /
[`bnns_graph_run!`](@ref), or to the low-level [`bnns_graph_execute!`](@ref).
"""
mutable struct BNNSGraphContext
    ctx::LA.bnns_graph_context_t
    graph::BNNSGraph
    function BNNSGraphContext(c::LA.bnns_graph_context_t, g::BNNSGraph)
        c.data == C_NULL && error("BNNSGraphContextMake returned NULL")
        Threads.atomic_add!(g.refs, 1)
        ctx = new(c, g)
        finalizer(ctx) do x
            _graph_context_destroy_v2(x.ctx)
            _graph_release(x.graph)
        end
        return ctx
    end
end

BNNSGraphContext(g::BNNSGraph) = BNNSGraphContext(LA.BNNSGraphContextMake(g.graph), g)

"""
    BNNSGraphContext(g, func, initial_states::Vector{BNNSTensor}) -> BNNSGraphContext

A streaming context (`BNNSGraphContextMakeStreaming`) seeded with recurrent
`initial_states`. Same finalizer as the plain constructor.
"""
function BNNSGraphContext(g::BNNSGraph, func, initial_states::Vector{LA.BNNSTensor})
    keep, fptr = _fnarg(func)
    c = GC.@preserve keep initial_states LA.BNNSGraphContextMakeStreaming(
        g.graph, fptr, Csize_t(length(initial_states)), pointer(initial_states))
    return BNNSGraphContext(c, g)
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

"""
    bnns_graph_context_workspace_size(c, func=nothing) -> Int

Workspace size (bytes) required to execute `func`
(`BNNSGraphContextGetWorkspaceSize`). Query it again after changing the batch size
or dynamic shapes. Allocate a suitable buffer with [`bnns_graph_workspace`](@ref).
"""
function bnns_graph_context_workspace_size(c::BNNSGraphContext, func = nothing)
    keep, fptr = _fnarg(func)
    sz = GC.@preserve keep _graph_context_workspace_size_v2(c.ctx, fptr)
    sz == typemax(Csize_t) && error("BNNSGraphContextGetWorkspaceSize failed")
    return Int(sz)
end

const _BNNS_PAGE = 16384   # ≥ the VM page size on both Apple Silicon (16K) and Intel (4K)

"""
    bnns_graph_workspace(c; func=nothing) -> Vector{UInt8}

A **page-aligned** scratch buffer of the size `func` currently needs
([`bnns_graph_context_workspace_size`](@ref)), for the `workspace` keyword of
[`bnns_graph_run!`](@ref) / [`bnns_graph_execute!`](@ref). `BNNSGraphContextExecute`
requires page alignment, which an ordinary `Vector{UInt8}` does not guarantee.
Reusing one buffer across calls makes execution allocation-free on the BNNS side;
without one BNNS allocates its own scratch on every call. The memory is released
when the vector is garbage collected.
"""
function bnns_graph_workspace(c::BNNSGraphContext; func = nothing)
    n = max(bnns_graph_context_workspace_size(c, func), 1)
    n = cld(n, _BNNS_PAGE) * _BNNS_PAGE
    p = Ref{Ptr{Cvoid}}(C_NULL)
    rc = @ccall posix_memalign(p::Ptr{Ptr{Cvoid}}, _BNNS_PAGE::Csize_t, n::Csize_t)::Cint
    rc == 0 || throw(OutOfMemoryError())
    return unsafe_wrap(Vector{UInt8}, Ptr{UInt8}(p[]), n; own = true)
end

"""
    bnns_graph_context_set_batch_size!(c, n; func=nothing) -> c

Set the batch size of a graph whose only dynamic dimension is a shared leading
(MIL-order) batch dimension (`BNNSGraphContextSetBatchSize`). For anything more
general use [`bnns_graph_context_set_dynamic_shapes!`](@ref).
"""
function bnns_graph_context_set_batch_size!(c::BNNSGraphContext, n::Integer; func = nothing)
    n >= 1 || throw(ArgumentError("BNNS: batch size must be ≥ 1"))
    keep, fptr = _fnarg(func)
    GC.@preserve keep _bnns_check(_graph_context_set_batch_size_v2(c.ctx, fptr, UInt64(n)),
                                  "BNNSGraphContextSetBatchSize")
    return c
end

"""
    bnns_graph_context_set_dynamic_shapes!(c, shapes; func=nothing) -> Vector{Pair{String,Dims}}

Bind concrete input shapes to a graph compiled with dynamic dimensions
(`BNNSGraphContextSetDynamicShapes`). `shapes` is a collection of
`name => dims` pairs for (some of) the graph's **inputs**; `dims` is given the
way Julia sees the array, i.e. `size(A)` of the array you will pass to
[`bnns_graph_run`](@ref) (the reverse of the MIL shape — see
[`bnns_graph_arguments`](@ref)). Inputs that are not mentioned keep the model's
default shape.

Returns the resulting shape of every argument, in execute order, again in Julia
order. A `0` in an *output* shape means BNNS cannot bound that dimension from the
input shapes alone (it depends on input values).
"""
function bnns_graph_context_set_dynamic_shapes!(c::BNNSGraphContext, shapes; func = nothing)
    names = bnns_graph_argument_names(c.graph, func)
    intents = bnns_graph_argument_intents(c.graph, func)
    bufs = Vector{UInt64}[UInt64[] for _ in names]
    for (k, dims) in _name_pairs(shapes)
        i = findfirst(==(k), names)
        i === nothing && throw(ArgumentError("BNNS: graph has no argument named $(repr(k)); arguments are $(names)"))
        intents[i] === :out && throw(ArgumentError("BNNS: $(repr(k)) is an output; only input shapes can be set"))
        bufs[i] = UInt64[reverse(collect(dims))...]
    end
    # Outputs get a rank-sized buffer so BNNS can report their deduced shape back.
    for i in eachindex(names)
        intents[i] === :out && (bufs[i] = zeros(UInt64, _graph_rank(c, names[i], func)))
    end
    keep, fptr = _fnarg(func)
    GC.@preserve keep bufs begin
        cs = [_BNNSGraphShape(Csize_t(length(b)), isempty(b) ? Ptr{UInt64}(C_NULL) : pointer(b)) for b in bufs]
        st = _graph_context_set_dynamic_shapes_v2(c.ctx, fptr, Csize_t(length(cs)), cs)
        st < 0 && error("BNNSGraphContextSetDynamicShapes failed with status $st")
    end
    info = bnns_graph_arguments(c; func)
    return [names[i] => (intents[i] === :out ? Tuple(Int.(reverse(bufs[i]))) : info[i].size)
            for i in eachindex(names)]
end

_graph_rank(c, name, func) = Int(bnns_graph_context_get_tensor(c, name; func, fill_shapes = false).rank)

# Normalise `name => value` collections (a pair, a tuple/vector of pairs, a Dict,
# a NamedTuple) to `String`-keyed pairs.
_name_pairs(x::Pair) = (String(first(x)) => last(x),)
_name_pairs(x::Union{AbstractDict,NamedTuple}) = [String(k) => v for (k, v) in pairs(x)]
_name_pairs(x::Union{Tuple,AbstractVector}) =
    all(p -> p isa Pair, x) ? [String(first(p)) => last(p) for p in x] :
        throw(ArgumentError("BNNS: expected `name => value` pairs"))

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

Low-level execute (`BNNSGraphContextExecute`): run `func` with raw argument
buffers, ordered as [`bnns_graph_argument_names`](@ref) reports them (outputs
first). The caller must keep the memory behind every argument alive for the call.
Prefer [`bnns_graph_run`](@ref) / [`bnns_graph_run!`](@ref), which build and
validate the arguments from Julia arrays.

`workspace` must be page-aligned — get one from [`bnns_graph_workspace`](@ref);
leave it empty to let BNNS allocate its own scratch.
"""
function bnns_graph_execute!(c::BNNSGraphContext, arguments::Vector{LA.bnns_graph_argument_t};
                             func = nothing, workspace::Vector{UInt8} = UInt8[])
    keep, fptr = _fnarg(func)
    n = bnns_graph_argument_count(c.graph, func)
    length(arguments) == n || throw(DimensionMismatch(
        "BNNS: graph function takes $n arguments, got $(length(arguments))"))
    if !isempty(workspace)
        UInt(pointer(workspace)) % _BNNS_PAGE == 0 || throw(ArgumentError(
            "BNNS: workspace must be page-aligned; allocate it with bnns_graph_workspace"))
        need = bnns_graph_context_workspace_size(c, func)
        length(workspace) >= need || throw(DimensionMismatch(
            "BNNS: workspace has $(length(workspace)) bytes, $need required"))
    end
    GC.@preserve keep arguments workspace begin
        wptr = isempty(workspace) ? Ptr{Cchar}(C_NULL) : Ptr{Cchar}(pointer(workspace))
        _bnns_check(_graph_context_execute_v2(c.ctx, fptr, Csize_t(length(arguments)),
            pointer(arguments), Csize_t(length(workspace)), wptr), "BNNSGraphContextExecute")
    end
    return c
end

# --- End-to-end inference with Julia arrays -------------------------------------

"""
    BNNSGraphArgument

Description of one argument of a graph function, as returned by
[`bnns_graph_arguments`](@ref):

  * `name::String`
  * `intent::Symbol` — `:in`, `:out` or `:inout`
  * `eltype::DataType` — Julia element type (`Float16`, `Float32`, `Int32`, `Bool`, …)
  * `shape::Dims` — the shape as written in the model (MIL / row-major order)
  * `size::Dims` — `size` of the Julia `Array` to pass for it: `reverse(shape)`

A dimension of `0` (in either tuple) is dynamic and not bound yet.
"""
struct BNNSGraphArgument
    name::String
    intent::Symbol
    eltype::DataType
    shape::Dims
    size::Dims
end

function Base.show(io::IO, a::BNNSGraphArgument)
    print(io, "BNNSGraphArgument(", repr(a.name), ", :", a.intent, ", ", a.eltype,
          ", shape=", a.shape, ", size=", a.size, ")")
end

"""
    bnns_graph_arguments(c::BNNSGraphContext; func=nothing) -> Vector{BNNSGraphArgument}
    bnns_graph_arguments(g::BNNSGraph; func=nothing)

Names, intents, element types and shapes of every argument of `func`, in execute
order (outputs first). Given a context, shapes reflect any batch size / dynamic
shapes already set on it.

## Memory layout: reverse the dimensions

Core ML / MIL tensors are **row-major**; Julia arrays are **column-major**. A MIL
tensor of shape `[N, C, H, W]` therefore has exactly the memory layout of a Julia
`Array` of size `(W, H, C, N)`. The graph functions here use that correspondence
— they pass Julia's memory to BNNS as is, with no copy — so every argument is a
Julia array whose `size` is the **reverse** of the model's shape
(`BNNSGraphArgument.size`). In particular a MIL matrix `[rows, cols]` is a Julia
`(cols, rows)` matrix, i.e. the transpose; use `permutedims` (or pass
`mil_order=true` to [`bnns_graph_run`](@ref)) when you want model index order.
BNNS graphs assume contiguous storage and ignore custom strides, so this is the
only zero-copy mapping.
"""
function bnns_graph_arguments(c::BNNSGraphContext; func = nothing)
    g = c.graph
    names = bnns_graph_argument_names(g, func)
    intents = bnns_graph_argument_intents(g, func)
    return map(names, intents) do name, intent
        t = bnns_graph_context_get_tensor(c, name; func, fill_shapes = true)
        shape = ntuple(i -> max(Int(t.shape[i]), 0), Int(t.rank))
        BNNSGraphArgument(name, intent, _julia_type(t.data_type), shape, reverse(shape))
    end
end

bnns_graph_arguments(g::BNNSGraph; func = nothing) =
    bnns_graph_arguments(BNNSGraphContext(g); func)

"""
    bnns_graph_run!(c, outputs, inputs; func=nothing, workspace=UInt8[]) -> outputs

Run `func` of the graph behind context `c`, reading `inputs` and writing into the
preallocated `outputs`. Both are collections of `name => Array` (a pair, a vector
or tuple of pairs, a `Dict`, or a `NamedTuple`); together they must supply every
argument of the function exactly once. Each array must be a dense `Array` whose
element type and `size` match [`bnns_graph_arguments`](@ref) — note the
**reversed-dimension layout** described there. Nothing is copied: BNNS reads and
writes the arrays' memory directly.

Pass a reusable page-aligned `workspace` from [`bnns_graph_workspace`](@ref) to
keep repeated calls allocation-free. A context must not be run from two threads at
once.
"""
function bnns_graph_run!(c::BNNSGraphContext, outputs, inputs; func = nothing,
                         workspace::Vector{UInt8} = UInt8[])
    info = bnns_graph_arguments(c; func)
    given = Dict{String,Array}()
    for (src, isout) in ((outputs, true), (inputs, false)), (k, A) in _name_pairs(src)
        i = findfirst(a -> a.name == k, info)
        i === nothing && throw(ArgumentError(
            "BNNS: graph has no argument named $(repr(k)); arguments are $([a.name for a in info])"))
        a = info[i]
        (isout ? a.intent !== :in : a.intent !== :out) || throw(ArgumentError(
            "BNNS: $(repr(k)) is an $(a.intent === :in ? "input" : "output"), passed as an $(isout ? "output" : "input")"))
        haskey(given, k) && throw(ArgumentError("BNNS: argument $(repr(k)) supplied more than once"))
        A isa Array || throw(ArgumentError(
            "BNNS: argument $(repr(k)) must be a dense Array, got $(typeof(A)); collect it first"))
        eltype(A) === a.eltype || throw(ArgumentError(
            "BNNS: argument $(repr(k)) must have element type $(a.eltype), got $(eltype(A))"))
        (ndims(A) == length(a.size) && all(d -> a.size[d] == 0 || a.size[d] == size(A, d), 1:ndims(A))) ||
            throw(DimensionMismatch("BNNS: argument $(repr(k)) must have size $(a.size) " *
                "(model shape $(a.shape), reversed for column-major storage), got $(size(A))"))
        given[k] = A
    end
    missing_args = [a.name for a in info if !haskey(given, a.name)]
    isempty(missing_args) || throw(ArgumentError("BNNS: missing graph arguments $(missing_args)"))
    arrays = Array[given[a.name] for a in info]
    GC.@preserve arrays begin
        args = [LA.bnns_graph_argument_t(Ptr{Cvoid}(pointer(A)), Csize_t(sizeof(A))) for A in arrays]
        _bnns_check(LA.BNNSGraphContextSetArgumentType(c.ctx,
            LA.BNNSGraphArgumentType(LA.BNNSGraphArgumentTypePointer)), "BNNSGraphContextSetArgumentType")
        bnns_graph_execute!(c, args; func, workspace)
    end
    return outputs
end

"""
    bnns_graph_run(c, inputs...; func=nothing, mil_order=false, workspace=UInt8[]) -> Dict{String,Array}

Run `func` on `inputs` (`name => Array` pairs, or a single `Dict`/`NamedTuple`),
allocating the outputs, and return them keyed by output name.

```julia
g = AppleAccelerate.BNNSGraph("classifier.mlmodelc")
c = AppleAccelerate.BNNSGraphContext(g)
AppleAccelerate.bnns_graph_arguments(c)          # names, eltypes, sizes
out = AppleAccelerate.bnns_graph_run(c, "image" => img)
probs = out["probabilities"]
```

By default arrays use the zero-copy **reversed-dimension** layout
([`bnns_graph_arguments`](@ref)): pass `size == BNNSGraphArgument.size`. With
`mil_order=true` inputs and outputs instead have the model's own shape
(`BNNSGraphArgument.shape`) and index order — `out[n, c, h, w]` means what it
means in the model — at the cost of one `permutedims` copy per array.

If the graph has dynamic dimensions, bind them first with
[`bnns_graph_context_set_batch_size!`](@ref) or
[`bnns_graph_context_set_dynamic_shapes!`](@ref); an output whose shape is still
unknown throws.
"""
function bnns_graph_run(c::BNNSGraphContext, inputs...; func = nothing, mil_order::Bool = false,
                        workspace::Vector{UInt8} = UInt8[])
    ins = Pair{String,Any}[]
    for x in inputs
        append!(ins, _name_pairs(x))
    end
    if mil_order
        ins = Pair{String,Any}[k => (A isa Array ? _reverse_dims(A) : A) for (k, A) in ins]
    end
    outs = Pair{String,Array}[]
    for a in bnns_graph_arguments(c; func)
        a.intent === :out || continue
        all(>(0), a.size) || error("BNNS: output $(repr(a.name)) has unresolved dynamic shape " *
            "$(a.shape); set the batch size / dynamic shapes on the context first")
        push!(outs, a.name => Array{a.eltype}(undef, a.size))
    end
    bnns_graph_run!(c, outs, ins; func, workspace)
    return Dict{String,Array}(k => (mil_order ? _reverse_dims(A) : A) for (k, A) in outs)
end

_reverse_dims(A::Array) = ndims(A) <= 1 ? A : permutedims(A, ndims(A):-1:1)

# =============================================================================
# Remaining DirectApply kernels
# =============================================================================

"""
    bnns_in_topk(input::Array, targets::Array{Int32}, K; dim=1) -> Array{Bool}

For each batch column, test whether the `targets` class index is among the top-`K`
scores of `input` along Julia dimension `dim` (`BNNSDirectApplyInTopK`). `input`
may be `Float32` or `Float16`.
"""
function bnns_in_topk(input::Array{<:Union{Float32,Float16}}, targets::Array{Int32}, K::Integer; dim::Integer = 1)
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
