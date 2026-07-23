# BNNS (Basic Neural Network Subroutines) idiomatic wrappers.
#
# This module wraps the Apple-CURRENT (non-deprecated) slice of BNNS exposed by
# the raw `LibAccelerate` layer: dense array descriptors (`BNNSArray`), the
# stateless tensor ops that remain current (`bnns_transpose`, `bnns_copy!`), the
# DirectApply reduction / top-k kernels (`bnns_reduce`, `bnns_topk`,
# `bnns_in_topk`), random generation, nearest neighbors, layout/size queries, and
# the full BNNS Graph compiler pipeline. The numerically important helpers are
# cross-validated against plain-Julia references.
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
    BNNSDataTypeFloat32, BNNSDataTypeInt32,
    BNNSDataLayoutVector, BNNSDataLayoutColumnMajorMatrix

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

const _BNNS_REDUCE = Dict(
    :max => LA.BNNSReduceFunctionMax, :min => LA.BNNSReduceFunctionMin,
    :sum => LA.BNNSReduceFunctionSum, :mean => LA.BNNSReduceFunctionMean,
    :sumsquare => LA.BNNSReduceFunctionSumSquare, :l1 => LA.BNNSReduceFunctionL1Norm,
    :l2 => LA.BNNSReduceFunctionL2Norm, :product => LA.BNNSReduceFunctionProduct,
    :logsumexp => LA.BNNSReduceFunctionLogSumExp,
)

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
