## fftplan.jl — a common plan interface on top of the vDSP FFT / DFT / DCT setups ##
#
# `dsp.jl` exposes vDSP's setup objects directly (`FFTSetup`, `DFTSetup`,
# `InterleavedDFTSetup`) and threads them through `fft(x, setup)`-style calls. A setup
# knows neither its transform direction nor its shape, so it cannot be *applied*. An
# `FFTPlan` pairs a setup with the transform kind, direction, shape and output scaling,
# which is what makes `plan * x`, `mul!(y, plan, x)`, `inv(plan)` and `plan \ x`
# well-defined.
#
# Everything here is a method on an AppleAccelerate-owned type. Nothing is exported and
# no other package's function is extended for foreign types, so loading AppleAccelerate
# still changes nothing about `FFTW.fft`/`AbstractFFTs.plan_fft` (see
# docs/src/extensions.md and issue #139).
#
# ## Execution is zero-copy where vDSP allows it
#
# vDSP's FFTs take *split*-complex operands (separate real/imag arrays) with a stride,
# while Julia stores `Complex{T}` interleaved. An interleaved buffer *is* a
# split-complex operand with stride 2: `realp = p`, `imagp = p + sizeof(T)`. The plans
# use that to transform straight out of `x` into `y` with no packing buffers, so `mul!`
# does not allocate and a plan holds no mutable scratch. Backends, in preference order:
#
#   * `vDSP_DFT_Interleaved_Execute` — native interleaved, fastest, but only lengths
#     ≤ 4096 (≤ 1024 for the mixed-radix f·2^k lengths); the create call returns NULL
#     otherwise and the next backend is used.
#   * `vDSP_fft_zop` / `fft2d_zop` / `fftm_zop` / `fft_zrop` with stride 2 — any power
#     of two, 1-D, 2-D and batched. Stride 2 misses vDSP's unit-stride kernels, so
#     above ~4096 points this is slower than `fft(x, setup)`'s copy-in/copy-out
#     (≈105 µs vs ≈83 µs at 16384, M-series) — the price of a scratch-free plan.
#   * `vDSP_DFT_Execute` — large mixed-radix lengths. This call has no stride argument,
#     so these plans pack into temporaries and `mul!` allocates. The same holds for the
#     2-D real plans, whose packed DC/Nyquist row cannot be unpacked in place.
#
# ## Thread safety
#
# vDSP setups are read-only during execution (vDSP.h, "Multithreading"), plans are
# immutable and scratch-free, and setups are taken from the shared caches under
# `_SETUP_CACHE_LOCK`. One plan can therefore be applied from many tasks at once.

const _PlanReal = Union{Float32,Float64}

"""
    FFTPlan{T,K,N,S}

A reusable transform plan: a vDSP setup together with the transform kind `K`
(`:fft`, `:fftm`, `:rfft`, `:brfft` or `:dct`), the direction, the input/output shape
(`N` dimensions) and the output scaling. `T` is the real precision and `S` the vDSP
setup type backing the plan.

Create one with [`fftplan`](@ref), [`bfftplan`](@ref), [`ifftplan`](@ref),
[`rfftplan`](@ref), [`brfftplan`](@ref), [`irfftplan`](@ref) or [`dctplan`](@ref), then
apply it with `plan * x` or `mul!(y, plan, x)`, and invert it with `inv(plan)` or
`plan \\ y`. `size(plan)` is the input size, [`output_size`](@ref) the output size, and
`eltype(plan)` the input element type.

Plans are immutable and hold no scratch memory, so a single plan may be applied
concurrently from several tasks. The underlying vDSP setup is freed by its own
finalizer once no plan references it.
"""
struct FFTPlan{T<:_PlanReal,K,N,S}
    setup::S
    insize::NTuple{N,Int}
    outsize::NTuple{N,Int}
    direction::Int   # FFT_FORWARD / FFT_INVERSE
    scale::T         # the transform output is multiplied by this (1 ⇒ skipped)
    dims::Int        # :fftm — the transformed dimension; otherwise 0
    dcttype::Int     # :dct — 2, 3 or 4; otherwise 0

    function FFTPlan{T,K}(setup::S, insize::NTuple{N,Int}, outsize::NTuple{N,Int},
                          direction::Integer, scale::Real;
                          dims::Integer=0, dcttype::Integer=0) where {T<:_PlanReal,K,N,S}
        new{T,K,N,S}(setup, insize, outsize, Int(direction), T(scale), Int(dims), Int(dcttype))
    end
end

# --- Introspection ---

_plankind(::FFTPlan{T,K}) where {T,K} = K

Base.size(p::FFTPlan) = p.insize
Base.size(p::FFTPlan, d::Integer) = d <= length(p.insize) ? p.insize[d] : 1
Base.ndims(::FFTPlan{T,K,N}) where {T,K,N} = N
Base.length(p::FFTPlan) = prod(p.insize)

Base.eltype(::Type{<:FFTPlan{T,K}}) where {T,K} = K in (:rfft, :dct) ? T : Complex{T}
Base.eltype(p::FFTPlan) = eltype(typeof(p))

_output_eltype(::FFTPlan{T,K}) where {T,K} = K in (:brfft, :dct) ? T : Complex{T}

"""
    output_size(plan::FFTPlan) -> Dims

Size of the array produced by `plan * x` (and required of `y` in `mul!(y, plan, x)`).
`size(plan)` is the matching *input* size; the two differ only for the real transforms.
"""
output_size(p::FFTPlan) = p.outsize

# Number of samples the transform sums over; the forward/backward round trip scales by this.
function _transform_length(p::FFTPlan)
    K = _plankind(p)
    K === :fftm && return p.insize[p.dims]
    K === :brfft && return prod(p.outsize)
    return prod(p.insize)
end

_backend_name(::InterleavedDFTSetup) = "vDSP_DFT_Interleaved"
_backend_name(::DFTSetup) = "vDSP_DFT"
_backend_name(::FFTSetup) = "vDSP_fft"

function _transform_name(p::FFTPlan)
    K = _plankind(p)
    K === :dct && return "DCT-" * ("II", "III", "IV")[p.dcttype - 1]
    fwd = p.direction == FFT_FORWARD
    K === :fft && return (fwd ? "forward" : "backward") * (ndims(p) == 2 ? " 2-D FFT" : " FFT")
    K === :fftm && return (fwd ? "forward" : "backward") * " FFT along dims=$(p.dims)"
    return (K === :rfft ? "forward" : "backward") * (ndims(p) == 2 ? " 2-D real FFT" : " real FFT")
end

function Base.show(io::IO, p::FFTPlan{T}) where {T}
    print(io, "AppleAccelerate.FFTPlan{", T, "}: ", _transform_name(p), ", ",
          join(p.insize, "×"), " ", eltype(p), " → ",
          join(p.outsize, "×"), " ", _output_eltype(p))
    isone(p.scale) || print(io, ", scaled by ", p.scale)
    print(io, " (", _backend_name(p.setup), ")")
end

# --- Setup acquisition ---

# Interleaved-DFT setups are fixed-length and direction-specific; cache them like the
# others. `nothing` records a length the create call rejected, so it is probed only once.
const _IDFT_SETUP_CACHE = Dict{Tuple{DataType,Int,Int},Union{Nothing,InterleavedDFTSetup}}()

for (T, createfn) in ((Float32, :vDSP_DFT_Interleaved_CreateSetup),
                      (Float64, :vDSP_DFT_Interleaved_CreateSetupD))
    @eval function _cached_idftsetup(::Type{$T}, n::Int, direction::Int)
        return lock(_SETUP_CACHE_LOCK) do
            get!(_IDFT_SETUP_CACHE, ($T, n, direction)) do
                ptr = LibAccelerate.$createfn(C_NULL, n, Cint(direction), Cint(0))
                ptr == C_NULL ? nothing : InterleavedDFTSetup{$T}(Ptr{Cvoid}(ptr), n, direction)
            end
        end::Union{Nothing,InterleavedDFTSetup{$T}}
    end
end

@noinline _plan_pow2_error(what, sz) = throw(ArgumentError(string(
    what, " requires every transformed dimension to be a power of two; got size ", sz,
    ". vDSP has no 2-D or batched mixed-radix FFT; use FFTW.jl for other sizes.")))

# --- Constructors: complex transforms ---

function _complex_plan(::Type{T}, sz::NTuple{N,Int}, direction::Int, normalized::Bool,
                       dims) where {T<:_PlanReal,N}
    if N == 1
        dims === nothing || dims == 1 || throw(ArgumentError("dims must be 1 for a vector; got $dims"))
        n = sz[1]
        is_supported_fft_length(n) || _unsupported_fft_length(n)
        scale = normalized ? inv(T(n)) : one(T)
        setup = _cached_idftsetup(T, n, direction)
        if setup === nothing
            setup = ispow2(n) ? _cached_fftsetup(T, n) : _cached_dftsetup(T, n, direction)
        end
        return FFTPlan{T,:fft}(setup, sz, sz, direction, scale)
    elseif N == 2
        if dims === nothing
            all(ispow2, sz) || _plan_pow2_error("a 2-D FFT plan", sz)
            scale = normalized ? inv(T(prod(sz))) : one(T)
            return FFTPlan{T,:fft}(_cached_fftsetup(T, max(sz...)), sz, sz, direction, scale)
        end
        dims == 1 || dims == 2 || throw(ArgumentError("dims must be 1 or 2; got $dims"))
        n = sz[dims]
        ispow2(n) || _plan_pow2_error("a batched FFT plan", sz)
        scale = normalized ? inv(T(n)) : one(T)
        return FFTPlan{T,:fftm}(_cached_fftsetup(T, n), sz, sz, direction, scale; dims=dims)
    end
    throw(ArgumentError("vDSP FFT plans support 1-D and 2-D arrays only; got $N dimensions"))
end

_plan_dims(sz::Integer) = (Int(sz),)
_plan_dims(sz::Tuple{Vararg{Integer}}) = map(Int, sz)

"""
    fftplan(x::StridedVecOrMat{Complex{T}}, [dims])
    fftplan(Complex{T}, size, [dims])

Plan a forward complex FFT for arrays shaped like `x` (only the element type and size of
`x` are used; its contents are untouched). Returns an [`FFTPlan`](@ref).

- Vector: any length accepted by [`is_supported_fft_length`](@ref) — a power of two, or
  `f·2^k` with `f ∈ {3, 5, 15}`, `k ≥ 3`.
- Matrix: a full 2-D transform; both dimensions must be powers of two.
- Matrix with `dims = 1` or `2`: independent 1-D transforms of every column or row
  (power-of-two length), like FFTW's `plan_fft(x, dims)`.

```jldoctest
julia> using LinearAlgebra: mul!

julia> x = ComplexF64[1, 2, 3, 4];

julia> p = AppleAccelerate.fftplan(x)
AppleAccelerate.FFTPlan{Float64}: forward FFT, 4 ComplexF64 → 4 ComplexF64 (vDSP_fft)

julia> y = p * x                  # allocate the result
4-element Vector{ComplexF64}:
 10.0 + 0.0im
 -2.0 + 2.0im
 -2.0 + 0.0im
 -2.0 - 2.0im

julia> mul!(similar(x), p, x) == y   # no allocation; `mul!(x, p, x)` transforms in place
true

julia> p \\ y ≈ x                  # same as inv(p) * y
true
```

See also [`bfftplan`](@ref), [`ifftplan`](@ref), [`rfftplan`](@ref).
"""
fftplan(x::StridedVecOrMat{Complex{T}}, dims=nothing) where {T<:_PlanReal} =
    _complex_plan(T, size(x), FFT_FORWARD, false, dims)
fftplan(::Type{Complex{T}}, sz, dims=nothing) where {T<:_PlanReal} =
    _complex_plan(T, _plan_dims(sz), FFT_FORWARD, false, dims)

"""
    bfftplan(x::StridedVecOrMat{Complex{T}}, [dims])
    bfftplan(Complex{T}, size, [dims])

Plan an unnormalized backward (inverse) complex FFT: `bfftplan(x) * (fftplan(x) * x)`
is `x` times the transform length. See [`fftplan`](@ref) for the supported shapes and
[`ifftplan`](@ref) for the normalized inverse.
"""
bfftplan(x::StridedVecOrMat{Complex{T}}, dims=nothing) where {T<:_PlanReal} =
    _complex_plan(T, size(x), FFT_INVERSE, false, dims)
bfftplan(::Type{Complex{T}}, sz, dims=nothing) where {T<:_PlanReal} =
    _complex_plan(T, _plan_dims(sz), FFT_INVERSE, false, dims)

"""
    ifftplan(x::StridedVecOrMat{Complex{T}}, [dims])
    ifftplan(Complex{T}, size, [dims])

Plan a normalized inverse complex FFT, so that `ifftplan(x) * (fftplan(x) * x) ≈ x`.
This is the plan `inv(fftplan(x))` returns. See [`fftplan`](@ref) for the supported
shapes.
"""
ifftplan(x::StridedVecOrMat{Complex{T}}, dims=nothing) where {T<:_PlanReal} =
    _complex_plan(T, size(x), FFT_INVERSE, true, dims)
ifftplan(::Type{Complex{T}}, sz, dims=nothing) where {T<:_PlanReal} =
    _complex_plan(T, _plan_dims(sz), FFT_INVERSE, true, dims)

# --- Constructors: real transforms ---

# `rsz` is the size of the *real* array; the spectrum has size (rsz[1]÷2+1, rsz[2:end]...).
function _real_plan(::Type{T}, rsz::NTuple{N,Int}, forward::Bool, scale::Real) where {T<:_PlanReal,N}
    N == 1 || N == 2 || throw(ArgumentError(
        "vDSP FFT plans support 1-D and 2-D arrays only; got $N dimensions"))
    n = rsz[1]
    if N == 2
        (all(ispow2, rsz) && all(>=(2), rsz)) || _plan_pow2_error("a 2-D real FFT plan (each dimension ≥ 2)", rsz)
        setup = _cached_fftsetup(T, max(rsz...))
    elseif ispow2(n)
        n >= 2 || throw(ArgumentError("a real FFT plan needs length ≥ 2; got $n"))
        setup = _cached_fftsetup(T, n)
    else
        (iseven(n) && is_supported_fft_length(n)) || _unsupported_rfft_length(n)
        # Throws the same ArgumentError if vDSP rejects the length (it needs k ≥ 4).
        setup = _cached_rdftsetup(T, n, forward ? DFT_FORWARD : DFT_INVERSE)
    end
    csz = (n >> 1 + 1, Base.tail(rsz)...)
    return forward ? FFTPlan{T,:rfft}(setup, rsz, csz, FFT_FORWARD, scale) :
                     FFTPlan{T,:brfft}(setup, csz, rsz, FFT_INVERSE, scale)
end

function _real_size(X::StridedVecOrMat, n::Integer)
    rsz = (Int(n), Base.tail(size(X))...)
    size(X, 1) == n >> 1 + 1 || throw(DimensionMismatch(
        "a real signal with first dimension $n has a spectrum with first dimension " *
        "$(n >> 1 + 1); got $(size(X, 1))"))
    return rsz
end

"""
    rfftplan(x::StridedVecOrMat{T})
    rfftplan(T, size)

Plan a forward FFT of a real array, producing the non-redundant half spectrum in the
same layout as FFTW's `rfft`: length `n÷2+1` for a vector of length `n`, size
`(n1÷2+1)×n2` for an `n1×n2` matrix.

Vector lengths may be a power of two (≥ 2) or a mixed-radix `f·2^k` length supported by
Apple's real-input DFT (`f ∈ {3, 5, 15}`, `k ≥ 4`); matrix dimensions must be powers of
two. `mul!` is allocation-free for power-of-two vectors; the mixed-radix and 2-D paths
pack through temporaries.

`inv(plan)` / `plan \\ X` give the normalized inverse ([`irfftplan`](@ref)).
"""
rfftplan(x::StridedVecOrMat{T}) where {T<:_PlanReal} = _real_plan(T, size(x), true, 1)
rfftplan(::Type{T}, sz) where {T<:_PlanReal} = _real_plan(T, _plan_dims(sz), true, 1)

"""
    brfftplan(X::StridedVecOrMat{Complex{T}}, n::Integer)

Plan the unnormalized inverse of [`rfftplan`](@ref): maps a half spectrum shaped like
`X` back to a real array whose first dimension is `n` (`size(X, 1)` must be `n÷2+1`).
The result is the original signal times the number of samples; use
[`irfftplan`](@ref) for the normalized inverse.
"""
brfftplan(X::StridedVecOrMat{Complex{T}}, n::Integer) where {T<:_PlanReal} =
    _real_plan(T, _real_size(X, n), false, 1)

"""
    irfftplan(X::StridedVecOrMat{Complex{T}}, n::Integer)

Plan the normalized inverse real FFT, so that
`irfftplan(X, size(x, 1)) * (rfftplan(x) * x) ≈ x`. See [`brfftplan`](@ref).
"""
function irfftplan(X::StridedVecOrMat{Complex{T}}, n::Integer) where {T<:_PlanReal}
    rsz = _real_size(X, n)
    return _real_plan(T, rsz, false, inv(T(prod(rsz))))
end

# --- Constructors: DCT ---

"""
    dctplan(x::StridedVector{Float32}, [dct_type = 2])
    dctplan(Float32, n, [dct_type = 2])

Plan an (unnormalized) discrete cosine transform of type II, III or IV
(`dct_type = 2, 3, 4`), matching [`dct`](@ref). The length must be `f·2^k` with
`f ∈ {1, 3, 5, 15}` and `k ≥ 4`. vDSP provides the DCT in single precision only.

`inv(plan)` uses the DCT-II/DCT-III duality (and DCT-IV's self-inverse property), scaled
by `2/n`, so `plan \\ (plan * x) ≈ x`.
"""
dctplan(x::StridedVector{Float32}, dct_type::Integer=2) = _dct_plan(length(x), Int(dct_type), 1f0)
dctplan(::Type{Float32}, n::Integer, dct_type::Integer=2) = _dct_plan(Int(n), Int(dct_type), 1f0)
dctplan(::StridedVector{Float64}, ::Integer=2) = _dct_no_float64()
dctplan(::Type{Float64}, ::Integer, ::Integer=2) = _dct_no_float64()

function _dct_plan(n::Int, dct_type::Int, scale::Float32)
    2 <= dct_type <= 4 || throw(ArgumentError(
        "DCT type $dct_type is not supported; vDSP provides types 2, 3 and 4"))
    (trailing_zeros(n) >= 4 && _odd_cofactor(n) in (1, 3, 5, 15)) || throw(ArgumentError(
        "unsupported DCT length $n: vDSP requires f*2^k with f ∈ {1, 3, 5, 15} and k ≥ 4"))
    return FFTPlan{Float32,:dct}(plan_dct(n, dct_type), (n,), (n,), FFT_FORWARD, scale;
                                 dcttype=dct_type)
end

# --- Inversion ---

"""
    inv(plan::FFTPlan)
    plan \\ y

`inv(plan)` is the plan of the inverse transform, normalized so that
`inv(plan) * (plan * x) ≈ x`; `plan \\ y` is `inv(plan) * y`. The inverse of a forward
plan is the matching `ifftplan`/`irfftplan`; the inverse of an unnormalized backward
plan is a forward plan scaled by `1/n`. Constructing the inverse only looks up a cached
vDSP setup, so it is cheap; hold on to it if you apply it in a hot loop.
"""
function Base.inv(p::FFTPlan{T,K}) where {T,K}
    scale = inv(_transform_length(p) * p.scale)
    if K === :fft || K === :fftm
        q = _complex_plan(T, p.insize, -p.direction, false, K === :fftm ? p.dims : nothing)
        return FFTPlan{T,K}(q.setup, q.insize, q.outsize, q.direction, scale; dims=q.dims)
    elseif K === :rfft
        return _real_plan(T, p.insize, false, scale)
    elseif K === :brfft
        return _real_plan(T, p.outsize, true, scale)
    else # :dct — II ↔ III, IV ↔ IV; vDSP's round trip scales by n/2
        scale isa Float32 || _dct_no_float64()
        return _dct_plan(p.insize[1], (3, 2, 4)[p.dcttype - 1], 2 * scale)
    end
end

Base.:\(p::FFTPlan, y::StridedArray) = inv(p) * y

# --- Application ---

function Base.:*(p::FFTPlan, x::StridedArray)
    y = Array{_output_eltype(p)}(undef, p.outsize)
    return LinearAlgebra.mul!(y, p, x)
end

"""
    mul!(y, plan::FFTPlan, x)

Apply `plan` to `x`, writing the result into `y` and returning `y`. `x` and `y` must be
contiguous (unit-stride) arrays of size `size(plan)` and [`output_size`](@ref)`(plan)`.
For the complex plans `y` may be `x` itself, which transforms in place; partially
overlapping `x` and `y` are not supported.

`mul!` does not allocate, except for the large mixed-radix, mixed-radix real and 2-D real
plans (see [`FFTPlan`](@ref)).
"""
function LinearAlgebra.mul!(y::StridedArray, p::FFTPlan, x::StridedArray)
    eltype(x) === eltype(p) || throw(ArgumentError(
        "plan expects input elements of type $(eltype(p)); got $(eltype(x))"))
    eltype(y) === _output_eltype(p) || throw(ArgumentError(
        "plan produces output elements of type $(_output_eltype(p)); got $(eltype(y))"))
    size(x) == p.insize || throw(DimensionMismatch(
        "plan expects input of size $(p.insize); got $(size(x))"))
    size(y) == p.outsize || throw(DimensionMismatch(
        "plan produces output of size $(p.outsize); got $(size(y))"))
    _check_contiguous(x)   # separately: a heterogeneous vararg tuple would allocate
    _check_contiguous(y)
    _execute!(y, p, x)
    isone(p.scale) || LinearAlgebra.rmul!(y, p.scale)
    return y
end

# Large mixed-radix complex: vDSP_DFT_Execute has no strides, so pack via `dft`.
function _execute!(y, p::FFTPlan{T,:fft,1,DFTSetup{T}}, x) where {T}
    copyto!(y, dft(x isa Vector ? x : Vector(x), p.setup))
end

# DCT (Float32 only).
function _execute!(y, p::FFTPlan{Float32,:dct,1}, x)
    setup = p.setup
    GC.@preserve x y setup begin
        LibAccelerate.vDSP_DCT_Execute(setup.setup, pointer(x), pointer(y))
    end
    return y
end

# Mixed-radix real: pack/unpack through the existing DFT helpers.
_execute!(y, ::FFTPlan{T,:rfft,1,DFTSetup{T}}, x) where {T} =
    copyto!(y, _rfft1d_dft(x isa Vector ? x : Vector(x)))
_execute!(y, p::FFTPlan{T,:brfft,1,DFTSetup{T}}, x) where {T} =
    copyto!(y, _brfft1d_dft(x isa Vector ? x : Vector(x), p.outsize[1]))

# 2-D real: the packed k0 = 0 / k0 = n1÷2 row cannot be unpacked in place.
_execute!(y, p::FFTPlan{T,:rfft,2}, x) where {T} =
    copyto!(y, _rfft2d(x isa Matrix ? x : Matrix(x), p.setup))
_execute!(y, p::FFTPlan{T,:brfft,2}, x) where {T} =
    copyto!(y, _brfft2d(x isa Matrix ? x : Matrix(x), p.outsize[1], p.setup))

# The out-of-place vDSP calls take two `const DSPSplitComplex *` operands. Two separate
# `Ref`s are not reliably stack-allocated (Julia 1.11 heap-allocates both), so keep the
# pair in one `Ref` and hand vDSP pointers into it; callers `GC.@preserve` the `Ref`.
@inline function _pair_ptrs(ops::Base.RefValue{NTuple{2,SC}}) where {SC}
    pa = Ptr{SC}(Base.unsafe_convert(Ptr{NTuple{2,SC}}, ops))
    return pa, pa + sizeof(SC)
end

for (T, SC, CT, idft_exec, zop, zip, zop2d, zip2d, zopm, zipm, zrop, zrip) in (
        (Float64, :DSPDoubleSplitComplex, :DSPDoubleComplex, :vDSP_DFT_Interleaved_ExecuteD,
         :vDSP_fft_zopD, :vDSP_fft_zipD, :vDSP_fft2d_zopD, :vDSP_fft2d_zipD,
         :vDSP_fftm_zopD, :vDSP_fftm_zipD, :vDSP_fft_zropD, :vDSP_fft_zripD),
        (Float32, :DSPSplitComplex, :DSPComplex, :vDSP_DFT_Interleaved_Execute,
         :vDSP_fft_zop, :vDSP_fft_zip, :vDSP_fft2d_zop, :vDSP_fft2d_zip,
         :vDSP_fftm_zop, :vDSP_fftm_zip, :vDSP_fft_zrop, :vDSP_fft_zrip))
    @eval begin
        # View an interleaved complex (or even/odd real) buffer as a stride-2 split-complex.
        @inline _split2(p::Ptr{<:Union{$T,Complex{$T}}}) = $SC(Ptr{$T}(p), Ptr{$T}(p) + sizeof($T))

        # 1-D complex, native interleaved DFT. `Ori` may equal `Iri`.
        function _execute!(y, p::FFTPlan{$T,:fft,1,InterleavedDFTSetup{$T}}, x)
            setup = p.setup
            GC.@preserve x y setup begin
                LibAccelerate.$idft_exec(setup.setup, Ptr{LibAccelerate.$CT}(pointer(x)),
                                         Ptr{LibAccelerate.$CT}(pointer(y)))
            end
            return y
        end

        # 1-D complex, power of two.
        function _execute!(y, p::FFTPlan{$T,:fft,1,FFTSetup{$T}}, x)
            setup = p.setup
            logn = trailing_zeros(p.insize[1])
            GC.@preserve x y setup begin
                if pointer(x) == pointer(y)
                    LibAccelerate.$zip(setup.plan, Ref(_split2(pointer(y))), 2, logn, p.direction)
                else
                    ops = Ref((_split2(pointer(x)), _split2(pointer(y))))
                    GC.@preserve ops begin
                        pa, pc = _pair_ptrs(ops)
                        LibAccelerate.$zop(setup.plan, pa, 2, pc, 2, logn, p.direction)
                    end
                end
            end
            return y
        end

        # 2-D complex. Julia's first (contiguous) dimension is vDSP's N0.
        function _execute!(y, p::FFTPlan{$T,:fft,2,FFTSetup{$T}}, x)
            setup = p.setup
            n1, n2 = p.insize
            log2n1 = trailing_zeros(n1); log2n2 = trailing_zeros(n2)
            GC.@preserve x y setup begin
                if pointer(x) == pointer(y)
                    LibAccelerate.$zip2d(setup.plan, Ref(_split2(pointer(y))), 2, 2n1,
                                         log2n1, log2n2, p.direction)
                else
                    ops = Ref((_split2(pointer(x)), _split2(pointer(y))))
                    GC.@preserve ops begin
                        pa, pc = _pair_ptrs(ops)
                        LibAccelerate.$zop2d(setup.plan, pa, 2, 2n1, pc, 2, 2n1,
                                             log2n1, log2n2, p.direction)
                    end
                end
            end
            return y
        end

        # Batched 1-D complex along `p.dims`.
        function _execute!(y, p::FFTPlan{$T,:fftm,2,FFTSetup{$T}}, x)
            setup = p.setup
            n1, n2 = p.insize
            logn = trailing_zeros(p.insize[p.dims])
            elstride, batchstride, m = p.dims == 1 ? (2, 2n1, n2) : (2n1, 2, n1)
            GC.@preserve x y setup begin
                if pointer(x) == pointer(y)
                    LibAccelerate.$zipm(setup.plan, Ref(_split2(pointer(y))), elstride, batchstride,
                                        logn, m, p.direction)
                else
                    ops = Ref((_split2(pointer(x)), _split2(pointer(y))))
                    GC.@preserve ops begin
                        pa, pc = _pair_ptrs(ops)
                        LibAccelerate.$zopm(setup.plan, pa, elstride, batchstride,
                                            pc, elstride, batchstride, logn, m, p.direction)
                    end
                end
            end
            return y
        end

        # 1-D real forward, power of two. The real input is read as even/odd
        # split-complex (stride 2) and the packed spectrum lands in y[1:half]; the
        # DC/Nyquist pair vDSP packs into element 1 is then split, and the forward
        # scaling by 2 removed.
        function _execute!(y, p::FFTPlan{$T,:rfft,1,FFTSetup{$T}}, x)
            setup = p.setup
            n = p.insize[1]
            half = n >> 1
            GC.@preserve x y setup begin
                ops = Ref((_split2(pointer(x)), _split2(pointer(y))))
                GC.@preserve ops begin
                    pa, pc = _pair_ptrs(ops)
                    LibAccelerate.$zrop(setup.plan, pa, 2, pc, 2, trailing_zeros(n), FFT_FORWARD)
                end
            end
            @inbounds begin
                packed = y[1]
                y[half + 1] = complex(imag(packed) / 2)
                y[1] = complex(real(packed) / 2)
                @simd for k in 2:half
                    y[k] /= 2
                end
            end
            return y
        end

        # 1-D real backward, power of two. The output buffer doubles as the packed
        # input (n reals = n÷2 interleaved pairs), so the spectrum `x` is left intact.
        function _execute!(y, p::FFTPlan{$T,:brfft,1,FFTSetup{$T}}, x)
            setup = p.setup
            n = p.outsize[1]
            half = n >> 1
            @inbounds begin
                y[1] = real(x[1])
                y[2] = real(x[half + 1])
                for k in 2:half
                    y[2k - 1] = real(x[k])
                    y[2k] = imag(x[k])
                end
            end
            GC.@preserve y setup begin
                LibAccelerate.$zrip(setup.plan, Ref(_split2(pointer(y))), 2,
                                    trailing_zeros(n), FFT_INVERSE)
            end
            return y
        end
    end
end

# --- Explicit escape hatch for the one-shot API ---
#
# `fft(x; fallback = FFTW.fft)` and friends: when vDSP cannot transform `x`, call the
# user-supplied function instead of throwing. Opt-in per call site — AppleAccelerate
# never picks another backend on its own (issue #139).
# Whether vDSP's real-input mixed-radix DFT has a setup for length `n`. The accepted
# set is irregular (it differs between Float32 and Float64), so ask the create call.
function _rdft_supported(::Type{T}, n::Int) where {T<:_PlanReal}
    try
        _cached_rdftsetup(T, n, DFT_FORWARD)
        return true
    catch err
        err isa ArgumentError || rethrow()
        return false
    end
end

_fft_fallback(::Nothing, err, n, _args...) = err(n)
_fft_fallback(fallback, _err, _n, args...) = fallback(args...)
