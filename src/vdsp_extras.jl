## vdsp_extras.jl — closing the last vDSP coverage gaps ##
#
# An exact audit of the raw layer (`gen/coverage_audit.jl`, which walks the lowered IR
# of every method rather than grepping source text) shows the idiomatic layer reaching
# 459 of the 469 generated `vDSP_*` functions. This file wraps the remaining ones that
# have a natural Julia surface, adds a type-dispatched front end over the existing
# `vfix*` / `vflt*` families, and records the final tally in `VDSP_COVERAGE`.

# =====================================================================
# int ↔ float conversion: type-dispatched front end over vfix* / vflt*
# =====================================================================
#
# `array.jl` wraps all 48 `vDSP_vfix*` / `vDSP_vflt*` kernels under their C names
# (`vfixr16`, `vfltu8`, …), which encode signedness, rounding and width in the name.
# `float_to_int` / `int_to_float` select the same kernels from the Julia *types*.

const _VDSPInt = Union{Int8,Int16,Int32,UInt8,UInt16,UInt32}

for (intT, width, u) in ((Int8, 8, ""), (Int16, 16, ""), (Int32, 32, ""),
                         (UInt8, 8, "u"), (UInt16, 16, "u"), (UInt32, 32, "u"))
    trunc! = Symbol("vfix", u, width, "!")
    round! = Symbol("vfixr", u, width, "!")
    flt!   = Symbol("vflt", u, width, "!")
    @eval begin
        _float_to_int!(C::StridedVector{$intT}, A, ::Val{:trunc})   = $trunc!(C, A)
        _float_to_int!(C::StridedVector{$intT}, A, ::Val{:nearest}) = $round!(C, A)
        _int_to_float!(C, A::StridedVector{$intT}) = $flt!(C, A)
    end
end

@noinline _bad_rounding(r) =
    throw(ArgumentError("rounding must be :trunc or :nearest; got $(repr(r))"))

"""
    float_to_int!(C::StridedVector{I}, A::StridedVector{T}; rounding = :trunc) -> C
    float_to_int(I, A::StridedVector{T}; rounding = :trunc) -> Vector{I}

Convert the floating-point vector `A` (`T` is `Float32` or `Float64`) to the integer
type `I ∈ (Int8, Int16, Int32, UInt8, UInt16, UInt32)`. Both arguments may be strided
views. `C` must satisfy `length(C) ≥ length(A)`.

`rounding` selects the vDSP kernel:

- `:trunc` — round toward zero, like `unsafe_trunc(I, x)` (`vDSP_vfix*` / `vDSP_vfixu*`);
- `:nearest` — round to nearest, ties to even, like `round(I, x)` (`vDSP_vfixr*` /
  `vDSP_vfixru*`).

!!! warning "Out-of-range input is unspecified"
    These kernels do **not** saturate and do **not** throw. For any element whose
    rounded value is not representable in `I` — including `NaN`, `±Inf`, and negative
    values when `I` is unsigned — the result is unspecified: depending on the width and
    on `T` it may wrap, saturate, or be zero, and it differs from both Julia's
    `round(I, x)` (which throws) and `clamp`. Clip first (e.g. with
    [`vclip`](@ref)) if the input is not already known to be in range, and note that
    `Float32` cannot represent `typemax(Int32)`/`typemax(UInt32)` exactly, so clip
    those to the largest `Float32` *below* the limit.

This is a type-dispatched front end over the C-named wrappers (`vfix16`,
`vfixru8`, …). Wraps
[`vDSP_vfix16`](https://developer.apple.com/documentation/accelerate/vdsp_vfix16) and
the rest of its family.
"""
function float_to_int!(C::StridedVector{I}, A::StridedVector{T};
                       rounding::Symbol = :trunc) where {I<:_VDSPInt, T<:Union{Float32,Float64}}
    rounding === :trunc   && return _float_to_int!(C, A, Val(:trunc))
    rounding === :nearest && return _float_to_int!(C, A, Val(:nearest))
    _bad_rounding(rounding)
end

function float_to_int(::Type{I}, A::StridedVector{T};
                      rounding::Symbol = :trunc) where {I<:_VDSPInt, T<:Union{Float32,Float64}}
    float_to_int!(Vector{I}(undef, length(A)), A; rounding)
end

@doc (@doc float_to_int!) float_to_int

"""
    int_to_float!(C::StridedVector{T}, A::StridedVector{I}) -> C
    int_to_float(T, A::StridedVector{I}) -> Vector{T}

Convert the integer vector `A` (`I ∈ (Int8, Int16, Int32, UInt8, UInt16, UInt32)`) to
the floating-point type `T` (`Float32` or `Float64`), equivalent to `T.(A)`. Both
arguments may be strided views. `C` must satisfy `length(C) ≥ length(A)`.

Every conversion is exact except 32-bit integers to `Float32`, which round to nearest
exactly as `Float32(x)` does.

This is a type-dispatched front end over the C-named wrappers (`vflt16`,
`vfltu8`, …). Wraps
[`vDSP_vflt16`](https://developer.apple.com/documentation/accelerate/vdsp_vflt16) and
the rest of its family.
"""
function int_to_float!(C::StridedVector{T}, A::StridedVector{I}) where {T<:Union{Float32,Float64}, I<:_VDSPInt}
    _int_to_float!(C, A)
end

function int_to_float(::Type{T}, A::StridedVector{I}) where {T<:Union{Float32,Float64}, I<:_VDSPInt}
    int_to_float!(Vector{T}(undef, length(A)), A)
end

@doc (@doc int_to_float!) int_to_float

# =====================================================================
# Multi-channel biquad: cross-precision coefficient / target setters
# =====================================================================
#
# `dsp.jl` wraps the matching-precision setters (`SetCoefficientsSingle` on a Float32
# setup, `SetCoefficientsDoubleD` on a Float64 setup). vDSP also provides the two
# cross-precision variants, so that coefficients designed in Float64 can be pushed into
# a Float32 filter (and vice versa) without a temporary converted copy:
#
#     setup     coeffs     C function
#     Float32   Float64    vDSP_biquadm_SetCoefficientsDouble  / SetTargetsDouble
#     Float64   Float32    vDSP_biquadm_SetCoefficientsSingleD / SetTargetsSingleD
#
# The interpolation rate/threshold always take the *setup's* precision.

for (T, S, coefffn, targetfn) in
        ((Float32, Float64, :vDSP_biquadm_SetCoefficientsDouble,  :vDSP_biquadm_SetTargetsDouble),
         (Float64, Float32, :vDSP_biquadm_SetCoefficientsSingleD, :vDSP_biquadm_SetTargetsSingleD))
    @eval begin
        function biquadm_setcoefficients!(setup::BiquadMulti{$T}, coeffs::Vector{$S},
                                          start_sec::Integer, start_chn::Integer,
                                          nsec::Integer, nchn::Integer)
            _check_biquadm_block(setup, length(coeffs), :coeffs, start_sec, start_chn, nsec, nchn)
            GC.@preserve coeffs begin
                LibAccelerate.$coefffn(setup.setup, pointer(coeffs), start_sec, start_chn, nsec, nchn)
            end
            return setup
        end

        function biquadm_settargets!(setup::BiquadMulti{$T}, targets::Vector{$S},
                                     interp_rate::Real, interp_threshold::Real,
                                     start_sec::Integer, start_chn::Integer,
                                     nsec::Integer, nchn::Integer)
            _check_biquadm_block(setup, length(targets), :targets, start_sec, start_chn, nsec, nchn)
            GC.@preserve targets begin
                LibAccelerate.$targetfn(setup.setup, pointer(targets),
                                        $T(interp_rate), $T(interp_threshold),
                                        start_sec, start_chn, nsec, nchn)
            end
            return setup
        end
    end
end

function _check_biquadm_block(setup::BiquadMulti, len::Int, name::Symbol,
                              start_sec, start_chn, nsec, nchn)
    (nsec >= 0 && nchn >= 0) ||
        throw(ArgumentError("nsec and nchn must be non-negative"))
    len >= 5 * nsec * nchn ||
        throw(DimensionMismatch("$name must contain 5*nsec*nchn values (need $(5*nsec*nchn))"))
    (start_sec >= 0 && start_sec + nsec <= setup.sections) ||
        throw(ArgumentError("section range out of bounds for $(setup.sections) sections"))
    (start_chn >= 0 && start_chn + nchn <= setup.channels) ||
        throw(ArgumentError("channel range out of bounds for $(setup.channels) channels"))
    return nothing
end

# =====================================================================
# Fixed-size 16-/32-point FFT, in place on interleaved data (copv)
# =====================================================================
#
# `fft16`/`fft32` in `dsp.jl` go through the split-complex `zopv` kernels, which costs
# five temporary vectors per call. The `copv` kernels take *interleaved* complex data —
# exactly the memory layout of a `Vector{ComplexF32}` — so the mutating forms below run
# with no allocation at all. The header requires "vector-block aligned" (16-byte)
# buffers and allows `Output == Input`.

@inline function _check_fixed_fft(A::StridedVector{ComplexF32}, N::Int, name::Symbol)
    length(A) == N ||
        throw(DimensionMismatch("$name must have length $N; got $(length(A))"))
    _check_unit_stride(A, name)
    UInt(pointer(A)) % 16 == 0 ||
        throw(ArgumentError("$name must be 16-byte aligned for the fixed-size FFT kernels " *
                            "(a view starting at an odd element offset is not); copy it first"))
    return nothing
end

for (N, copv) in ((16, :vDSP_FFT16_copv), (32, :vDSP_FFT32_copv))
    _f! = Symbol("_fft", N, "!")
    @eval function $_f!(out::StridedVector{ComplexF32}, x::StridedVector{ComplexF32}, direction)
        _check_fixed_fft(out, $N, :out)
        _check_fixed_fft(x, $N, :x)
        GC.@preserve out x begin
            pout = pointer(out); px = pointer(x)
            # The kernel permits exact aliasing only; any other overlap is undefined.
            (pout == px || Base.abs(Int(pout) - Int(px)) >= $(N * sizeof(ComplexF32))) ||
                throw(ArgumentError("out and x may be the same array but must not partially overlap"))
            LibAccelerate.$copv(Ptr{Cfloat}(pout), Ptr{Cfloat}(px), direction)
        end
        return out
    end
    for (f!, dir) in ((Symbol("fft", N, "!"), :FFT_FORWARD), (Symbol("bfft", N, "!"), :FFT_INVERSE))
        @eval begin
            $f!(out::StridedVector{ComplexF32}, x::StridedVector{ComplexF32}) = $_f!(out, x, $dir)
            $f!(x::StridedVector{ComplexF32}) = $_f!(x, x, $dir)
        end
    end
end

"""
    fft16!(x)  ;  fft16!(out, x)  ;  bfft16!(x)  ;  bfft16!(out, x)
    fft32!(x)  ;  fft32!(out, x)  ;  bfft32!(x)  ;  bfft32!(out, x)

Allocation-free forms of the fixed-size [`fft16`](@ref) / [`fft32`](@ref) transforms.
`x` and `out` are length-16 (resp. 32) `ComplexF32` vectors; the one-argument form
transforms `x` in place. `fft*!` is the forward transform and `bfft*!` the unnormalized
inverse (`bfft16!(fft16!(x)) ≈ 16 .* x`).

These use vDSP's *interleaved*-complex kernels, which operate directly on the memory of
a `Vector{ComplexF32}` with no split-complex repacking. Unit-stride views are accepted
provided they are 16-byte aligned (any view starting at an even element offset of a
`Vector{ComplexF32}` is); `out` and `x` may be the same array but must not otherwise
overlap.

Wraps [`vDSP_FFT16_copv`](https://developer.apple.com/documentation/accelerate/vdsp_fft16_copv) /
[`vDSP_FFT32_copv`](https://developer.apple.com/documentation/accelerate/vdsp_fft32_copv).
"""
fft16!

@doc (@doc fft16!) bfft16!
@doc (@doc fft16!) fft32!
@doc (@doc fft16!) bfft32!

# =====================================================================
# Coverage note
# =====================================================================

"""
    AppleAccelerate.VDSP_COVERAGE

Documentation of idiomatic-wrapper coverage for `vDSP.h`. Companion to
`AppleAccelerate.VMATH_COVERAGE`.

The numbers come from `gen/coverage_audit.jl`, which inspects the lowered IR of every
method in the package and lists the generated `LibAccelerate.vDSP_*` functions that
nothing calls. Do **not** estimate coverage by grepping `src/` for C symbol names: most
wrappers assemble the name at macro-expansion time
(`Symbol(string("vDSP_vfix", intname, suff))`), so whole families that are fully wrapped
— the 48 `vfix*`/`vflt*` conversions, the fixed-point `_s1_15`/`_s8_24` kernels, every
`D`-suffixed Float64 variant — look "missing" to a text search.

## Tally (macOS 26 SDK raw layer: 469 `vDSP_*` functions)

| Status | Count | Functions |
|--------|------:|-----------|
| Reached by an idiomatic wrapper | 465 | everything not listed below |
| Intentionally left to the raw layer | 4 | `vDSP_DFT_CreateSetup`, `vDSP_DFT_zop`, `vDSP_fft2d_zrip`, `vDSP_fft2d_zripD` |
| Remaining | 0 | — |

The raw layer itself captures ~93% of `vDSP.h`; the functions libclang drops are the
overloads whose signatures use `arm_neon`/`simd` register types (see `gen/README.md`),
which have no array surface.

## Intentionally left to the raw layer

- **`vDSP_DFT_CreateSetup` / `vDSP_DFT_zop`** — the original DFT entry points. The
  header says "We recommend you use `vDSP_DFT_zop_CreateSetup` instead of this
  routine"; that replacement family (`vDSP_DFT_zop_CreateSetup(D)` +
  `vDSP_DFT_Execute(D)`) is what [`plan_dft`](@ref) / [`dft`](@ref) wrap. The old pair
  is Float32-only and computes the identical transform.
- **`vDSP_fft2d_zrip(D)`** — the in-place, no-temporary 2-D real FFT. The same
  transform is wrapped through `vDSP_fft2d_zrop(D)` (the matrix methods of [`rfft`](@ref)) and the
  caller-workspace variants `vDSP_fft2d_zropt(D)` / `vDSP_fft2d_zript(D)`. A Julia
  `Matrix{T}` must be repacked into vDSP's even/odd split-complex layout before any of
  these can run, so "in place on the packed buffer" has no Julian surface of its own,
  and the `zript` form (caller-owned scratch instead of an internal `malloc`) is the
  one worth exposing.

## Closed by this file

- `vDSP_biquadm_SetCoefficientsDouble`, `vDSP_biquadm_SetCoefficientsSingleD`,
  `vDSP_biquadm_SetTargetsDouble`, `vDSP_biquadm_SetTargetsSingleD` — cross-precision
  methods of [`biquadm_setcoefficients!`](@ref) / [`biquadm_settargets!`](@ref).
- `vDSP_FFT16_copv`, `vDSP_FFT32_copv` — allocation-free [`fft16!`](@ref),
  [`bfft16!`](@ref), [`fft32!`](@ref), [`bfft32!`](@ref).
"""
const VDSP_COVERAGE = nothing
