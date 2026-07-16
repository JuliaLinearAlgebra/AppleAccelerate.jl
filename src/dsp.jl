## dsp.jl ##
#
# RAW-LAYER MIGRATION NOTE: the vDSP wrappers here now call the generated
# `LibAccelerate` submodule instead of embedding `ccall` strings. This includes the
# real-valued routines (conv, biquad / biquadm, deq22, desamp, wiener, the *_window
# generators, DCT/DFT execute / destroy) AND, as of the complex/FFT migration, the
# split-complex spectral routines (zaspec, zcoher, ztrans, zcspec) and the full FFT
# family (fft_zop/zip/zrop, 1D and 2D, plus create/destroy_fftsetup). The
# split-complex routines now pass `Ref{DSPSplitComplex}` / `Ref{DSPDoubleSplitComplex}`
# built from the structs that `complexarray.jl` aliases directly to
# `LibAccelerate.DSPSplitComplex` / `LibAccelerate.DSPDoubleSplitComplex`, so the
# `Ref` converts cleanly to the `Ptr{DSPSplitComplex}` the generated wrappers expect.
#
# A few setup CREATORS are still direct `ccall`s on purpose: vDSP_DCT_CreateSetup and
# vDSP_DFT_zop_CreateSetup{,D} return an opaque `Ptr{Cvoid}` and chain a `previous`
# handle; they are kept as-is to avoid churn around the opaque handle types. (The FFT
# setup creators ARE migrated, converting the generated `FFTSetup` pointer result back
# to the `Ptr{Cvoid}` field type.)

mutable struct DFTSetup{T}
    setup::Ptr{Cvoid}
    direction::Int

    function DFTSetup(T::DataType, setup::Ptr{Cvoid}, direction::Int)
        dftsetup = new{T}(setup, direction)
        finalizer(plan_destroy, dftsetup)
        dftsetup
    end
end

mutable struct Biquad{T}
    setup::Ptr{Cvoid}
    sections::Int

    function Biquad(T::DataType, setup::Ptr{Cvoid}, sections::Int)
        biquadsetup = new{T}(setup, sections)
        finalizer(biquaddestroy, biquadsetup)
        biquadsetup
    end
end

## === FFT/DFT constants === ##

const     FFT_FORWARD         = 1
const     FFT_INVERSE         = -1
const     DFT_FORWARD         = 1
const     DFT_INVERSE         = -1
const     SIGNAL_STRIDE       = 1

## === FUNCTIONS == ##

for (T, suff) in ((Float64, "D"), (Float32, ""))

    """
    Convolution between an input Vector{T} 'X', and a kernel/filter Vector{T} 'K'.
    Result vector should have at least length(X) + length(K) - 1 elements

    Returns: 'result'. Computation result is also stored in 'result' argument.
    """
    @eval begin
        function conv!(result::Vector{$T}, X::Vector{$T}, K::Vector{$T})
            ksize = length(K)
            xsize = length(X)
            rsize = length(result)
            if (rsize < xsize + ksize - 1)
                error("'result' must have at least length(X) + length(K) - 1 elements")
            end
            # vDSP reads rsize + ksize - 1 input elements (not just the minimum),
            # so size the zero-padding from rsize, not the natural result length,
            # to keep reads in-bounds for an oversized `result` buffer. Total
            # padded length = rsize + 2*ksize - 1 (front pad ksize-1).
            xpadded::Vector{$T} = [zeros($T, ksize-1); X; zeros($T, rsize + ksize - xsize)]
            # `pointer(K, ksize)` is a bare pointer into K's buffer; root K across
            # the ccall so it cannot be collected mid-call.
            GC.@preserve K LibAccelerate.$(Symbol(string("vDSP_conv", suff)))(xpadded,1,pointer(K, ksize),-1,result,1,rsize,ksize)
            return result
        end
    end


    """
    Convolution between an input Vector{T} 'X', and a kernel/filter Vector{T} 'K'.

    Returns: Vector{T} with length = length(X) + length(K) - 1
    """
    @eval begin
        function conv(X::Vector{$T}, K::Vector{$T})
            result = Array{$T}(undef, length(X) + length(K) - 1)
            conv!(result, X, K)
        end
    end


    """
    In-place convolution between an input Vector{T} 'X', and a kernel/filter Vector{T} 'K'.

    Returns: 'X'. 'X' is overwritten with computation result.
    """
    @eval begin
        function conv!(X::Vector{$T}, K::Vector{$T})
            conv!(X, X, K)
        end
    end

    """
    Cross-correlation of two Vector{T}'s 'X' and 'Y'.
    Result vector should have at least length(X) + length(Y) - 1 elements

    Returns: 'result'. The result of the computation is also stored in 'result'
    """
    @eval begin
        function xcorr!(result::Vector{$T}, X::Vector{$T}, Y::Vector{$T})
            ysize = length(Y)
            xsize = length(X)
            rsize = length(result)
            if (rsize < xsize + ysize - 1)
                error("'result' must have at least length(X) + length(Y) - 1 elements")
            end
            # vDSP reads rsize + ysize - 1 input elements; size the zero-padding
            # from rsize so an oversized `result` buffer stays in-bounds. Total
            # padded length = rsize + 2*ysize - 1 (front pad ysize-1).
            xpadded::Vector{$T} = [zeros($T, ysize-1); X; zeros($T, rsize + ysize - xsize)]
            # Y backs the kernel argument; root it across the ccall.
            GC.@preserve Y LibAccelerate.$(Symbol(string("vDSP_conv", suff)))(xpadded,1,Y,1,result,1,rsize,ysize)
            return result
        end
    end


    """
    Cross-correlation of two Vector{T}'s 'X' and 'Y'.

    Returns: Vector{T} with length(X) + length(Y) - 1
    """
    @eval begin
        function xcorr(X::Vector{$T}, Y::Vector{$T})
            result = Array{$T}(undef, length(X) + length(Y) - 1)
            xcorr!(result, X, Y)
        end
    end


    """
    In-place cross-correlation of two Vector{T}'s 'X' and 'Y'.

    Returns: 'X'. 'X' is overwritten with the result of the cross-correlation.
    """
    @eval begin
        function xcorr!(X::Vector{$T}, Y::Vector{$T})
            xcorr!(X, X, Y)
        end
    end


    """
    Performs auto-correlation of a 1-D signal with itself.

    Returns: Vector{T} with length = 2*length(X) - 1
    """
    @eval begin
        function xcorr(X::Vector{$T})
            xcorr(X, X)
        end
    end
end

"""
    conv(X::Vector{T}, K::Vector{T}) where T <: Union{Float32, Float64}

Compute the convolution of signal `X` with kernel `K` via
[`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv).
Returns a vector of length `length(X) + length(K) - 1`.
"""
conv

"""
    conv!(result::Vector{T}, X::Vector{T}, K::Vector{T})
    conv!(X::Vector{T}, K::Vector{T})

In-place convolution. The 3-argument form stores the result in `result`
(which must have at least `length(X) + length(K) - 1` elements).
The 2-argument form overwrites `X`.
Wraps [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv).
"""
conv!

"""
    xcorr(X::Vector{T}, Y::Vector{T}) where T <: Union{Float32, Float64}
    xcorr(X::Vector{T})

Cross-correlation of `X` and `Y` via
[`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv)
(with reversed kernel). The single-argument form computes auto-correlation.
Returns a vector of length `length(X) + length(Y) - 1`.
"""
xcorr

"""
    xcorr!(result::Vector{T}, X::Vector{T}, Y::Vector{T})
    xcorr!(X::Vector{T}, Y::Vector{T})

In-place cross-correlation. The 3-argument form stores the result in `result`.
The 2-argument form overwrites `X`.
Wraps [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv).
"""
xcorr!

## == Biquadratic/IIR filtering

# Both vDSP_biquad_CreateSetup (Float32) and vDSP_biquad_CreateSetupD (Float64)
# accept Float64 coefficients. The processing type differs: vDSP_biquad processes
# Float32 data, vDSP_biquadD processes Float64 data.

for (T, suff, Dsuff) in ((Float64, "D", "D"), (Float32, "", ""))

    @eval begin
        """
        Initializes a vDSP_biquad_setup for use with vDSP_biquad. A multi-section filter
        can be initialized with a single call to biquad_create_setup. coefficients must
        contain 5 coefficients for each section (as Float64). The three feed-forward
        coefficients are specified first, followed by the two feedback coefficients.

        Returns: Biquad{$($T)}
        """
        function biquadcreate(coefficients::Vector{Float64}, sections::Int, ::Type{$T})
            if length(coefficients) < 5*sections
                error("Incomplete biquad specification provided - coefficients must
                            contain 5 elements for each filter section")
            end
            setup = ccall(($(string("vDSP_biquad_CreateSetup", Dsuff), libacc)),  Ptr{Cvoid},
                          (Ptr{Float64}, UInt64),
                          coefficients, sections)
            setup == C_NULL && error("vDSP_biquad_CreateSetup$($Dsuff) failed to create a setup for $sections section(s)")
            return Biquad($T, setup, sections)
        end
    end

    @eval begin
        """
        Filters an input array X with the specified Biquad filter and filter delay values provided
        in delays; only numelem elements are filtered. After execution, delays contains the final
        state data of the filter.

        Returns: Vector{$($T)}
        """
        function biquad(X::Vector{$T}, delays::Vector{$T}, numelem::Int, biquad::Biquad{$T})
            if length(delays) < (2*(biquad.sections)+2)
                error("Incomplete delay specification provided - delays must contain 2*M + 2
                                values where M is the number of sections in the biquad")
            end
            if numelem > length(X)
                error("numelem = $numelem exceeds the input length $(length(X))")
            end
            result::Vector{$T} = similar(X)
            LibAccelerate.$(Symbol(string("vDSP_biquad", suff)))(biquad.setup,delays,X,1,result,1,numelem)
            return result
        end
    end

    @eval begin
        """
        Frees all resources associated with a particular Biquad previously
        created through a call to biquad_create_setup. This is called automatically
        when the setup object is no longer visible to the garbage collector.

        Returns: Cvoid
        """
        function biquaddestroy(biquad::Biquad{$T})
            LibAccelerate.$(Symbol(string("vDSP_biquad_DestroySetup", Dsuff)))(biquad.setup)
        end
    end
end

# Backward-compatible method: biquadcreate(::Vector{Float64}, ::Int) defaults to Float64
biquadcreate(coefficients::Vector{Float64}, sections::Int) = biquadcreate(coefficients, sections, Float64)

"""
    biquadcreate(coefficients::Vector{Float64}, sections::Int, [T=Float64])

Create a biquad IIR filter setup. `coefficients` must contain 5 values per section
(3 feed-forward + 2 feedback). Returns a `Biquad{T}` setup object.
Wraps [`vDSP_biquad_CreateSetup`](https://developer.apple.com/documentation/accelerate/vdsp_biquad_createsetup).
"""
biquadcreate

"""
    biquad(X, delays, numelem, bq::Biquad{T})

Apply a cascaded biquad IIR filter to input `X`. `delays` holds the filter state
and is updated in-place. Returns the filtered output vector.
Wraps [`vDSP_biquad`](https://developer.apple.com/documentation/accelerate/vdsp_biquad).
"""
biquad

"""
    biquaddestroy(bq::Biquad{T})

Free resources associated with a biquad setup. Called automatically by the finalizer.
Wraps [`vDSP_biquad_DestroySetup`](https://developer.apple.com/documentation/accelerate/vdsp_biquad_destroysetup).
"""
biquaddestroy

## == Multi-channel Biquadratic/IIR filtering

mutable struct BiquadMulti{T}
    setup::Ptr{Cvoid}
    channels::Int
    sections::Int

    function BiquadMulti(T::DataType, setup::Ptr{Cvoid}, channels::Int, sections::Int)
        obj = new{T}(setup, channels, sections)
        finalizer(biquadm_destroy, obj)
        obj
    end
end

# Both vDSP_biquadm_CreateSetup (Float32) and vDSP_biquadm_CreateSetupD (Float64)
# accept Float64 coefficients.
for (T, suff, Dsuff) in ((Float32, "", ""), (Float64, "D", "D"))
    @eval begin
        """
        Creates a multi-channel biquad IIR filter setup. `coefficients` must contain
        `5 * channels * sections` Float64 values. Delay values are initialized to zero.

        The coefficients are laid out section-major: for each section (outer), the
        five values `[b0, b1, b2, a1, a2]` are given for every channel (inner). That
        is, element `5 * (section * channels + channel)` begins the block for a given
        `(section, channel)`. With a single section this is just the per-channel
        blocks back to back.

        Returns: BiquadMulti{$($T)}
        """
        function biquadm_create(coefficients::Vector{Float64}, channels::Int, sections::Int, ::Type{$T})
            if length(coefficients) < 5 * channels * sections
                error("Incomplete biquadm specification - coefficients must contain 5*channels*sections elements")
            end
            setup = ccall(($(string("vDSP_biquadm_CreateSetup", Dsuff)), libacc), Ptr{Cvoid},
                          (Ptr{Float64}, UInt64, UInt64),
                          coefficients, sections, channels)
            setup == C_NULL && error("vDSP_biquadm_CreateSetup$($Dsuff) failed to create a setup for $channels channel(s) and $sections section(s)")
            return BiquadMulti($T, setup, channels, sections)
        end

        """
        Apply a multi-channel biquad IIR filter. `X` is a vector of input channel vectors,
        each of length `numelem`. Returns a vector of output channel vectors.

        Returns: Vector{Vector{$($T)}}
        """
        function biquadm(X::Vector{Vector{$T}}, numelem::Int, setup::BiquadMulti{$T})
            M = setup.channels
            length(X) == M || error("Expected $M channels, got $(length(X))")
            all(x -> length(x) >= numelem, X) ||
                error("Each input channel must contain at least numelem = $numelem elements")
            Y = [similar(x) for x in X]
            # Build the pointer arrays inside GC.@preserve so the channel buffers
            # (and the Ptr arrays passed as Ptr{Ptr{T}}) stay rooted for the whole
            # ccall. Taking the pointers before the preserve region left them
            # unprotected and intermittently produced an unwritten channel.
            xptrs = Vector{Ptr{$T}}(undef, M)
            yptrs = Vector{Ptr{$T}}(undef, M)
            GC.@preserve X Y xptrs yptrs begin
                for i in 1:M
                    xptrs[i] = pointer(X[i])
                    yptrs[i] = pointer(Y[i])
                end
                LibAccelerate.$(Symbol(string("vDSP_biquadm", suff)))(setup.setup,xptrs,1,yptrs,1,numelem)
            end
            return Y
        end

        function biquadm_destroy(setup::BiquadMulti{$T})
            LibAccelerate.$(Symbol(string("vDSP_biquadm_DestroySetup", Dsuff)))(setup.setup)
        end
    end
end

# Convenience: default to Float32
biquadm_create(coefficients::Vector{Float64}, channels::Int, sections::Int) =
    biquadm_create(coefficients, channels, sections, Float32)

"""
    biquadm_create(coefficients, channels, sections, [T=Float32])

Create a multi-channel biquad IIR filter setup. `coefficients` must contain
`5 * channels * sections` Float64 values, laid out section-major: for each
section, the `[b0, b1, b2, a1, a2]` block is given for every channel (so element
`5 * (section * channels + channel)` starts a given `(section, channel)` block).
Returns a `BiquadMulti{T}` setup object.
Wraps [`vDSP_biquadm_CreateSetup`](https://developer.apple.com/documentation/accelerate/vdsp_biquadm_createsetup).
"""
biquadm_create

"""
    biquadm(X, numelem, setup::BiquadMulti{T})

Apply a multi-channel biquad IIR filter. `X` is a vector of per-channel input vectors.
Returns a vector of per-channel output vectors.
Wraps [`vDSP_biquadm`](https://developer.apple.com/documentation/accelerate/vdsp_biquadm).
"""
biquadm

## == Spectral Analysis == ##

for (T, suff, SC) in ((Float32, "", :DSPSplitComplex), (Float64, "D", :DSPDoubleSplitComplex))
    @eval begin
        """
        Accumulating autospectrum: `C[n] += |A[n]|^2`. `A` is a complex vector,
        `C` is a real vector that accumulates the power spectrum.

        Wraps [`vDSP_zaspec`](https://developer.apple.com/documentation/accelerate/vdsp_zaspec).
        """
        function zaspec!(C::Vector{$T}, A::Vector{Complex{$T}})
            n = length(A)
            length(C) >= n || throw(DimensionMismatch(
                "zaspec!: output `C` must have at least length(A) = $n elements; got $(length(C))"))
            realp = $T.(real.(A))
            imagp = $T.(imag.(A))
            GC.@preserve realp imagp begin
                sc = $SC(pointer(realp), pointer(imagp))
                LibAccelerate.$(Symbol(string("vDSP_zaspec", suff)))(Ref(sc), C, n)
            end
            return C
        end

        function zaspec(A::Vector{Complex{$T}})
            C = zeros($T, length(A))
            zaspec!(C, A)
        end

        """
        Coherence function: `D[n] = |C[n]|^2 / (A[n] * B[n])`. `A` and `B` are real
        power spectra, `C` is a complex cross-spectrum, `D` is the output coherence.

        Wraps [`vDSP_zcoher`](https://developer.apple.com/documentation/accelerate/vdsp_zcoher).
        """
        function zcoher!(D::Vector{$T}, A::Vector{$T}, B::Vector{$T}, C::Vector{Complex{$T}})
            n = length(A)
            (length(B) >= n && length(C) >= n && length(D) >= n) || throw(DimensionMismatch(
                "zcoher!: `B`, `C`, and `D` must each have at least length(A) = $n elements; " *
                "got B=$(length(B)), C=$(length(C)), D=$(length(D))"))
            cr = $T.(real.(C))
            ci = $T.(imag.(C))
            GC.@preserve cr ci begin
                sc = $SC(pointer(cr), pointer(ci))
                LibAccelerate.$(Symbol(string("vDSP_zcoher", suff)))(A, B, Ref(sc), D, n)
            end
            return D
        end

        function zcoher(A::Vector{$T}, B::Vector{$T}, C::Vector{Complex{$T}})
            D = Vector{$T}(undef, length(A))
            zcoher!(D, A, B, C)
        end

        """
        Transfer function: `C[n] = B[n] / A[n]`. `A` is a real power spectrum,
        `B` is a complex cross-spectrum, `C` is the output complex transfer function.

        Wraps [`vDSP_ztrans`](https://developer.apple.com/documentation/accelerate/vdsp_ztrans).
        """
        function ztrans!(C::Vector{Complex{$T}}, A::Vector{$T}, B::Vector{Complex{$T}})
            n = length(A)
            (length(B) >= n && length(C) >= n) || throw(DimensionMismatch(
                "ztrans!: `B` and `C` must each have at least length(A) = $n elements; " *
                "got B=$(length(B)), C=$(length(C))"))
            br = $T.(real.(B))
            bi = $T.(imag.(B))
            cr = Vector{$T}(undef, n)
            ci = Vector{$T}(undef, n)
            GC.@preserve br bi cr ci begin
                scb = $SC(pointer(br), pointer(bi))
                scc = $SC(pointer(cr), pointer(ci))
                LibAccelerate.$(Symbol(string("vDSP_ztrans", suff)))(A, Ref(scb), Ref(scc), n)
            end
            @inbounds for i in 1:n
                C[i] = complex(cr[i], ci[i])
            end
            return C
        end

        function ztrans(A::Vector{$T}, B::Vector{Complex{$T}})
            C = Vector{Complex{$T}}(undef, length(A))
            ztrans!(C, A, B)
        end

        """
        Accumulating cross-spectrum: `C[n] += conj(A[n]) * B[n]`. `A` and `B` are
        complex vectors, `C` is the complex output that accumulates.

        Wraps [`vDSP_zcspec`](https://developer.apple.com/documentation/accelerate/vdsp_zcspec).
        """
        function zcspec!(C::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            n = length(A)
            (length(B) >= n && length(C) >= n) || throw(DimensionMismatch(
                "zcspec!: `B` and `C` must each have at least length(A) = $n elements; " *
                "got B=$(length(B)), C=$(length(C))"))
            ar = $T.(real.(A))
            ai = $T.(imag.(A))
            br = $T.(real.(B))
            bi = $T.(imag.(B))
            cr = $T.(real.(C))
            ci = $T.(imag.(C))
            GC.@preserve ar ai br bi cr ci begin
                sca = $SC(pointer(ar), pointer(ai))
                scb = $SC(pointer(br), pointer(bi))
                scc = $SC(pointer(cr), pointer(ci))
                LibAccelerate.$(Symbol(string("vDSP_zcspec", suff)))(Ref(sca), Ref(scb), Ref(scc), n)
            end
            @inbounds for i in 1:n
                C[i] = complex(cr[i], ci[i])
            end
            return C
        end

        function zcspec(A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            C = zeros(Complex{$T}, length(A))
            zcspec!(C, A, B)
        end
    end
end

"""
    zaspec(A::Vector{Complex{T}}) -> Vector{T}

Autospectrum (power spectrum): returns a real vector `C` where `C[n] = |A[n]|^2`.
Wraps [`vDSP_zaspec`](https://developer.apple.com/documentation/accelerate/vdsp_zaspec).

See also [`zaspec!`](@ref) for the accumulating in-place variant.
"""
zaspec

"""
    zcspec(A::Vector{Complex{T}}, B::Vector{Complex{T}}) -> Vector{Complex{T}}

Cross-spectrum: returns a complex vector `C` where `C[n] = conj(A[n]) * B[n]`.
Wraps [`vDSP_zcspec`](https://developer.apple.com/documentation/accelerate/vdsp_zcspec).

See also [`zcspec!`](@ref) for the accumulating in-place variant.
"""
zcspec

"""
    zcoher(A::Vector{T}, B::Vector{T}, C::Vector{Complex{T}}) -> Vector{T}

Coherence function: returns a real vector `D` where `D[n] = |C[n]|^2 / (A[n] * B[n])`.
`A` and `B` are real power spectra, `C` is a complex cross-spectrum.
Wraps [`vDSP_zcoher`](https://developer.apple.com/documentation/accelerate/vdsp_zcoher).

See also [`zcoher!`](@ref) for the in-place variant.
"""
zcoher

"""
    ztrans(A::Vector{T}, B::Vector{Complex{T}}) -> Vector{Complex{T}}

Transfer function: returns a complex vector `C` where `C[n] = B[n] / A[n]`.
`A` is a real power spectrum, `B` is a complex cross-spectrum.
Wraps [`vDSP_ztrans`](https://developer.apple.com/documentation/accelerate/vdsp_ztrans).

See also [`ztrans!`](@ref) for the in-place variant.
"""
ztrans

## == Recursive Filter, FIR Decimation, Wiener-Levinson == ##

for (T, suff) in ((Float32, ""), (Float64, "D"))

    @eval begin
        """
            deq22!(C::Vector{$($T)}, A::Vector{$($T)}, B::Vector{$($T)})

        Second-order (two-pole two-zero) recursive filter using `vDSP_deq22`.
        `A` has N+2 elements (A[1:2] are initial state), `B` has 5 coefficients,
        `C` has N+2 elements (C[1:2] must be preset as initial output state).
        Computes: `C[n] = A[n]*B[1] + A[n-1]*B[2] + A[n-2]*B[3] - C[n-1]*B[4] - C[n-2]*B[5]`

        Returns: `C`
        """
        function deq22!(C::Vector{$T}, A::Vector{$T}, B::Vector{$T})
            length(B) == 5 || error("B must have exactly 5 coefficients")
            length(A) >= 3 || error("A must have at least 3 elements (2 state + 1 sample)")
            length(C) == length(A) || error("C must have the same length as A")
            N = UInt64(length(A) - 2)
            LibAccelerate.$(Symbol(string("vDSP_deq22", suff)))(A,1,B,C,1,N)
            return C
        end

        """
            deq22(A::Vector{$($T)}, B::Vector{$($T)})

        Allocating version of `deq22!`. Pads `A` with 2 leading zeros and returns
        only the N output samples (without the 2-element state prefix).

        Returns: `Vector{$($T)}` of length `length(A)`
        """
        function deq22(A::Vector{$T}, B::Vector{$T})
            N = length(A)
            Apad = [$T(0); $T(0); A]
            C = zeros($T, N + 2)
            deq22!(C, Apad, B)
            return C[3:end]
        end
    end

    @eval begin
        """
            desamp!(C::Vector{$($T)}, A::Vector{$($T)}, DF::Int, F::Vector{$($T)})

        FIR decimation filter using `vDSP_desamp`. Filters input `A` with FIR
        coefficients `F` (P taps) and decimation factor `DF`.
        `C` must have at least `div(length(A) - P, DF) + 1` elements.
        Computes: `C[n] = sum(A[n*DF+p] * F[p] for p in 0:P-1)` (0-indexed)

        Returns: `C`
        """
        function desamp!(C::Vector{$T}, A::Vector{$T}, DF::Int, F::Vector{$T})
            # Compute the output count in SIGNED arithmetic; converting to
            # UInt64 first would underflow when length(F) > length(A).
            length(A) >= length(F) ||
                error("length(A) ($(length(A))) must be >= length(F) ($(length(F)))")
            P = length(F)
            nout = div(length(A) - P, DF) + 1
            length(C) >= nout || error("C must have at least $(nout) elements")
            LibAccelerate.$(Symbol(string("vDSP_desamp", suff)))(A, DF, F, C, UInt64(nout), UInt64(P))
            return C
        end

        """
            desamp(A::Vector{$($T)}, DF::Int, F::Vector{$($T)})

        Allocating version of `desamp!`. Returns a vector of length
        `div(length(A) - length(F), DF) + 1`.

        Returns: `Vector{$($T)}`
        """
        function desamp(A::Vector{$T}, DF::Int, F::Vector{$T})
            length(A) >= length(F) ||
                error("length(A) ($(length(A))) must be >= length(F) ($(length(F)))")
            P = length(F)
            Nout = div(length(A) - P, DF) + 1
            C = Vector{$T}(undef, Nout)
            desamp!(C, A, DF, F)
        end
    end

    @eval begin
        """
            wiener!(F::Vector{$($T)}, P::Vector{$($T)}, A::Vector{$($T)}, C::Vector{$($T)}; flag::Int=0)

        Wiener-Levinson filter using `vDSP_wiener`. Solves the Wiener-Hopf equation
        `R * F = C` where `R` is the Toeplitz autocorrelation matrix formed from `A`.
        `A` = autocorrelation coefficients (L elements), `C` = cross-correlation (L),
        `F` = output filter coefficients (L), `P` = workspace (L).

        Returns: `(F, error_code)` where error_code is 0 on success.
        """
        function wiener!(F::Vector{$T}, P::Vector{$T}, A::Vector{$T}, C::Vector{$T}; flag::Int=0)
            L = UInt64(length(A))
            length(C) == L || error("C must have the same length as A")
            length(F) >= L || error("F must have at least L elements")
            length(P) >= L || error("P (workspace) must have at least L elements")
            err = Ref{Cint}(0)
            LibAccelerate.$(Symbol(string("vDSP_wiener", suff)))(L,A,C,F,P,Cint(flag),err)
            return (F, Int(err[]))
        end

        """
            wiener(A::Vector{$($T)}, C::Vector{$($T)}; flag::Int=0)

        Allocating version of `wiener!`. Returns `(F, error_code)`.
        """
        function wiener(A::Vector{$T}, C::Vector{$T}; flag::Int=0)
            L = length(A)
            F = Vector{$T}(undef, L)
            P = Vector{$T}(undef, L)
            wiener!(F, P, A, C; flag=flag)
        end
    end
end

## == WINDOW GENERATION == ##

"""
    blackman(length, [rtype=Float64])

Generate a Blackman window of the given length.
Wraps [`vDSP_blkman_window`](https://developer.apple.com/documentation/accelerate/vdsp_blkman_window).
"""
function blackman(length::Int, rtype::DataType=Float64)
    result::Vector{rtype} = Array{rtype}(undef, length)
    blackman!(result, length, 0)
end

"""
    hamming(length, [rtype=Float64])

Generate a Hamming window of the given length.
Wraps [`vDSP_hamm_window`](https://developer.apple.com/documentation/accelerate/vdsp_hamm_window).
"""
function hamming(length::Int, rtype::DataType=Float64)
    result::Vector{rtype} = Array{rtype}(undef, length)
    hamming!(result, length, 0)
end

"""
    hanning(length, [rtype=Float64])

Generate a Hanning window of the given length.
Wraps [`vDSP_hann_window`](https://developer.apple.com/documentation/accelerate/vdsp_hann_window).
"""
function hanning(length::Int, rtype::DataType=Float64)
    result::Vector{rtype} = Array{rtype}(undef, length)
    hanning!(result, length, 0)
end

for (T, suff) in ((Float32, ""), (Float64, "D"))

    """
    Generates a Blackman window of length 'length' and stores it in `result'. By
    default, a full window is created, but if flag=1, only the first (n+1)/2 points
    will be calculated.

    Returns: Vector{$T}
    """
    @eval begin
        function blackman!(result::Vector{$T},  length::Int, flag::Int=0)
            LibAccelerate.$(Symbol(string("vDSP_blkman_window", suff)))(result,length,flag)
            return result
        end
    end


    """
    Generates a Hamming window of length 'length' and stores it in `result'. By
    default, a full window is created, but if flag=1, only the first (n+1)/2 points
    will be calculated.

    Returns: Vector{$T}
    """
    @eval begin
        function hamming!(result::Vector{$T},  length::Int, flag::Int=0)
            LibAccelerate.$(Symbol(string("vDSP_hamm_window", suff)))(result,length,flag)
            return result
        end
    end

    """
    Generates a Hanning window of length 'length' and stores it in `result'.
    Different window options can be set using the flag argument; default value
    is zero (denormalized Hanning window)

    flag = 0: Returns a denormalized Hanning window
    flag = 1: Returns a denormalized Hanning window with only (n+1)/2 points
    flag = 2: Returns a normalized Hanning window
    flag = 3: Returns a normalized Hanning window with only (n+1)/2 points

    Returns: Vector{$T}
    """
    @eval begin
        function hanning!(result::Vector{$T},  length::Int, flag::Int=0)
            LibAccelerate.$(Symbol(string("vDSP_hann_window", suff)))(result,length,flag)
            return result
        end
    end

end


## == Discrete Cosine Transform (DCT) & Discrete Fourier Transform (DFT) == ##
"""
    plan_dct(length, dct_type, [previous])

Create a DCT setup object. `dct_type` must be 2, 3, or 4 (Type II, III, IV).
Length must be `f * 2^n` where `f ∈ {1,3,5,15}` and `n ≥ 4`.
Wraps [`vDSP_DCT_CreateSetup`](https://developer.apple.com/documentation/accelerate/vdsp_dct_createsetup).
"""
function plan_dct(length::Int,  dct_type::Int, previous=C_NULL)
    n = trailing_zeros(length)
    f = length >> n
    if dct_type < 2 ||  dct_type > 4
        error("DCT type ", dct_type, " is not supported. Only DCT types 2, 3 and 4 are supported")
    elseif !(n >= 4 && f in (1,3,5,15))
        error("Invalid DCT length. Length must be equal to f*(2^n) where f = 1,3,5,15 and n >= 4")
    end
    setup::Ptr{Cvoid} = ccall(("vDSP_DCT_CreateSetup", libacc), Ptr{Cvoid},
                             (Ptr{Cvoid}, UInt64, UInt64),
                             previous, length, dct_type)
    setup == C_NULL && error("vDSP_DCT_CreateSetup failed to create a setup for length $length and DCT type $dct_type")
    return DFTSetup(Float32, setup, 0)
end


"""
    dct(X::Vector{Float32}, setup::DFTSetup)
    dct(X::Vector{Float32}, [dct_type=2])

Compute the Discrete Cosine Transform of `X`.
Wraps [`vDSP_DCT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dct_execute).

!!! note
    Apple's Accelerate framework provides only a single-precision DCT
    (`vDSP_DCT_CreateSetup`/`vDSP_DCT_Execute`); there is no `Float64` variant
    (no `vDSP_DCT_CreateSetupD`/`vDSP_DCT_ExecuteD` exists in vDSP). Calling
    `dct` or `idct` with a `Vector{Float64}` therefore throws an `ArgumentError`.
    Convert the input to `Float32`, or use FFTW.jl for a double-precision DCT.
"""
function dct(X::Vector{Float32}, setup::DFTSetup)
    result = similar(X)
    LibAccelerate.vDSP_DCT_Execute(setup.setup,X,result)
    return result
end


function dct(X::Vector{Float32}, dct_type::Int=2)
    setup = plan_dct(length(X), dct_type)
    return dct(X, setup)
end

@noinline _dct_no_float64() = throw(ArgumentError(string(
    "Apple's Accelerate framework provides only a single-precision DCT ",
    "(vDSP_DCT_CreateSetup/vDSP_DCT_Execute); no Float64 variant exists in vDSP. ",
    "Convert the input to Float32 (e.g. dct(Float32.(x))) or use FFTW.jl for a ",
    "double-precision DCT.")))

dct(X::Vector{Float64}, setup::DFTSetup) = _dct_no_float64()
dct(X::Vector{Float64}, dct_type::Int=2) = _dct_no_float64()

"""
    idct(X::Vector{Float32})

Compute the inverse of the (unnormalized) DCT-II computed by [`dct`](@ref), so that
`idct(dct(x)) ≈ x`. Uses the DCT-II/DCT-III duality: vDSP's DCT-III applied to a
DCT-II spectrum reproduces the input scaled by `N/2`, so this returns
`dct(X, 3) .* (2/N)`.
Wraps [`vDSP_DCT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dct_execute).

Not available for `Float64` inputs: see the note in [`dct`](@ref).
"""
idct(X::Vector{Float32}) = dct(X, 3) .* (2.0f0 / length(X))
idct(X::Vector{Float64}) = _dct_no_float64()


"""
    plan_destroy(setup::DFTSetup)

Destroy a DCT/DFT setup object, freeing its resources.
Wraps [`vDSP_DFT_DestroySetup`](https://developer.apple.com/documentation/accelerate/vdsp_dft_destroysetup).
"""
function plan_destroy(setup::DFTSetup{Float32})
    LibAccelerate.vDSP_DFT_DestroySetup(setup.setup)
end

function plan_destroy(setup::DFTSetup{Float64})
    LibAccelerate.vDSP_DFT_DestroySetupD(setup.setup)
end


# --- Complex DFT (non-power-of-2 support) ---

"""
    plan_dft(length::Int, direction::Int, ::Type{T}=Float32; previous=C_NULL) where T

Create a DFT setup for complex-to-complex DFT of the given `length` and `direction`
(`DFT_FORWARD` or `DFT_INVERSE`). Length must be `f * 2^n` where `f ∈ {1, 3, 5, 15}`
and `n ≥ 3`. Optionally pass a `previous` setup to share underlying data.
Wraps [`vDSP_DFT_zop_CreateSetup`](https://developer.apple.com/documentation/accelerate/vdsp_dft_zop_createsetup).

Returns: `DFTSetup{T}`
"""
function plan_dft(length::Int, direction::Int, ::Type{Float32}=Float32; previous=C_NULL)
    setup = ccall(("vDSP_DFT_zop_CreateSetup", libacc), Ptr{Cvoid},
                  (Ptr{Cvoid}, UInt64, Cint),
                  previous isa DFTSetup ? previous.setup : previous, length, Cint(direction))
    setup == C_NULL && error("Invalid DFT length $length. Length must be f*(2^n) where f ∈ {1,3,5,15} and n ≥ 3")
    return DFTSetup(Float32, setup, direction)
end

function plan_dft(length::Int, direction::Int, ::Type{Float64}; previous=C_NULL)
    setup = ccall(("vDSP_DFT_zop_CreateSetupD", libacc), Ptr{Cvoid},
                  (Ptr{Cvoid}, UInt64, Cint),
                  previous isa DFTSetup ? previous.setup : previous, length, Cint(direction))
    setup == C_NULL && error("Invalid DFT length $length. Length must be f*(2^n) where f ∈ {1,3,5,15} and n ≥ 3")
    return DFTSetup(Float64, setup, direction)
end

"""
    dft(Ir::Vector{T}, Ii::Vector{T}, setup::DFTSetup{T})

Execute the complex DFT defined by `setup` on split-complex input (`Ir`, `Ii`).
Wraps [`vDSP_DFT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dft_execute).

Returns: `(Or, Oi)` — real and imaginary parts of the output.
"""
function dft(Ir::Vector{Float32}, Ii::Vector{Float32}, setup::DFTSetup{Float32})
    n = length(Ir)
    Or = Vector{Float32}(undef, n)
    Oi = Vector{Float32}(undef, n)
    LibAccelerate.vDSP_DFT_Execute(setup.setup,Ir,Ii,Or,Oi)
    return (Or, Oi)
end

function dft(Ir::Vector{Float64}, Ii::Vector{Float64}, setup::DFTSetup{Float64})
    n = length(Ir)
    Or = Vector{Float64}(undef, n)
    Oi = Vector{Float64}(undef, n)
    LibAccelerate.vDSP_DFT_ExecuteD(setup.setup,Ir,Ii,Or,Oi)
    return (Or, Oi)
end

"""
    dft(X::Vector{Complex{T}}, setup::DFTSetup{T})

Execute the complex DFT on an interleaved complex vector.
Wraps [`vDSP_DFT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dft_execute).

Returns: `Vector{Complex{T}}`
"""
function dft(X::Vector{Complex{T}}, setup::DFTSetup{T}) where {T<:Union{Float32,Float64}}
    Ir = T.(real.(X))
    Ii = T.(imag.(X))
    Or, Oi = dft(Ir, Ii, setup)
    return complex.(Or, Oi)
end

"""
    dft(X::Vector{Complex{T}}, direction::Int) where T
    dft(X::Vector{Complex{T}}) where T

Compute the DFT of `X`, auto-creating a setup. Default direction is forward.
Wraps [`vDSP_DFT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dft_execute).

Returns: `Vector{Complex{T}}`
"""
function dft(X::Vector{Complex{T}}, direction::Int) where {T<:Union{Float32,Float64}}
    setup = _cached_dftsetup(T, length(X), direction)
    return dft(X, setup)
end

function dft(X::Vector{Complex{T}}) where {T<:Union{Float32,Float64}}
    dft(X, DFT_FORWARD)
end

"""
    idft(X::Vector{Complex{T}}, setup::DFTSetup{T})
    idft(X::Vector{Complex{T}})

Compute the normalized inverse DFT of `X`. The setup must have been created
with `DFT_INVERSE` direction. Returns `dft(X, setup) ./ length(X)`.
Wraps [`vDSP_DFT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dft_execute).

Returns: `Vector{Complex{T}}`
"""
function idft(X::Vector{Complex{T}}, setup::DFTSetup{T}) where {T<:Union{Float32,Float64}}
    return dft(X, setup) ./ length(X)
end

function idft(X::Vector{Complex{T}}) where {T<:Union{Float32,Float64}}
    setup = _cached_dftsetup(T, length(X), DFT_INVERSE)
    return dft(X, setup) ./ length(X)
end


mutable struct FFTSetup{T}
    plan::Ptr{Cvoid}

    function FFTSetup{Float64}(n::Integer, radix::Integer = 2)
        @assert ispow2(n) "n must be a power of 2"
        logn = trailing_zeros(n)
        plan = Ptr{Cvoid}(LibAccelerate.vDSP_create_fftsetupD(logn, radix))
        setup = new{Float64}(plan)
        finalizer(destroy_fftsetup, setup)
        setup
    end

    function FFTSetup{Float32}(n::Integer, radix::Integer = 2)
        @assert ispow2(n) "n must be a power of 2"
        logn = trailing_zeros(n)
        plan = Ptr{Cvoid}(LibAccelerate.vDSP_create_fftsetup(logn, radix))
        setup = new{Float32}(plan)
        finalizer(destroy_fftsetup, setup)
        setup
    end
end

"""
    plan_fft(x::VecOrMat{Complex{T}}) where T <: Union{Float32, Float64}
    plan_fft(n::Integer, [T=Float64], [radix=2])

Create a reusable FFT setup object for repeated transforms of the same size.
Wraps [`vDSP_create_fftsetup`](https://developer.apple.com/documentation/accelerate/vdsp_create_fftsetup) /
[`vDSP_create_fftsetupD`](https://developer.apple.com/documentation/accelerate/vdsp_create_fftsetupd).
The setup is automatically destroyed when garbage collected.
"""
function plan_fft(n::Integer, ::Type{T}=Float64, radix::Integer = 2) where T <: Union{Float32, Float64}
    FFTSetup{T}(n, radix)
end

function plan_fft(x::Vector{Complex{T}}) where T <: Union{Float32, Float64}
    FFTSetup{T}(length(x))
end

function plan_fft(x::Matrix{Complex{T}}) where T <: Union{Float32, Float64}
    nrows, ncols = size(x)
    FFTSetup{T}(max(nrows, ncols))
end

function destroy_fftsetup(setup::FFTSetup{Float64})
    LibAccelerate.vDSP_destroy_fftsetupD(setup.plan)
end

function destroy_fftsetup(setup::FFTSetup{Float32})
    LibAccelerate.vDSP_destroy_fftsetup(setup.plan)
end

# --- Global setup caches for the convenience (no-plan) API ---
#
# The convenience entry points (`fft(x)`, `ifft(x)`, `rfft(x)`, ...) used to create a
# fresh setup on every call, which is expensive: `vDSP_create_fftsetup` builds sizable
# twiddle-factor tables and can dominate the cost of small transforms. Apple documents
# FFT setup structures as reusable for any number of transforms of equal or smaller
# size, and as read-only during transform execution, so a single cached setup per
# (precision, log2 size, radix) can safely be shared — including across threads.
# Cache growth is bounded by the number of distinct transform sizes used, mirroring
# FFTW's plan caching behavior. The explicit `plan_fft`/`plan_rfft`/`plan_dft` API is
# unchanged and still returns fresh, independently finalized setups.

const _FFT_SETUP_CACHE = Dict{Tuple{DataType,Int,Int},FFTSetup}()
const _DFT_SETUP_CACHE = Dict{Tuple{DataType,Int,Int},DFTSetup}()
const _SETUP_CACHE_LOCK = ReentrantLock()

# Shared, cached FFTSetup{T} for power-of-2 transforms of size n (keyed by log2n, radix).
function _cached_fftsetup(::Type{T}, n::Integer, radix::Integer=2) where {T<:Union{Float32,Float64}}
    @assert ispow2(n) "n must be a power of 2"
    key = (T, trailing_zeros(n), Int(radix))
    return lock(_SETUP_CACHE_LOCK) do
        get!(() -> FFTSetup{T}(n, radix), _FFT_SETUP_CACHE, key)
    end::FFTSetup{T}
end

_cached_fftsetup(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} =
    _cached_fftsetup(T, length(x))
# 2D transforms only need a setup large enough for the bigger dimension (setups
# support all transforms of equal or smaller size), matching `plan_fft(::Matrix)`.
_cached_fftsetup(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} =
    _cached_fftsetup(T, max(size(x)...))

# Shared, cached DFTSetup{T} for mixed-radix (f*2^k) transforms, keyed by length and direction.
function _cached_dftsetup(::Type{T}, n::Integer, direction::Int) where {T<:Union{Float32,Float64}}
    key = (T, Int(n), direction)
    return lock(_SETUP_CACHE_LOCK) do
        get!(() -> plan_dft(Int(n), direction, T), _DFT_SETUP_CACHE, key)
    end::DFTSetup{T}
end

# --- Mixed-radix length support (for non-power-of-2 1D complex FFT) ---
#
# Apple's vDSP DFT (vDSP_DFT_zop_CreateSetup) supports lengths of the form
# `f * 2^n` where `f ∈ {1, 3, 5, 15}`. `fft`/`ifft`/`bfft` on a 1D vector with no
# explicit `FFTSetup` transparently route to this mixed-radix DFT path when the
# length is not a power of two but is a supported `f * 2^n`; power-of-2 lengths
# keep using the (faster, in-place-capable) `vDSP_fft_*` path unchanged.

# Return the odd cofactor `f` of `n` (i.e. `n` with all factors of 2 removed).
_odd_cofactor(n::Integer) = n >> trailing_zeros(n)

"""
    is_supported_fft_length(n::Integer) -> Bool

Return `true` if a 1D complex FFT of length `n` is supported by Apple vDSP via the
idiomatic `fft`/`ifft`/`bfft` API. Supported lengths are `f * 2^k` with
`f ∈ {1, 3, 5, 15}` (this includes all powers of two).
"""
is_supported_fft_length(n::Integer) = n >= 1 && _odd_cofactor(n) in (1, 3, 5, 15)

# Throw an informative error for an unsupported 1D FFT length.
@noinline function _unsupported_fft_length(n::Integer)
    throw(ArgumentError(string(
        "unsupported FFT length $n: AppleAccelerate only supports 1D complex FFT ",
        "lengths of the form f*2^k with f ∈ {1, 3, 5, 15}. ",
        "For arbitrary (e.g. prime) lengths, use FFTW.jl instead.")))
end

# --- Internal 1D FFT (direction-based) ---

function _fft1d(r::Vector{ComplexF64}, setup::FFTSetup{Float64}, direction::Int)
    n = length(r)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)

    realp = real(r)
    imagp = imag(r)
    retr = Vector{Float64}(undef, n)
    reti = Vector{Float64}(undef, n)

    GC.@preserve realp imagp retr reti begin
        input = DSPDoubleSplitComplex(realp, imagp)
        output = DSPDoubleSplitComplex(retr, reti)
        LibAccelerate.vDSP_fft_zopD(setup.plan, Ref(input), SIGNAL_STRIDE,
                                    Ref(output), SIGNAL_STRIDE, logn, direction)
    end

    return complex.(retr, reti)
end

function _fft1d(r::Vector{ComplexF32}, setup::FFTSetup{Float32}, direction::Int)
    n = length(r)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)

    realp = Float32.(real(r))
    imagp = Float32.(imag(r))
    retr = Vector{Float32}(undef, n)
    reti = Vector{Float32}(undef, n)

    GC.@preserve realp imagp retr reti begin
        input = DSPSplitComplex(realp, imagp)
        output = DSPSplitComplex(retr, reti)
        LibAccelerate.vDSP_fft_zop(setup.plan, Ref(input), SIGNAL_STRIDE,
                                   Ref(output), SIGNAL_STRIDE, logn, direction)
    end

    return complex.(retr, reti)
end

# --- Internal 2D FFT (direction-based) ---

function _fft2d(r::Matrix{ComplexF64}, setup::FFTSetup{Float64}, direction::Int)
    nrows, ncols = size(r)
    @assert ispow2(nrows) && ispow2(ncols) "dimensions must be powers of 2"
    log2nr = trailing_zeros(nrows)
    log2nc = trailing_zeros(ncols)

    realp = real.(r)
    imagp = imag.(r)
    retr = Matrix{Float64}(undef, nrows, ncols)
    reti = Matrix{Float64}(undef, nrows, ncols)

    GC.@preserve realp imagp retr reti begin
        input = DSPDoubleSplitComplex(pointer(realp), pointer(imagp))
        output = DSPDoubleSplitComplex(pointer(retr), pointer(reti))
        LibAccelerate.vDSP_fft2d_zopD(setup.plan, Ref(input), SIGNAL_STRIDE, 0,
                                      Ref(output), SIGNAL_STRIDE, 0,
                                      log2nr, log2nc, direction)
    end

    return complex.(retr, reti)
end

function _fft2d(r::Matrix{ComplexF32}, setup::FFTSetup{Float32}, direction::Int)
    nrows, ncols = size(r)
    @assert ispow2(nrows) && ispow2(ncols) "dimensions must be powers of 2"
    log2nr = trailing_zeros(nrows)
    log2nc = trailing_zeros(ncols)

    realp = Float32.(real.(r))
    imagp = Float32.(imag.(r))
    retr = Matrix{Float32}(undef, nrows, ncols)
    reti = Matrix{Float32}(undef, nrows, ncols)

    GC.@preserve realp imagp retr reti begin
        input = DSPSplitComplex(pointer(realp), pointer(imagp))
        output = DSPSplitComplex(pointer(retr), pointer(reti))
        LibAccelerate.vDSP_fft2d_zop(setup.plan, Ref(input), SIGNAL_STRIDE, 0,
                                     Ref(output), SIGNAL_STRIDE, 0,
                                     log2nr, log2nc, direction)
    end

    return complex.(retr, reti)
end

# --- Public API: fft (forward FFT) ---

# Mixed-radix 1D complex transform via the vDSP DFT (non-power-of-2). `direction` is
# DFT_FORWARD or DFT_INVERSE. Result is unnormalized (matching bfft/dft conventions).
function _fft1d_dft(x::Vector{Complex{T}}, direction::Int) where {T<:Union{Float32,Float64}}
    setup = _cached_dftsetup(T, length(x), direction)  # shared cached setup
    return dft(x, setup)
end

"""
    fft(x::VecOrMat{Complex{T}}, [setup::FFTSetup{T}])

Compute the forward FFT of `x` via Apple vDSP. Supports 1D vectors and 2D matrices
with `ComplexF32` or `ComplexF64` elements.

For a 1D vector with no explicit `setup`, non-power-of-2 lengths of the form
`f * 2^k` with `f ∈ {1, 3, 5, 15}` are also supported and transparently use Apple's
mixed-radix DFT (see [`is_supported_fft_length`](@ref)); an unsupported length throws
an `ArgumentError`. When a `setup::FFTSetup` is supplied, or for 2D inputs, all
dimensions must be powers of 2. If `setup` is omitted, a temporary plan is created
automatically.
Wraps [`vDSP_fft_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zop) (1D) /
[`vDSP_fft2d_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zop) (2D).
"""
fft(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d(x, setup, FFT_FORWARD)
fft(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d(x, setup, FFT_FORWARD)
function fft(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}}
    n = length(x)
    ispow2(n) && return fft(x, _cached_fftsetup(x))
    is_supported_fft_length(n) || _unsupported_fft_length(n)
    return _fft1d_dft(x, DFT_FORWARD)
end
fft(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = fft(x, _cached_fftsetup(x))

# --- Public API: bfft (backward/unnormalized inverse FFT) ---

"""
    bfft(x::VecOrMat{Complex{T}}, [setup::FFTSetup{T}])

Compute the unnormalized inverse (backward) FFT of `x` via Apple vDSP.
The result is *not* divided by `length(x)`; use [`ifft`](@ref) for the normalized version.

For a 1D vector with no explicit `setup`, non-power-of-2 lengths of the form
`f * 2^k` with `f ∈ {1, 3, 5, 15}` are also supported (via Apple's mixed-radix DFT);
an unsupported length throws an `ArgumentError`. With an explicit `setup::FFTSetup`,
or for 2D inputs, all dimensions must be powers of 2.
Wraps [`vDSP_fft_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zop) (1D) /
[`vDSP_fft2d_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zop) (2D).
"""
bfft(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d(x, setup, FFT_INVERSE)
bfft(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d(x, setup, FFT_INVERSE)
function bfft(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}}
    n = length(x)
    ispow2(n) && return bfft(x, _cached_fftsetup(x))
    is_supported_fft_length(n) || _unsupported_fft_length(n)
    return _fft1d_dft(x, DFT_INVERSE)
end
bfft(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft(x, _cached_fftsetup(x))

# --- Public API: ifft (normalized inverse FFT) ---

"""
    ifft(x::VecOrMat{Complex{T}}, [setup::FFTSetup{T}])

Compute the normalized inverse FFT of `x` via Apple vDSP.
Equivalent to `bfft(x) / length(x)`. Satisfies `ifft(fft(x)) ≈ x`.

For a 1D vector with no explicit `setup`, non-power-of-2 lengths of the form
`f * 2^k` with `f ∈ {1, 3, 5, 15}` are also supported (via Apple's mixed-radix DFT);
an unsupported length throws an `ArgumentError`. With an explicit `setup::FFTSetup`,
or for 2D inputs, all dimensions must be powers of 2.
Wraps [`vDSP_fft_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zop) (1D) /
[`vDSP_fft2d_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zop) (2D).
"""
ifft(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = bfft(x, setup) ./ length(x)
ifft(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = bfft(x, setup) ./ length(x)
# No-setup 1D routes through bfft, which handles both power-of-2 and mixed-radix lengths.
ifft(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft(x) ./ length(x)
ifft(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft(x, _cached_fftsetup(x))

# --- Internal in-place 1D complex FFT ---

function _fft1d!(x::Vector{ComplexF64}, setup::FFTSetup{Float64}, direction::Int)
    n = length(x)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)

    realp = real.(x)
    imagp = imag.(x)

    GC.@preserve realp imagp begin
        c = DSPDoubleSplitComplex(realp, imagp)
        LibAccelerate.vDSP_fft_zipD(setup.plan, Ref(c), SIGNAL_STRIDE, logn, direction)
    end

    @inbounds for i in eachindex(x)
        x[i] = complex(realp[i], imagp[i])
    end
    return x
end

function _fft1d!(x::Vector{ComplexF32}, setup::FFTSetup{Float32}, direction::Int)
    n = length(x)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)

    realp = Float32.(real.(x))
    imagp = Float32.(imag.(x))

    GC.@preserve realp imagp begin
        c = DSPSplitComplex(realp, imagp)
        LibAccelerate.vDSP_fft_zip(setup.plan, Ref(c), SIGNAL_STRIDE, logn, direction)
    end

    @inbounds for i in eachindex(x)
        x[i] = complex(realp[i], imagp[i])
    end
    return x
end

# --- Internal in-place 2D complex FFT ---

function _fft2d!(x::Matrix{ComplexF64}, setup::FFTSetup{Float64}, direction::Int)
    nrows, ncols = size(x)
    @assert ispow2(nrows) && ispow2(ncols) "dimensions must be powers of 2"
    log2nr = trailing_zeros(nrows)
    log2nc = trailing_zeros(ncols)

    realp = real.(x)
    imagp = imag.(x)

    GC.@preserve realp imagp begin
        c = DSPDoubleSplitComplex(pointer(realp), pointer(imagp))
        LibAccelerate.vDSP_fft2d_zipD(setup.plan, Ref(c), SIGNAL_STRIDE, 0,
                                      log2nr, log2nc, direction)
    end

    @inbounds for i in eachindex(x)
        x[i] = complex(realp[i], imagp[i])
    end
    return x
end

function _fft2d!(x::Matrix{ComplexF32}, setup::FFTSetup{Float32}, direction::Int)
    nrows, ncols = size(x)
    @assert ispow2(nrows) && ispow2(ncols) "dimensions must be powers of 2"
    log2nr = trailing_zeros(nrows)
    log2nc = trailing_zeros(ncols)

    realp = Float32.(real.(x))
    imagp = Float32.(imag.(x))

    GC.@preserve realp imagp begin
        c = DSPSplitComplex(pointer(realp), pointer(imagp))
        LibAccelerate.vDSP_fft2d_zip(setup.plan, Ref(c), SIGNAL_STRIDE, 0,
                                     log2nr, log2nc, direction)
    end

    @inbounds for i in eachindex(x)
        x[i] = complex(realp[i], imagp[i])
    end
    return x
end

# --- Public API: fft! (in-place forward FFT) ---

"""
    fft!(x::VecOrMat{Complex{T}}, [setup::FFTSetup{T}])

Compute the forward FFT of `x` in-place via Apple vDSP, overwriting `x` with the result.
Wraps [`vDSP_fft_zip`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zip) (1D) /
[`vDSP_fft2d_zip`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zip) (2D).
"""
fft!(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d!(x, setup, FFT_FORWARD)
fft!(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d!(x, setup, FFT_FORWARD)
fft!(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = fft!(x, _cached_fftsetup(x))
fft!(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = fft!(x, _cached_fftsetup(x))

# --- Public API: bfft! (in-place backward/unnormalized inverse FFT) ---

"""
    bfft!(x::VecOrMat{Complex{T}}, [setup::FFTSetup{T}])

Compute the unnormalized inverse FFT of `x` in-place via Apple vDSP.
Wraps [`vDSP_fft_zip`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zip) (1D) /
[`vDSP_fft2d_zip`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zip) (2D).
"""
bfft!(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d!(x, setup, FFT_INVERSE)
bfft!(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d!(x, setup, FFT_INVERSE)
bfft!(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft!(x, _cached_fftsetup(x))
bfft!(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft!(x, _cached_fftsetup(x))

# --- Public API: ifft! (in-place normalized inverse FFT) ---

"""
    ifft!(x::VecOrMat{Complex{T}}, [setup::FFTSetup{T}])

Compute the normalized inverse FFT of `x` in-place via Apple vDSP.
Satisfies `ifft!(fft!(copy(x))) ≈ x`.
Wraps [`vDSP_fft_zip`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zip) (1D) /
[`vDSP_fft2d_zip`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zip) (2D).
"""
function ifft!(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}}
    bfft!(x, setup)
    x ./= length(x)
end
function ifft!(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}}
    bfft!(x, setup)
    x ./= length(x)
end
ifft!(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft!(x, _cached_fftsetup(x))
ifft!(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft!(x, _cached_fftsetup(x))

# --- Internal 1D real FFT ---
# vDSP real FFT uses packed split-complex format:
#   Forward output: realp[0]=DC, imagp[0]=Nyquist, realp[k]+i*imagp[k] for k=1..N/2-1
#   Forward scale: output is 2x the mathematical DFT
# We unpack to standard format: Complex vector of length N/2+1

function _rfft1d(x::Vector{Float64}, setup::FFTSetup{Float64})
    n = length(x)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)
    half = n >> 1

    # Pack real input into split-complex: realp[k]=x[2k-1], imagp[k]=x[2k] (1-based)
    inp_realp = x[1:2:end]
    inp_imagp = x[2:2:end]
    out_realp = Vector{Float64}(undef, half)
    out_imagp = Vector{Float64}(undef, half)

    GC.@preserve inp_realp inp_imagp out_realp out_imagp begin
        input = DSPDoubleSplitComplex(inp_realp, inp_imagp)
        output = DSPDoubleSplitComplex(out_realp, out_imagp)
        LibAccelerate.vDSP_fft_zropD(setup.plan, Ref(input), SIGNAL_STRIDE,
                                     Ref(output), SIGNAL_STRIDE, logn, FFT_FORWARD)
    end

    # Unpack and divide by 2 (vDSP scales forward real FFT by 2)
    result = Vector{ComplexF64}(undef, half + 1)
    result[1] = complex(out_realp[1] / 2, 0.0)
    @inbounds for k in 2:half
        result[k] = complex(out_realp[k] / 2, out_imagp[k] / 2)
    end
    result[half + 1] = complex(out_imagp[1] / 2, 0.0)
    return result
end

function _rfft1d(x::Vector{Float32}, setup::FFTSetup{Float32})
    n = length(x)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)
    half = n >> 1

    inp_realp = Float32.(x[1:2:end])
    inp_imagp = Float32.(x[2:2:end])
    out_realp = Vector{Float32}(undef, half)
    out_imagp = Vector{Float32}(undef, half)

    GC.@preserve inp_realp inp_imagp out_realp out_imagp begin
        input = DSPSplitComplex(inp_realp, inp_imagp)
        output = DSPSplitComplex(out_realp, out_imagp)
        LibAccelerate.vDSP_fft_zrop(setup.plan, Ref(input), SIGNAL_STRIDE,
                                    Ref(output), SIGNAL_STRIDE, logn, FFT_FORWARD)
    end

    result = Vector{ComplexF32}(undef, half + 1)
    result[1] = complex(out_realp[1] / 2, Float32(0))
    @inbounds for k in 2:half
        result[k] = complex(out_realp[k] / 2, out_imagp[k] / 2)
    end
    result[half + 1] = complex(out_imagp[1] / 2, Float32(0))
    return result
end

# --- Internal 1D inverse real FFT (unnormalized) ---
# Pack standard complex input back into vDSP format, run inverse, unpack real output.

function _brfft1d(X::Vector{ComplexF64}, n::Int, setup::FFTSetup{Float64})
    @assert ispow2(n) "output length must be a power of 2"
    logn = trailing_zeros(n)
    half = n >> 1
    @assert length(X) == half + 1 "input must have length n÷2+1"

    # Pack standard complex input into vDSP format
    inp_realp = Vector{Float64}(undef, half)
    inp_imagp = Vector{Float64}(undef, half)
    inp_realp[1] = real(X[1])
    inp_imagp[1] = real(X[half + 1])
    @inbounds for k in 2:half
        inp_realp[k] = real(X[k])
        inp_imagp[k] = imag(X[k])
    end

    out_realp = Vector{Float64}(undef, half)
    out_imagp = Vector{Float64}(undef, half)

    GC.@preserve inp_realp inp_imagp out_realp out_imagp begin
        input = DSPDoubleSplitComplex(inp_realp, inp_imagp)
        output = DSPDoubleSplitComplex(out_realp, out_imagp)
        LibAccelerate.vDSP_fft_zropD(setup.plan, Ref(input), SIGNAL_STRIDE,
                                     Ref(output), SIGNAL_STRIDE, logn, FFT_INVERSE)
    end

    # Unpack interleaved real output
    result = Vector{Float64}(undef, n)
    @inbounds for k in 1:half
        result[2k - 1] = out_realp[k]
        result[2k] = out_imagp[k]
    end
    return result
end

function _brfft1d(X::Vector{ComplexF32}, n::Int, setup::FFTSetup{Float32})
    @assert ispow2(n) "output length must be a power of 2"
    logn = trailing_zeros(n)
    half = n >> 1
    @assert length(X) == half + 1 "input must have length n÷2+1"

    inp_realp = Vector{Float32}(undef, half)
    inp_imagp = Vector{Float32}(undef, half)
    inp_realp[1] = real(X[1])
    inp_imagp[1] = real(X[half + 1])
    @inbounds for k in 2:half
        inp_realp[k] = real(X[k])
        inp_imagp[k] = imag(X[k])
    end

    out_realp = Vector{Float32}(undef, half)
    out_imagp = Vector{Float32}(undef, half)

    GC.@preserve inp_realp inp_imagp out_realp out_imagp begin
        input = DSPSplitComplex(inp_realp, inp_imagp)
        output = DSPSplitComplex(out_realp, out_imagp)
        LibAccelerate.vDSP_fft_zrop(setup.plan, Ref(input), SIGNAL_STRIDE,
                                    Ref(output), SIGNAL_STRIDE, logn, FFT_INVERSE)
    end

    result = Vector{Float32}(undef, n)
    @inbounds for k in 1:half
        result[2k - 1] = out_realp[k]
        result[2k] = out_imagp[k]
    end
    return result
end

# --- Internal 2D real FFT ---
#
# vDSP 2D real FFT (vDSP_fft2d_zrop) conventions, mapped to Julia's column-major
# layout: Julia's first (contiguous) dimension is vDSP's N0 and the second is N1.
# For an n1×n2 real input the transform works on an (n1/2)×n2 split-complex array.
#
# Input packing (same even/odd interleave as the 1D real FFT, along dim 1):
#   realp[k, j] = x[2k-1, j],  imagp[k, j] = x[2k, j]
#
# Output packing (from Apple's vDSP_fft2d_zrip documentation, translated so that
# H[k0, k1] is the full 2D DFT with k0 the dim-1 frequency and k1 the dim-2
# frequency, both 0-based; h1 = n1/2, h2 = n2/2):
#   * For 0 < k0 < h1 and all k1 (the "regular" region):
#       C[k0+1, k1+1] = H[k0, k1]
#   * Row 1 of C packs the two real-Hermitian columns k0 = 0 (in realp) and
#     k0 = h1 (in imagp), each in 1D real-FFT packed order with Re/Im
#     interleaved along dim 2 (verified empirically against FFTW):
#       realp[1, 1] = H[0, 0]  (DC, real)      imagp[1, 1] = H[h1, 0]   (real)
#       realp[1, 2] = H[0, h2] (Nyquist, real) imagp[1, 2] = H[h1, h2]  (real)
#       for 0 < k1 < h2:
#         realp[1, 2k1+1] = Re(H[0, k1])    imagp[1, 2k1+1] = Re(H[h1, k1])
#         realp[1, 2k1+2] = Im(H[0, k1])    imagp[1, 2k1+2] = Im(H[h1, k1])
#     The remaining H[0, k1] / H[h1, k1] for k1 > h2 follow from Hermitian
#     symmetry: H[0, n2-k1] = conj(H[0, k1]), H[h1, n2-k1] = conj(H[h1, k1]).
# Like the 1D real FFT, the forward output is scaled by 2 relative to the
# mathematical DFT, and the inverse of the (unscaled) DFT coefficients returns
# n1*n2 times the original signal.
#
# The unpacked result matches FFTW's `rfft(A)`: an (n1÷2+1)×n2 complex matrix.

for (T, SC, fft2d_zrop) in ((Float64, :DSPDoubleSplitComplex, :vDSP_fft2d_zropD),
                            (Float32, :DSPSplitComplex, :vDSP_fft2d_zrop))
    @eval begin
        function _rfft2d(x::Matrix{$T}, setup::FFTSetup{$T})
            n1, n2 = size(x)
            @assert ispow2(n1) && ispow2(n2) "dimensions must be powers of 2"
            @assert n1 >= 2 && n2 >= 2 "each dimension must be at least 2"
            log2n1 = trailing_zeros(n1)
            log2n2 = trailing_zeros(n2)
            h1 = n1 >> 1
            h2 = n2 >> 1

            # Even/odd split of the contiguous (first) dimension.
            inp_realp = x[1:2:end, :]
            inp_imagp = x[2:2:end, :]
            out_realp = Matrix{$T}(undef, h1, n2)
            out_imagp = Matrix{$T}(undef, h1, n2)

            GC.@preserve inp_realp inp_imagp out_realp out_imagp begin
                input = $SC(pointer(inp_realp), pointer(inp_imagp))
                output = $SC(pointer(out_realp), pointer(out_imagp))
                LibAccelerate.$fft2d_zrop(setup.plan, Ref(input), SIGNAL_STRIDE, 0,
                                          Ref(output), SIGNAL_STRIDE, 0,
                                          log2n1, log2n2, FFT_FORWARD)
            end

            # Unpack into FFTW `rfft` layout, dividing by 2 (vDSP forward scaling).
            result = Matrix{Complex{$T}}(undef, h1 + 1, n2)
            @inbounds for j in 1:n2, k in 2:h1
                result[k, j] = complex(out_realp[k, j], out_imagp[k, j]) / 2
            end
            # Purely real corner points (DC/Nyquist combinations).
            result[1, 1]        = complex(out_realp[1, 1] / 2)
            result[h1+1, 1]     = complex(out_imagp[1, 1] / 2)
            result[1, h2+1]     = complex(out_realp[1, 2] / 2)
            result[h1+1, h2+1]  = complex(out_imagp[1, 2] / 2)
            # Packed k0 = 0 and k0 = n1/2 rows; upper half via Hermitian symmetry.
            @inbounds for k1 in 1:h2-1
                result[1, k1+1]       = complex(out_realp[1, 2k1+1], out_realp[1, 2k1+2]) / 2
                result[h1+1, k1+1]    = complex(out_imagp[1, 2k1+1], out_imagp[1, 2k1+2]) / 2
                result[1, n2-k1+1]    = conj(result[1, k1+1])
                result[h1+1, n2-k1+1] = conj(result[h1+1, k1+1])
            end
            return result
        end

        function _brfft2d(X::Matrix{Complex{$T}}, n1::Int, setup::FFTSetup{$T})
            n2 = size(X, 2)
            @assert ispow2(n1) && ispow2(n2) "dimensions must be powers of 2"
            @assert n1 >= 2 && n2 >= 2 "each dimension must be at least 2"
            log2n1 = trailing_zeros(n1)
            log2n2 = trailing_zeros(n2)
            h1 = n1 >> 1
            h2 = n2 >> 1
            @assert size(X, 1) == h1 + 1 "input must have size (n1÷2+1)×n2"

            # Pack the FFTW-layout coefficients back into the vDSP packed format
            # (inverse of the unpacking above; the redundant upper halves of the
            # k0 = 0 and k0 = n1/2 rows are dropped).
            inp_realp = Matrix{$T}(undef, h1, n2)
            inp_imagp = Matrix{$T}(undef, h1, n2)
            @inbounds for j in 1:n2, k in 2:h1
                inp_realp[k, j] = real(X[k, j])
                inp_imagp[k, j] = imag(X[k, j])
            end
            inp_realp[1, 1] = real(X[1, 1])
            inp_imagp[1, 1] = real(X[h1+1, 1])
            inp_realp[1, 2] = real(X[1, h2+1])
            inp_imagp[1, 2] = real(X[h1+1, h2+1])
            @inbounds for k1 in 1:h2-1
                inp_realp[1, 2k1+1] = real(X[1, k1+1])
                inp_realp[1, 2k1+2] = imag(X[1, k1+1])
                inp_imagp[1, 2k1+1] = real(X[h1+1, k1+1])
                inp_imagp[1, 2k1+2] = imag(X[h1+1, k1+1])
            end

            out_realp = Matrix{$T}(undef, h1, n2)
            out_imagp = Matrix{$T}(undef, h1, n2)

            GC.@preserve inp_realp inp_imagp out_realp out_imagp begin
                input = $SC(pointer(inp_realp), pointer(inp_imagp))
                output = $SC(pointer(out_realp), pointer(out_imagp))
                LibAccelerate.$fft2d_zrop(setup.plan, Ref(input), SIGNAL_STRIDE, 0,
                                          Ref(output), SIGNAL_STRIDE, 0,
                                          log2n1, log2n2, FFT_INVERSE)
            end

            # Unpack even/odd interleaved real output.
            result = Matrix{$T}(undef, n1, n2)
            @inbounds for j in 1:n2, k in 1:h1
                result[2k - 1, j] = out_realp[k, j]
                result[2k, j] = out_imagp[k, j]
            end
            return result
        end
    end
end

# --- Internal 1D real DFT (mixed-radix, non-power-of-2 lengths) ---
#
# vDSP_DFT_zrop uses the same even/odd input packing and DC/Nyquist output
# packing as vDSP_fft_zrop, including the forward scaling by 2, but supports
# lengths of the form f * 2^k with f ∈ {1, 3, 5, 15}. Setup creation returns
# NULL for unsupported lengths.

# Throw an informative error for an unsupported 1D real FFT length.
@noinline function _unsupported_rfft_length(n::Integer)
    throw(ArgumentError(string(
        "unsupported real FFT length $n: AppleAccelerate supports power-of-2 real FFT ",
        "lengths, plus mixed-radix lengths of the form f*2^k with f ∈ {3, 5, 15} ",
        "(k ≥ 4) via Apple's real-input DFT. For other lengths, use FFTW.jl instead.")))
end

# Real-input DFT setups are cached like the complex ones (see the global setup
# caches above), in a separate dict since zrop setups are distinct from zop ones.
const _RDFT_SETUP_CACHE = Dict{Tuple{DataType,Int,Int},DFTSetup}()

# Shared, cached real-input DFTSetup{T}, keyed by length and direction.
function _cached_rdftsetup(::Type{T}, n::Int, direction::Int) where {T<:Union{Float32,Float64}}
    key = (T, n, direction)
    return lock(_SETUP_CACHE_LOCK) do
        get!(() -> _plan_rdft(n, direction, T), _RDFT_SETUP_CACHE, key)
    end::DFTSetup{T}
end

for (T, createfn, execfn) in ((Float64, :vDSP_DFT_zrop_CreateSetupD, :vDSP_DFT_ExecuteD),
                              (Float32, :vDSP_DFT_zrop_CreateSetup, :vDSP_DFT_Execute))
    @eval begin
        # Create a real-input DFT setup; throws for unsupported lengths.
        function _plan_rdft(n::Int, direction::Int, ::Type{$T})
            setup = LibAccelerate.$createfn(C_NULL, n, Cint(direction))
            setup == C_NULL && _unsupported_rfft_length(n)
            return DFTSetup($T, Ptr{Cvoid}(setup), direction)
        end

        function _rfft1d_dft(x::Vector{$T})
            n = length(x)
            half = n >> 1
            setup = _cached_rdftsetup($T, n, DFT_FORWARD)
            Ir = x[1:2:end]
            Ii = x[2:2:end]
            Or = Vector{$T}(undef, half)
            Oi = Vector{$T}(undef, half)
            GC.@preserve setup begin
                LibAccelerate.$execfn(setup.setup, Ir, Ii, Or, Oi)
            end
            result = Vector{Complex{$T}}(undef, half + 1)
            result[1] = complex(Or[1] / 2)
            @inbounds for k in 2:half
                result[k] = complex(Or[k], Oi[k]) / 2
            end
            result[half + 1] = complex(Oi[1] / 2)
            return result
        end

        function _brfft1d_dft(X::Vector{Complex{$T}}, n::Int)
            half = n >> 1
            length(X) == half + 1 || throw(DimensionMismatch(
                "brfft: input must have length n÷2+1 = $(half + 1); got $(length(X))"))
            setup = _cached_rdftsetup($T, n, DFT_INVERSE)
            Ir = Vector{$T}(undef, half)
            Ii = Vector{$T}(undef, half)
            Ir[1] = real(X[1])
            Ii[1] = real(X[half + 1])
            @inbounds for k in 2:half
                Ir[k] = real(X[k])
                Ii[k] = imag(X[k])
            end
            Or = Vector{$T}(undef, half)
            Oi = Vector{$T}(undef, half)
            GC.@preserve setup begin
                LibAccelerate.$execfn(setup.setup, Ir, Ii, Or, Oi)
            end
            result = Vector{$T}(undef, n)
            @inbounds for k in 1:half
                result[2k - 1] = Or[k]
                result[2k] = Oi[k]
            end
            return result
        end
    end
end

# --- Public API: rfft (forward real FFT) ---

"""
    plan_rfft(x::VecOrMat{T}) where T <: Union{Float32, Float64}

Create a reusable FFT setup for real-input forward transforms of the same size as `x`.
For a matrix, the setup covers both dimensions (and can also be reused for the
matching inverse transforms).
Wraps [`vDSP_create_fftsetup`](https://developer.apple.com/documentation/accelerate/vdsp_create_fftsetup).
"""
plan_rfft(x::Vector{T}) where {T<:Union{Float32,Float64}} = FFTSetup{T}(length(x))
plan_rfft(x::Matrix{T}) where {T<:Union{Float32,Float64}} = FFTSetup{T}(max(size(x)...))

"""
    rfft(x::VecOrMat{T}, [setup::FFTSetup{T}])

Compute the forward FFT of a real-valued vector or matrix `x` via Apple vDSP.
Returns the non-redundant complex coefficients: length `N÷2+1` for a vector of
length `N`, and size `(n1÷2+1)×n2` for an `n1×n2` matrix (matching FFTW's `rfft`).

For a 1D vector with no explicit `setup`, non-power-of-2 lengths supported by
Apple's real-input mixed-radix DFT (`f * 2^k`, `f ∈ {3, 5, 15}`, `k ≥ 4`) are also
accepted; an unsupported length throws an `ArgumentError`. With an explicit
`setup::FFTSetup`, or for 2D inputs, all dimensions must be powers of 2.
Wraps [`vDSP_fft_zrop`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zrop) (1D) /
[`vDSP_fft2d_zrop`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zrop) (2D) /
[`vDSP_DFT_zrop_CreateSetup`](https://developer.apple.com/documentation/accelerate/vdsp_dft_zrop_createsetup) (1D mixed-radix).
"""
rfft(x::Vector{T}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _rfft1d(x, setup)
rfft(x::Matrix{T}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _rfft2d(x, setup)
function rfft(x::Vector{T}) where {T<:Union{Float32,Float64}}
    n = length(x)
    ispow2(n) && return rfft(x, _cached_fftsetup(T, n))
    (iseven(n) && is_supported_fft_length(n)) || _unsupported_rfft_length(n)
    return _rfft1d_dft(x)
end
rfft(x::Matrix{T}) where {T<:Union{Float32,Float64}} = rfft(x, _cached_fftsetup(T, max(size(x)...)))

# --- Public API: brfft (backward/unnormalized inverse real FFT) ---

"""
    brfft(X::Vector{Complex{T}}, n::Int, [setup::FFTSetup{T}])
    brfft(X::Matrix{Complex{T}}, n1::Int, [setup::FFTSetup{T}])

Compute the unnormalized inverse real FFT, returning a real vector of length `n`
(`X` must have length `n÷2+1`) or a real matrix of size `n1×n2` where
`n2 = size(X, 2)` (`X` must have size `(n1÷2+1)×n2`). The result is *not* divided
by the output length; use [`irfft`](@ref) for the normalized version.

For a 1D input with no explicit `setup`, non-power-of-2 lengths supported by
Apple's real-input mixed-radix DFT are also accepted (see [`rfft`](@ref)).
Wraps [`vDSP_fft_zrop`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zrop) (1D) /
[`vDSP_fft2d_zrop`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zrop) (2D).
"""
brfft(X::Vector{Complex{T}}, n::Int, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _brfft1d(X, n, setup)
brfft(X::Matrix{Complex{T}}, n1::Int, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _brfft2d(X, n1, setup)
function brfft(X::Vector{Complex{T}}, n::Int) where {T<:Union{Float32,Float64}}
    ispow2(n) && return brfft(X, n, _cached_fftsetup(T, n))
    (iseven(n) && is_supported_fft_length(n)) || _unsupported_rfft_length(n)
    return _brfft1d_dft(X, n)
end
brfft(X::Matrix{Complex{T}}, n1::Int) where {T<:Union{Float32,Float64}} =
    brfft(X, n1, _cached_fftsetup(T, max(n1, size(X, 2))))

# --- Public API: irfft (normalized inverse real FFT) ---

"""
    irfft(X::Vector{Complex{T}}, n::Int, [setup::FFTSetup{T}])
    irfft(X::Matrix{Complex{T}}, n1::Int, [setup::FFTSetup{T}])

Compute the normalized inverse real FFT, returning a real vector of length `n`
or a real matrix of size `n1×n2` where `n2 = size(X, 2)`.
Satisfies `irfft(rfft(x), size(x, 1)) ≈ x`.

For a 1D input with no explicit `setup`, non-power-of-2 lengths supported by
Apple's real-input mixed-radix DFT are also accepted (see [`rfft`](@ref)).
Wraps [`vDSP_fft_zrop`](https://developer.apple.com/documentation/accelerate/vdsp_fft_zrop) (1D) /
[`vDSP_fft2d_zrop`](https://developer.apple.com/documentation/accelerate/vdsp_fft2d_zrop) (2D).
"""
irfft(X::Vector{Complex{T}}, n::Int, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = brfft(X, n, setup) ./ n
irfft(X::Vector{Complex{T}}, n::Int) where {T<:Union{Float32,Float64}} = brfft(X, n) ./ n
irfft(X::Matrix{Complex{T}}, n1::Int, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} =
    brfft(X, n1, setup) ./ (n1 * size(X, 2))
irfft(X::Matrix{Complex{T}}, n1::Int) where {T<:Union{Float32,Float64}} =
    brfft(X, n1) ./ (n1 * size(X, 2))

# --- Batched 1D complex FFT (vDSP_fftm_*) ---

for (T, SC, fftm_zop) in ((Float64, :DSPDoubleSplitComplex, :vDSP_fftm_zopD),
                          (Float32, :DSPSplitComplex, :vDSP_fftm_zop))
    @eval begin
        function _fftm(x::Matrix{Complex{$T}}, dims::Integer, setup::FFTSetup{$T}, direction::Int)
            dims == 1 || dims == 2 || throw(ArgumentError("dims must be 1 or 2; got $dims"))
            nrows, ncols = size(x)
            n = size(x, dims)                      # transform length
            m = dims == 1 ? ncols : nrows          # number of transforms
            @assert ispow2(n) "transform length must be a power of 2"
            logn = trailing_zeros(n)
            # Julia is column-major: for column transforms (dims=1) the element
            # stride is 1 and consecutive transforms are nrows apart; for row
            # transforms (dims=2) it is the other way around.
            elstride, batchstride = dims == 1 ? (1, nrows) : (nrows, 1)

            realp = $T.(real.(x))
            imagp = $T.(imag.(x))
            retr = Matrix{$T}(undef, nrows, ncols)
            reti = Matrix{$T}(undef, nrows, ncols)

            GC.@preserve realp imagp retr reti begin
                input = $SC(pointer(realp), pointer(imagp))
                output = $SC(pointer(retr), pointer(reti))
                LibAccelerate.$fftm_zop(setup.plan, Ref(input), elstride, batchstride,
                                        Ref(output), elstride, batchstride,
                                        logn, m, direction)
            end

            return complex.(retr, reti)
        end
    end
end

"""
    fft(x::Matrix{Complex{T}}, dims::Integer, [setup::FFTSetup{T}])

Batched 1D forward FFT: transform each column (`dims = 1`) or each row
(`dims = 2`) of `x` independently, like FFTW's `fft(x, dims)`. The transform
length (`size(x, dims)`) must be a power of 2.
Wraps [`vDSP_fftm_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fftm_zop).
"""
fft(x::Matrix{Complex{T}}, dims::Integer, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} =
    _fftm(x, dims, setup, FFT_FORWARD)
fft(x::Matrix{Complex{T}}, dims::Integer) where {T<:Union{Float32,Float64}} =
    fft(x, dims, _cached_fftsetup(T, size(x, dims)))

"""
    bfft(x::Matrix{Complex{T}}, dims::Integer, [setup::FFTSetup{T}])

Batched 1D unnormalized inverse (backward) FFT along dimension `dims`
(1 = columns, 2 = rows). Use [`ifft`](@ref) for the normalized version.
Wraps [`vDSP_fftm_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fftm_zop).
"""
bfft(x::Matrix{Complex{T}}, dims::Integer, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} =
    _fftm(x, dims, setup, FFT_INVERSE)
bfft(x::Matrix{Complex{T}}, dims::Integer) where {T<:Union{Float32,Float64}} =
    bfft(x, dims, _cached_fftsetup(T, size(x, dims)))

"""
    ifft(x::Matrix{Complex{T}}, dims::Integer, [setup::FFTSetup{T}])

Batched 1D normalized inverse FFT along dimension `dims` (1 = columns, 2 = rows).
Satisfies `ifft(fft(x, dims), dims) ≈ x`.
Wraps [`vDSP_fftm_zop`](https://developer.apple.com/documentation/accelerate/vdsp_fftm_zop).
"""
ifft(x::Matrix{Complex{T}}, dims::Integer, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} =
    bfft(x, dims, setup) ./ size(x, dims)
ifft(x::Matrix{Complex{T}}, dims::Integer) where {T<:Union{Float32,Float64}} =
    bfft(x, dims) ./ size(x, dims)
