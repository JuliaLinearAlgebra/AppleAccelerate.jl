## DSP.jl ##

mutable struct DFTSetup{T}
    setup::Ptr{Cvoid}
    direction::Int

    function DFTSetup{T}(setup::Ptr{Cvoid}, direction::Int) where T
        dftsetup = new(setup, direction)
        finalizer(plan_destroy, dftsetup)
        dftsetup
    end
end

mutable struct Biquad{T}
    setup::Ptr{Cvoid}
    sections::Int

    function Biquad{T}(setup::Ptr{Cvoid}, sections::Int) where T
        biquadsetup = new(setup, sections)
        finalizer(biquaddestroy, biquadsetup)
        biquadsetup
    end
end

## === FUNCTIONS == ##

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
            xpadded::Vector{$T} = [zeros($T, ksize-1); X; zeros($T, ksize)]
            ccall(($(string("vDSP_conv", suff), libacc)),  Cvoid,
                  (Ptr{$T}, Int64,  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64, UInt64),
                  xpadded, 1, pointer(K, ksize), -1, result, 1,  rsize, ksize)
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
            xpadded::Vector{$T} = [zeros($T, ysize-1); X; zeros($T, ysize)]
            ccall(($(string("vDSP_conv", suff), libacc)),  Cvoid,
                  (Ptr{$T}, Int64,  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64, UInt64),
                  xpadded, 1, Y, 1, result, 1,  rsize, ysize)
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

## == Biquadratic/IIR filtering
for (T, suff) in ((Float64, "D"), )

    """
    Initializes a vDSP_biquad_setup for use with vDSP_biquad. A multi-section filter
    can be initialized with a single call to biquad_create_setup. coefficients must
    contain 5 coefficients for each section. The three feed-forward coefficients are
    specified first, followed by the two feedback coefficients. Returns a Biquad object.

    Returns: Biquad{T}
    """
    @eval begin
        function biquadcreate(coefficients::Vector{$T},  sections::Int)
            if length(coefficients) < 5*sections
                error("Incomplete biquad specification provided - coefficients must
                            contain 5 elements for each filter section")
            end
            setup = ccall(($(string("vDSP_biquad_CreateSetup", suff), libacc)),  Ptr{Cvoid},
                          (Ptr{$T}, UInt64),
                          coefficients, sections)
            return Biquad{$T}(setup, sections)
        end
    end


    """
    Filters an input array X with the specified Biquad filter and filter delay values provided
    in delays; only numelem elements are filtered. After execution, delays contains the final
    state data of the filter.

    Returns: Vector{T}
    """
    @eval begin
        function biquad(X::Vector{$T}, delays::Vector{$T}, numelem::Int, biquad::Biquad{$T})
            if length(delays) < (2*(biquad.sections)+2)
                error("Incomplete delay specification provided - delays must contain 2*M + 2
                                values where M is the number of sections in the biquad")
            end
            result::Vector{$T} = similar(X)
            ccall(($(string("vDSP_biquad", suff), libacc)),  Cvoid,
                  (Ptr{Cvoid},  Ptr{$T},  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64),
                  biquad.setup, delays, X,  1, result,  1,  numelem)
            return result
        end
    end


    """
    Frees all resources associated with a particular Biquad previously
    created through a call to biquad_create_setup. This is called automatically
    when the setup object is no longer visible to the garbage collector.

    Returns: Cvoid
    """
    @eval begin
        function biquaddestroy(biquad::Biquad{$T})
            ccall(($(string("vDSP_biquad_DestroySetup", suff), libacc)),  Cvoid,
                  (Ptr{Cvoid}, ),
                  biquad.setup)
        end
    end
end

## == WINDOW GENERATION == ##
for (T, suff) in ((Float32, ""), (Float64, "D"))

    """
    Generates a Blackman window of length 'length' and stores it in `result'. By
    default, a full window is created, but if flag=1, only the first (n+1)/2 points
    will be calculated.

    Returns: Vector{$T}
    """
    @eval begin
        function blackman!(result::Vector{$T},  length::Int, flag::Int=0)
            ccall(($(string("vDSP_blkman_window", suff), libacc)), Cvoid,
                  (Ptr{$T}, UInt64,  Int64),
                  result, length, flag)
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
            ccall(($(string("vDSP_hamm_window", suff), libacc)), Cvoid,
                  (Ptr{$T}, UInt64,  Int64),
                  result, length, flag)
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
            ccall(($(string("vDSP_hann_window", suff), libacc)), Cvoid,
                  (Ptr{$T}, UInt64,  Int64),
                  result, length, flag)
            return result
        end
    end

    """
    Alias function for hanning!

    Returns: Vector{$T}
    """
    @eval begin
        function hann!(result::Vector{$T}, length::Int, flag::Int=0)
            hanning!(result, length, flag)
        end
    end

end

"""
Generates a Blackman window of length 'length'. Default return type
is Vector{Float64}, but if rtype=Float32, Vector{Float32}
will be returned.

Returns: Vector{T}
"""
function blackman(length::Int, rtype::Type=Float64)
    result = Array{rtype}(undef, length)
    blackman!(result, length, 0)
end

"""
Generates a Hamming window of length 'length'. Default return type
is Vector{Float64}, but if rtype=Float32, Vector{Float32}
will be returned.

Returns: Vector{T}
"""
function hamming(length::Int, rtype::Type=Float64)
    result = Array{rtype}(undef, length)
    hamming!(result, length, 0)
end

"""
Generates a denormalized Hanning window of length 'length'. Default
return type is Vector{Float64}, but if rtype=Float32, Vector{Float32}
will be returned.

Returns: Vector{T}
"""
function hanning(length::Int, rtype::Type=Float64)
    result = Array{rtype}(undef, length)
    hanning!(result, length, 0)
end

"""
Alias function for hanning

Returns: Vector{T}
"""
function hann(length::Int, rtype::Type=Float64)
    hanning(length, rtype)
end


## == Discrete Cosine Transform (DCT) == ##
"""
Initializes a new DCT setup object. 'dct_type' must be 2, 3, 4 corresponding to type II, III and IV.
DCT 'length' must be equal to f*(2^n) where f = 1,3,5,15 and n >= 4. If you have a previous DCT setup
object, that can be passed in as 'previous'. The returned DCT setup will share the underlying data
storage of the previous setup object.

Returns: DFTSetup
"""
function plan_dct(length::Int,  dct_type::Int, previous=C_NULL)
    n = trailing_zeros(length)
    f = length >> n
    if dct_type < 2 &&  dct_type > 4
        error("DCT type ", dct_type, " is not supported. Only DCT types 2, 3 and 4 are supported")
    elseif !(n >= 4 && f in (1,3,5,15))
        error("Invalid DCT length. Length must be equal to f*(2^n) where f = 1,3,5,15 and n >= 4")
    end
    setup::Ptr{Cvoid} = ccall(("vDSP_DCT_CreateSetup", libacc), Ptr{Cvoid},
                             (Ptr{Cvoid}, UInt64, UInt64),
                             previous, length, dct_type)
    return DFTSetup{Float32}(setup, 0)
end


"""
Computes the DCT of a given input vector X using the parameters setup in
the DFTSetup object.

Returns: Vector{Float32}
"""
function dct(X::Vector{Float32}, setup::DFTSetup)
    result = similar(X)
    ccall(("vDSP_DCT_Execute", libacc),  Cvoid,
          (Ptr{Cvoid},  Ptr{Float32},  Ptr{Float32}),
          setup.setup,  X, result)
    return result
end


"""
Computes the DCT of a given input vector X using a DCT type 'dct_type' (defaults to type II).
This function does not require a separate call to dct_setup.

Returns: Vector{Float32}
"""
function dct(X::Vector{Float32}, dct_type::Int=2)
    setup = plan_dct(length(X), dct_type)
    return dct(X, setup)
end


"""
Deinitializes a DFTSetup object created by plan_dct
"""
function plan_destroy(setup::DFTSetup)
    ccall(("vDSP_DFT_DestroySetup", libacc), Cvoid,
          (Ptr{Cvoid},),
          setup.setup)
end


# const     FFT_FORWARD         = +1
# const     FFT_INVERSE         = -1
# const     SIGNAL_STRIDE       = 1
#
# typealias vDSP_Stride Clong
# typealias vDSP_Length Culong
# typealias FFTDirection Cint
#
# immutable DSPDoubleSplitComplex
#     realp::Ptr{Float64}
#     imagp::Ptr{Float64}
# end
#
# DSPDoubleSplitComplex(realp::Vector{Float64},imagp::Vector{Float64})=DSPDoubleSplitComplex(pointer(realp),pointer(imagp))
#
#
# function plan_fft(n,radix::Integer=2)
#     @assert isinteger(log2(n))
#     logn=round(Int,log2(n))
#     ccall(("vDSP_create_fftsetupD",libacc),Ptr{Cvoid},(Cuint,Cint),logn,radix)
# end
#
# function newfft(r::Vector{Complex128},plan)
#     n=length(r)
#     @assert isinteger(log2(n))
#     logn=round(Int,log2(n))
#
#     realp,imagp=real(r),imag(r)         # keep references to aCvoid garbage collection
#     vals=DSPDoubleSplitComplex(realp,imagp)
#     retr,reti=Array(Float64,n),Array(Float64,n)        # keep references to aCvoid garbage collection
#     ret=DSPDoubleSplitComplex(retr,reti)
#
#
#     ccall(("vDSP_fft_zopD",libacc),Cvoid,
#           (Ptr{Cvoid},DSPDoubleSplitComplex,vDSP_Stride,     DSPDoubleSplitComplex,vDSP_Stride,  vDSP_Length,FFTDirection),
#            plan,     vals,                 SIGNAL_STRIDE,   ret,                  SIGNAL_STRIDE,logn,       FFT_FORWARD)
#
#     retr+im*reti
# End
