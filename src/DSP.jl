## DSP.jl ##

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
            return Biquad($T, setup, sections)
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

"""
    blackman(length, [rtype=Float64])

Generates a Blackman window of length 'length'. Default return type
is Vector{Float64}, but if rtype=Float32, Vector{Float32}
will be returned.
"""
function blackman(length::Int, rtype::DataType=Float64)
    result::Vector{rtype} = Array{rtype}(undef, length)
    blackman!(result, length, 0)
end

"""
    hamming(length, [rtype=Float64])

Generates a Hamming window of length 'length'. Default return type
is Vector{Float64}, but if rtype=Float32, Vector{Float32}
will be returned.
"""
function hamming(length::Int, rtype::DataType=Float64)
    result::Vector{rtype} = Array{rtype}(undef, length)
    hamming!(result, length, 0)
end

"""
    hanning(length, [rtype=Float64])

Generates a denormalized Hanning window of length 'length'. Default
return type is Vector{Float64}, but if rtype=Float32, Vector{Float32}
will be returned.
"""
function hanning(length::Int, rtype::DataType=Float64)
    result::Vector{rtype} = Array{rtype}(undef, length)
    hanning!(result, length, 0)
end

"""
Alias function for `hanning`
"""
function hann(length::Int, rtype::DataType=Float64)
    hanning(length, rtype)
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


## == Discrete Cosine Transform (DCT) == ##
"""
Initializes a new DCT setup object. 'dct_type' must be 2, 3, 4 corresponding to Type II, III and IV.
DCT 'length' must be equal to f*(2^n) where f = 1,3,5,15 and n >= 4. If you have a previous DCT setup
object, that can be passed in as 'previous'. The returned DCT setup will share the underlying data
storage of the previous setup object.

Returns: DFTSetup
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
    return DFTSetup(Float32, setup, 0)
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
Computes the DCT of a given input vector X using a DCT Type 'dct_type' (defaults to type II).
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


const     FFT_FORWARD         = 1
const     FFT_INVERSE         = -1
const     SIGNAL_STRIDE       = 1

mutable struct FFTSetup{T}
    plan::Ptr{Cvoid}

    function FFTSetup{Float64}(n::Integer, radix::Integer = 2)
        @assert ispow2(n) "n must be a power of 2"
        logn = trailing_zeros(n)
        plan = ccall(("vDSP_create_fftsetupD", libacc), Ptr{Cvoid}, (Culong, Cint), logn, radix)
        setup = new{Float64}(plan)
        finalizer(destroy_fftsetup, setup)
        setup
    end

    function FFTSetup{Float32}(n::Integer, radix::Integer = 2)
        @assert ispow2(n) "n must be a power of 2"
        logn = trailing_zeros(n)
        plan = ccall(("vDSP_create_fftsetup", libacc), Ptr{Cvoid}, (Culong, Cint), logn, radix)
        setup = new{Float32}(plan)
        finalizer(destroy_fftsetup, setup)
        setup
    end
end

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
    ccall(("vDSP_destroy_fftsetupD", libacc), Cvoid, (Ptr{Cvoid},), setup.plan)
end

function destroy_fftsetup(setup::FFTSetup{Float32})
    ccall(("vDSP_destroy_fftsetup", libacc), Cvoid, (Ptr{Cvoid},), setup.plan)
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
        ccall(("vDSP_fft_zopD", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPDoubleSplitComplex}, Clong,
               Ref{DSPDoubleSplitComplex}, Clong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, output, SIGNAL_STRIDE, logn, direction)
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
        ccall(("vDSP_fft_zop", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPSplitComplex}, Clong,
               Ref{DSPSplitComplex}, Clong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, output, SIGNAL_STRIDE, logn, direction)
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
        ccall(("vDSP_fft2d_zopD", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPDoubleSplitComplex}, Clong, Clong,
               Ref{DSPDoubleSplitComplex}, Clong, Clong, Culong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, 0, output, SIGNAL_STRIDE, 0,
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
        ccall(("vDSP_fft2d_zop", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPSplitComplex}, Clong, Clong,
               Ref{DSPSplitComplex}, Clong, Clong, Culong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, 0, output, SIGNAL_STRIDE, 0,
              log2nr, log2nc, direction)
    end

    return complex.(retr, reti)
end

# --- Public API: fft (forward FFT) ---

fft(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d(x, setup, FFT_FORWARD)
fft(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d(x, setup, FFT_FORWARD)
fft(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = fft(x, plan_fft(x))
fft(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = fft(x, plan_fft(x))

# --- Public API: bfft (backward/unnormalized inverse FFT) ---

bfft(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d(x, setup, FFT_INVERSE)
bfft(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d(x, setup, FFT_INVERSE)
bfft(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft(x, plan_fft(x))
bfft(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft(x, plan_fft(x))

# --- Public API: ifft (normalized inverse FFT) ---

ifft(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = bfft(x, setup) ./ length(x)
ifft(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = bfft(x, setup) ./ length(x)
ifft(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft(x, plan_fft(x))
ifft(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft(x, plan_fft(x))

# --- Internal in-place 1D complex FFT ---

function _fft1d!(x::Vector{ComplexF64}, setup::FFTSetup{Float64}, direction::Int)
    n = length(x)
    @assert ispow2(n) "length of input must be a power of 2"
    logn = trailing_zeros(n)

    realp = real.(x)
    imagp = imag.(x)

    GC.@preserve realp imagp begin
        c = DSPDoubleSplitComplex(realp, imagp)
        ccall(("vDSP_fft_zipD", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPDoubleSplitComplex}, Clong, Culong, Cint),
              setup.plan, c, SIGNAL_STRIDE, logn, direction)
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
        ccall(("vDSP_fft_zip", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPSplitComplex}, Clong, Culong, Cint),
              setup.plan, c, SIGNAL_STRIDE, logn, direction)
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
        ccall(("vDSP_fft2d_zipD", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPDoubleSplitComplex}, Clong, Clong, Culong, Culong, Cint),
              setup.plan, c, SIGNAL_STRIDE, 0, log2nr, log2nc, direction)
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
        ccall(("vDSP_fft2d_zip", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPSplitComplex}, Clong, Clong, Culong, Culong, Cint),
              setup.plan, c, SIGNAL_STRIDE, 0, log2nr, log2nc, direction)
    end

    @inbounds for i in eachindex(x)
        x[i] = complex(realp[i], imagp[i])
    end
    return x
end

# --- Public API: fft! (in-place forward FFT) ---

fft!(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d!(x, setup, FFT_FORWARD)
fft!(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d!(x, setup, FFT_FORWARD)
fft!(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = fft!(x, plan_fft(x))
fft!(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = fft!(x, plan_fft(x))

# --- Public API: bfft! (in-place backward/unnormalized inverse FFT) ---

bfft!(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft1d!(x, setup, FFT_INVERSE)
bfft!(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _fft2d!(x, setup, FFT_INVERSE)
bfft!(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft!(x, plan_fft(x))
bfft!(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = bfft!(x, plan_fft(x))

# --- Public API: ifft! (in-place normalized inverse FFT) ---

function ifft!(x::Vector{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}}
    bfft!(x, setup)
    x ./= length(x)
end
function ifft!(x::Matrix{Complex{T}}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}}
    bfft!(x, setup)
    x ./= length(x)
end
ifft!(x::Vector{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft!(x, plan_fft(x))
ifft!(x::Matrix{Complex{T}}) where {T<:Union{Float32,Float64}} = ifft!(x, plan_fft(x))

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
        ccall(("vDSP_fft_zropD", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPDoubleSplitComplex}, Clong,
               Ref{DSPDoubleSplitComplex}, Clong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, output, SIGNAL_STRIDE, logn, FFT_FORWARD)
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
        ccall(("vDSP_fft_zrop", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPSplitComplex}, Clong,
               Ref{DSPSplitComplex}, Clong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, output, SIGNAL_STRIDE, logn, FFT_FORWARD)
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
        ccall(("vDSP_fft_zropD", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPDoubleSplitComplex}, Clong,
               Ref{DSPDoubleSplitComplex}, Clong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, output, SIGNAL_STRIDE, logn, FFT_INVERSE)
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
        ccall(("vDSP_fft_zrop", libacc), Cvoid,
              (Ptr{Cvoid}, Ref{DSPSplitComplex}, Clong,
               Ref{DSPSplitComplex}, Clong, Culong, Cint),
              setup.plan, input, SIGNAL_STRIDE, output, SIGNAL_STRIDE, logn, FFT_INVERSE)
    end

    result = Vector{Float32}(undef, n)
    @inbounds for k in 1:half
        result[2k - 1] = out_realp[k]
        result[2k] = out_imagp[k]
    end
    return result
end

# --- Public API: rfft (forward real FFT) ---

plan_rfft(x::Vector{T}) where {T<:Union{Float32,Float64}} = FFTSetup{T}(length(x))

rfft(x::Vector{T}, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _rfft1d(x, setup)
rfft(x::Vector{T}) where {T<:Union{Float32,Float64}} = rfft(x, plan_rfft(x))

# --- Public API: brfft (backward/unnormalized inverse real FFT) ---

brfft(X::Vector{Complex{T}}, n::Int, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = _brfft1d(X, n, setup)
brfft(X::Vector{Complex{T}}, n::Int) where {T<:Union{Float32,Float64}} = brfft(X, n, FFTSetup{T}(n))

# --- Public API: irfft (normalized inverse real FFT) ---

irfft(X::Vector{Complex{T}}, n::Int, setup::FFTSetup{T}) where {T<:Union{Float32,Float64}} = brfft(X, n, setup) ./ n
irfft(X::Vector{Complex{T}}, n::Int) where {T<:Union{Float32,Float64}} = irfft(X, n, FFTSetup{T}(n))
