## DSP.jl ##

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
            ccall(($(string("vDSP_conv", suff), libacc)),  Void,
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
            result = Array($T, length(X) + length(K) - 1)
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
            ccall(($(string("vDSP_conv", suff), libacc)),  Void,
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
            result = Array($T, length(X) + length(Y) - 1)
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


## == WINDOW GENERATION == ##
for (T, suff) in ((Float32, ""), (Float64, "D"))

    """
    Generates a Blackman window of length 'length' and stores it in `result'. By
    default, a full window is created, but if flag=1, only the first (n+1)/2 points
    will be calculated.

    Returns: Vector{$T}
    """
    @eval begin
        function blackman!(result::Vector{$T},  length::Integer, flag::Integer=0)
            ccall(($(string("vDSP_blkman_window", suff), libacc)), Void,
                  (Ptr{$T}, UInt64,  Int64),
                  result, length, flag)
            return result
        end
    end

    """
    Generates a Blackman window of length 'length'. Default return type
    is Vector{Float64}, but if rtype=Float32, Vector{Float32}
    will be returned.

    Returns: Vector{$T}
    """
    @eval begin
        function blackman(length::Integer, rtype::DataType=Float64)
            result::Vector{rtype} = Array(rtype, length)
            blackman!(result, length, 0)
        end
    end

    """
    Generates a Hamming window of length 'length' and stores it in `result'. By
    default, a full window is created, but if flag=1, only the first (n+1)/2 points
    will be calculated.

    Returns: Vector{$T}
    """
    @eval begin
        function hamming!(result::Vector{$T},  length::Integer, flag::Integer=0)
            ccall(($(string("vDSP_hamm_window", suff), libacc)), Void,
                  (Ptr{$T}, UInt64,  Int64),
                  result, length, flag)
            return result
        end
    end

    """
    Generates a Hamming window of length 'length'. Default return type
    is Vector{Float64}, but if rtype=Float32, Vector{Float32}
    will be returned.

    Returns: Vector{$T}
    """
    @eval begin
        function hamming(length::Integer, rtype::DataType=Float64)
            result::Vector{rtype} = Array(rtype, length)
            hamming!(result, length, 0)
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
        function hanning!(result::Vector{$T},  length::Integer, flag::Integer=0)
            ccall(($(string("vDSP_hann_window", suff), libacc)), Void,
                  (Ptr{$T}, UInt64,  Int64),
                  result, length, flag)
            return result
        end
    end

    """
    Generates a denormalized Hanning window of length 'length'. Default
    return type is Vector{Float64}, but if rtype=Float32, Vector{Float32}
    will be returned.

    Returns: Vector{$T}
    """
    @eval begin
        function hanning(length::Integer, rtype::DataType=Float64)
            result::Vector{rtype} = Array(rtype, length)
            hanning!(result, length, 0)
        end
    end

    """
    Alias function for hanning!

    Returns: Vector{$T}
    """
    @eval begin
        function hann!(result::Vector{$T}, length::Integer, flag::Integer=0)
            hanning!(result, length, flag)
        end
    end

    """
    Alias function for hanning

    Returns: Vector{$T}
    """
    @eval begin
        function hann(length::Integer, rtype::DataType=Float64)
            hanning(length, rtype)
        end
    end

end


function plan_dct(n,k::Integer)
    @assert isinteger(Base.log2(n))
    @assert 2≤k≤4
    ccall(("vDSP_DCT_CreateSetup",libacc),Ptr{Void},(Ptr{Void},Cint,Cint),C_NULL,n,k)
end


function dct(r::Vector{Float32},plan)
    n=length(r)
    @assert isinteger(Base.log2(n))
    out=Array(Float32,n)
    ccall(("vDSP_DCT_Execute",libacc),Void,(Ptr{Void},Ptr{Float32},Ptr{Float32}),plan,r,out)
    out
end

dct(r::Vector{Float32},k::Integer=2)=dct(r,plan_dct(length(r),k))




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
#     ccall(("vDSP_create_fftsetupD",libacc),Ptr{Void},(Cuint,Cint),logn,radix)
# end
#
# function newfft(r::Vector{Complex128},plan)
#     n=length(r)
#     @assert isinteger(log2(n))
#     logn=round(Int,log2(n))
#
#     realp,imagp=real(r),imag(r)         # keep references to avoid garbage collection
#     vals=DSPDoubleSplitComplex(realp,imagp)
#     retr,reti=Array(Float64,n),Array(Float64,n)        # keep references to avoid garbage collection
#     ret=DSPDoubleSplitComplex(retr,reti)
#
#
#     ccall(("vDSP_fft_zopD",libacc),Void,
#           (Ptr{Void},DSPDoubleSplitComplex,vDSP_Stride,     DSPDoubleSplitComplex,vDSP_Stride,  vDSP_Length,FFTDirection),
#            plan,     vals,                 SIGNAL_STRIDE,   ret,                  SIGNAL_STRIDE,logn,       FFT_FORWARD)
#
#     retr+im*reti
# End
