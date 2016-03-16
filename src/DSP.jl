## DSP.jl ##

## === TYPES === ##
immutable Biquad{T}
    setup::Ptr{Void}
    sections::Integer

end


## === FUNCTIONS == ##
"""
Convolution/Cross-Correlation
"""
for (T, suff) in ((Float64, "D"), (Float32, ""))

    """
    Returns the convolution of a 1-D vector X with a 1-D kernel vector K.
    Returns vector of size length(X) + length(K) -1.
    """
    @eval begin
        function conv(X::Vector{$T}, K::Vector{$T})
            ksize = length(K)
            xsize = length(X)
            rsize = xsize + ksize - 1
            x_padded::Vector{$T} = [zeros($T, ksize-1); X; zeros($T, ksize)]
            result = Vector($T, rsize)
            ccall(($(string("vDSP_conv", suff), libacc)),  Void,
                  (Ptr{$T}, Int64,  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64, UInt64),
                  x_padded, 1, pointer(K, ksize), -1, result, 1,  rsize, ksize)
            return result
        end
    end


    """
    Cross-correlation of two 1-D signals X & Y. Returned vector is of size
    length(X) + length(Y) - -1.
    """
    @eval begin
        function xcorr(X::Vector{$T}, Y::Vector{$T})
            ysize = length(Y)
            xsize = length(X)
            rsize = xsize + ysize - 1
            x_padded::Vector{$T} = [zeros($T, ysize-1); X; zeros($T, ysize)]
            result = Vector($T, rsize)
            ccall(($(string("vDSP_conv", suff), libacc)),  Void,
                  (Ptr{$T}, Int64,  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64, UInt64),
                  x_padded, 1, Y, 1, result, 1,  rsize, ysize)
            return result
        end
    end


    """
    Performs auto-correlation of a 1-D signal with itself. Returned vector
    is of size 2*length(X) - 1
    """
    @eval begin
        function xcorr(X::Vector{$T})
            return xcorr(X, X)
        end
    end
end


"""
Biquadratic/IIR filtering
"""
for (T, suff) in ((Float64, "D"), (Float32, ""))

    """
    Initializes a vDSP_biquad_setup for use with vDSP_biquad. A multi-section filter
    can be initialized with a single call to biquad_create_setup. coefficients must
    contain 5 coefficients for each section. The three feed-forward coefficients are
    specified first, followed by the two feedback coefficients. Returns a Biquad object.
    """
    @eval begin
        function biquad_create_setup(coefficients::Vector{$T},  sections::Integer)
            if length(coefficients) != 5*sections
                error("Incomplete biquad specification provided - coefficients must
                            contain 5 elements for each filter section")
                return
            end
            setup = ccall(($(string("vDSP_biquad_CreateSetup", suff), libacc)),  Ptr{Void},
                          (Ptr{$T}, UInt64),
                          coefficients, sections)
            return Biquad{$T}(setup, sections)
        end
    end


    """
    Filters an input array X with the specified Biquad filter and filter delay values provided
    in delays; only numelem elements are filtered. After execution, delays contains the final
    state data of the filter.
    """
    @eval begin
        function biquad(X::Vector{$T}, delays::Vector{$T}, numelem::Integer, biquad::Biquad)
            if length(delays) != (2*(biquad.sections)+2)
                error("Incomplete delay specification provided - delays must contain 2*M +2
                                values where M is the number of sections in the biquad")
            end
            result::Vector{$T} = similar(X, $T, length(X))
            ccall(($(string("vDSP_biquad", suff), libacc)),  Void,
                  (Ptr{Void},  Ptr{$T},  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64),
                  biquad.setup, delays, X,  1, result,  1,  numelem)
        end
    end


    """
    Frees all resources associated with a particular Biquad previously
    created through a call to biquad_create_setup
    """
    @eval begin
        function biquad_destroy(biquad::Biquad{$T})
            ccall(($(string("vDSP_biquad_DestroySetup", suff), libacc)),  Void,
                  (Ptr{Void}, ),
                  biquad.setup)
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
# end
