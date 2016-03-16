## DSP.jl ##

immutable DFTSetup{T}
    setup::Ptr{Void}
    direction::Integer
end

for (T, suff) in ((Float64, "D"), (Float32, ""))

    @eval begin
        function conv(X::Array{$T}, K::Array{$T})
            ksize = length(K)
            xsize = length(X)
            rsize = xsize + ksize - 1
            x_padded::Array{$T} = [zeros($T, ksize-1); X; zeros($T, ksize)]
            result = Array($T, rsize)
            ccall(($(string("vDSP_conv", suff), libacc)),  Void,
                  (Ptr{$T}, Int64,  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64, UInt64),
                  x_padded, 1, pointer(K, ksize), -1, result, 1,  rsize, ksize)
            return result
        end
    end


    @eval begin
        function xcorr(X::Array{$T}, K::Array{$T})
            ksize = length(K)
            xsize = length(X)
            rsize = xsize + ksize - 1
            x_padded::Array{$T} = [zeros($T, ksize-1); X; zeros($T, ksize)]
            result = Array($T, rsize)
            ccall(($(string("vDSP_conv", suff), libacc)),  Void,
                  (Ptr{$T}, Int64,  Ptr{$T},  Int64,  Ptr{$T},  Int64, UInt64, UInt64),
                  x_padded, 1, K, 1, result, 1,  rsize, ksize)
            return result
        end
    end

    @eval begin
        function xcorr(X::Array{$T})
            return xcorr(X, X)
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
function plan_dct(length::Integer,  dct_type::Integer, previous=C_NULL)
    if dct_type < 2 &&  dct_type > 4
        error("DCT type ", dct_type, " is not supported. Only DCT types 2, 3 and 4 are supported")
    elseif !(isinteger(Base.log2(length)))
        error("Invalid DCT length. Length must be equal to f*(2^n) where f = 1,3,5,15 and n >= 4")
    end
    setup::Ptr{Void} = ccall(("vDSP_DCT_CreateSetup", libacc), Ptr{Void},
                             (Ptr{Void}, UInt64, UInt64),
                             previous, length, dct_type)
    return DFTSetup{Float32}(setup, 0)
end


"""
Computes the DCT of a given input vector X using the parameters setup in
the DFTSetup object.

Returns: Vector{Float32}
"""
function dct(X::Vector{Float32}, setup::DFTSetup)
    Y = similar(X)
    ccall(("vDSP_DCT_Execute", libacc),  Void,
          (Ptr{Void},  Ptr{Float32},  Ptr{Float32}),
          setup.setup,  X, Y)
    return Y
end


"""
Computes the DCT of a given input vector X using a DCT Type 'dct_type' (defaults to type II).
This function does not require a separate call to dct_setup.

Returns: Vector{Float32}
"""
function dct(X::Vector{Float32}, dct_type::Integer=2)
    setup = plan_dct(length(X), dct_type)
    return dct(X, setup)
end


"""
Deinitializes a DFTSetup object created by plan_dct
"""
function plan_destroy(setup::DFTSetup)
    ccall(("vDSP_DFT_DestroySetup", libacc), Void,
          (Ptr{Void},),
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
