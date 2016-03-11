## DSP.jl ##

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
