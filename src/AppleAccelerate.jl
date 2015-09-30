module AppleAccelerate

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

# TODO:
# remainder

for (T, suff) in ((Float64, ""), (Float32, "f"))

    # 1 arg functions
    for f in (:ceil,:floor,:sqrt,:rsqrt,:rec,
              :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
              :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
              :sinh,:cosh,:tanh,:asinh,:acosh,:atanh)
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = Array($T,size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
        end
    end

    # renamed 1 arg functions
    for (f,fa) in ((:trunc,:int),(:round,:nint),(:exponent,:logb),
                   (:abs,:fabs))
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = Array($T,size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
        end
    end

    # 2 arg functions
    for f in (:copysign,:atan2)
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array($T,size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Y,&length(X))
                out
            end
        end
    end

    # for some bizarre reason, vvpow/vvpowf reverse the order of arguments.
    for f in (:pow,)
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array($T,size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,Y,X,&length(X))
                out
            end
        end
    end


    # renamed 2 arg functions
    for (f,fa) in ((:rem,:fmod),(:fdiv,:div))
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array($T,size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Y,&length(X))
                out
            end
        end
    end

    # two-arg return
    for f in (:sincos,)
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out1 = Array($T,size(X))
                out2 = Array($T,size(X))
                ($f!)(out1, out2, X)
            end
            function ($f!)(out1::Array{$T}, out2::Array{$T}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out1,out2,X,&length(X))
                out1, out2
            end
        end
    end

    # complex return
    for f in (:cis,)
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = Array(Complex{$T},size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::Array{Complex{$T}}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{Complex{$T}},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
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



if VERSION>=v"0.4.0-rc1"
    tupletypelength(a)=length(a.parameters)
else
    tupletypelength(a)=length(a)
end


macro replaceBase(fs...)
    b = Expr(:block)
    for f in fs
        if f == :./
            fa = :fdiv
        elseif f == :.^
            fa = :pow
        else
            fa = f
        end
        e = quote
            if tupletypelength(methods($f).defs.sig) == 1
                (Base.$f)(X::Array{Float64}) = ($fa)(X)
                (Base.$f)(X::Array{Float32}) = ($fa)(X)
            else
                (Base.$f)(X::Array{Float64},Y::Array{Float64}) = ($fa)(X,Y)
                (Base.$f)(X::Array{Float32},Y::Array{Float32}) = ($fa)(X,Y)
            end
        end
        push!(b.args,e)
    end
    b
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



end # module
