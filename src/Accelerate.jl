module Accelerate

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

# TODO:
# div, fabs, fmod, remainder
# cosisin, sincos

for (T, suff) in ((Float64, ""), (Float32, "f"))

    # 1 arg functions
    for f in (:ceil,:floor,:sqrt,:rsqrt,:rec,
              :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
              :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
              :sinh,:cosh,:tanh,:asinh,:acosh,:atanh)
        @eval begin
            function ($f)(X::Array{$T})
                out = Array($T,size(X))
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
            function ($(symbol("$(f)!")))(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
        end
    end

    # renamed functions
    for (f,fa) in ((:trunc,:int),(:round,:nint),(:exponent,:logb))
        @eval begin
            function ($f)(X::Array{$T})
                out = Array($T,size(X))
                ccall(($(string("vv",fa,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
            function ($(symbol("$(f)!")))(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
        end
    end

    # 2 arg functions
    for f in (:copysign,:pow,:atan2)
        @eval begin
            function ($f)(X::Array{$T},Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array($T,size(X))
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Y,&length(X))
                out
            end
            function ($(symbol("$(f)!")))(X::Array{$T},Y::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Void,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Y,&length(X))
                out
            end

        end
    end
end


end # module
