module Accelerate

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

# div, fabs, fmod, remainder
# cosisin, sincos

for f in (:ceil,:floor,:sqrt,:rsqrt,:rec,
          :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,:logb,
          :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
          :sinh,:cosh,:tanh,:asinh,:acosh,:atanh),
    (t, s) in ((Float64, ""), (Float32, "f"))
    @eval begin
        function ($f)(X::Array{$t})
            n = length(X)
            Y = Array($t,n)
            ccall(($(string("vv",f,s)),libacc),Void,(Ptr{$t},Ptr{$t},Ptr{Cint}),Y,X,&n)
            Y
        end
        function ($(symbol("$(f)!")))(X::Array{$t})
            n = length(X)
            ccall(($(string("vv",f,s)),libacc),Void,(Ptr{$t},Ptr{$t},Ptr{Cint}),X,X,&n)
        end
    end
end

for (f,fa) in ((:trunc,:int),(:round,:nint))
    @eval begin
        function ($f)(X::Array{Float64})
            n = length(X)
            Y = Array(Float64,n)
            ccall(($(string("vv",fa)),libacc),Void,(Ptr{Float64},Ptr{Float64},Ptr{Cint}),Y,X,&n)
            Y
        end
        function ($f)(X::Array{Float32})
            n = length(X)
            Y = Array(Float32,n)
            ccall(($(string("vv",fa,"f")),libacc),Void,(Ptr{Float32},Ptr{Float32},Ptr{Cint}),Y,X,&n)
            Y
        end
    end
end



for f in (:copysign,:pow,:atan2)
    @eval begin
        function ($f)(X::Array{Float64},Y::Array{Float64})
            size(X) == size(Y) || error("dimensions must match")
            Z = Array(Float64,size(X))
            ccall(($(string("vv",f)),libacc),Void,
                  (Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Cint}),Z,X,Y,&length(X))
            Z
        end
        function ($f)(X::Array{Float32},Y::Array{Float32})
            size(X) == size(Y) || error("dimensions must match")
            Y = Array(Float64,size(X))
            ccall(($(string("vv",f,"f")),libacc),Void,
                  (Ptr{Float32},Ptr{Float32},Ptr{Float32},Ptr{Cint}),Z,X,Y,&length(X))
            Y
        end
    end
end


end # module
