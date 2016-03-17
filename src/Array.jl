## Array.jl ##

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
    for (f,fa) in ((:cis,:cosisin),)
        f! = symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = Array(Complex{$T},size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::Array{Complex{$T}}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Void,
                      (Ptr{Complex{$T}},Ptr{$T},Ptr{Cint}),out,X,&length(X))
                out
            end
        end
    end
end

for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, fa) in ((:maximum, :maxv), (:minimum, :minv), (:mean, :meanv),
                    (:meansqr, :measqv), (:meanmag,  :meamgv), (:sum, :sve))
        @eval begin
            function ($f)(X::Array{$T})
                val::$T = 0.0
                ccall(($(string("vDSP_", fa, suff), libacc)),  Void,
                      (Ptr{$T}, Int64,  Ref{$T}, UInt64),
                      X, 1, val, length(X))
                return val
            end
        end
    end

    for (f, fa) in ((:findmax, :maxvi), (:findmin, :minvi))
        @eval begin
            function ($f)(X::Array{$T})
                index::Int = 0
                val::$T = 0.0
                ccall(($(string("vDSP_", fa, suff), libacc)),  Void,
                      (Ptr{$T}, Int64,  Ref{$T}, Ref{$T}, UInt64),
                      X, 1, val, index, length(X))
                return (val, index)
            end
        end
    end
end
