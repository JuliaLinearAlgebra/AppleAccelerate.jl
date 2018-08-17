## Array.jl ##

for (T, suff) in ((Float64, ""), (Float32, "f"))

    # 1 arg functions
    for f in (:ceil,:floor,:sqrt,:rsqrt,:rec,
              :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
              :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
              :sinh,:cosh,:tanh,:asinh,:acosh,:atanh)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = similar(X)
                ($f!)(out, X)
            end
            function ($f!)(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Ref{Cint}(length(X)))
                out
            end
        end
    end

    # renamed 1 arg functions
    for (f,fa) in ((:trunc,:int),(:round,:nint),(:exponent,:logb),
                   (:abs,:fabs))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = similar(X)
                ($f!)(out, X)
            end
            function ($f!)(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Ref{Cint}(length(X)))
                out
            end
        end
    end

    # 2 arg functions
    for f in (:copysign,)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = similar(X)
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Y,Ref{Cint}(length(X)))
                out
            end
        end
    end

    # for some bizarre reason, vvpow/vvpowf reverse the order of arguments.
    for f in (:pow,)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = similar(X)
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,Y,X,Ref{Cint}(length(X)))
                out
            end
        end
    end


    # renamed 2 arg functions
    for (f,fa) in ((:rem,:fmod),(:fdiv,:div),(:atan,:atan2))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = similar(X)
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out,X,Y,Ref{Cint}(length(X)))
                out
            end
        end
    end

    # two-arg return
    for f in (:sincos,)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out1 = similar(X)
                out2 = similar(X)
                ($f!)(out1, out2, X)
            end
            function ($f!)(out1::Array{$T}, out2::Array{$T}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ptr{Cint}),out1,out2,X,Ref{Cint}(length(X)))
                out1, out2
            end
        end
    end

    # complex return
    for (f,fa) in ((:cis,:cosisin),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = Array{Complex{$T}}(undef, size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::Array{Complex{$T}}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{Complex{$T}},Ptr{$T},Ptr{Cint}),out,X,Ref{Cint}(length(X)))
                out
            end
        end
    end
end

# Functions over single vectors that return scalars/tuples
for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, fa) in ((:maximum, :maxv), (:minimum, :minv), (:mean, :meanv),
                    (:meansqr, :measqv), (:meanmag,  :meamgv), (:sum, :sve))
        @eval begin
            function ($f)(X::Vector{$T})
                val = Ref{$T}(0.0)
                ccall(($(string("vDSP_", fa, suff), libacc)),  Cvoid,
                      (Ptr{$T}, Int64,  Ref{$T}, UInt64),
                      X, 1, val, length(X))
                return val[]
            end
        end
    end

    for (f, fa) in ((:findmax, :maxvi), (:findmin, :minvi))
        @eval begin
            function ($f)(X::Vector{$T})
                index = Ref{Int}(0)
                val = Ref{$T}(0.0)
                ccall(($(string("vDSP_", fa, suff), libacc)),  Cvoid,
                      (Ptr{$T}, Int64,  Ref{$T}, Ref{Int}, UInt64),
                      X, 1, val, index, length(X))
                return (val[], index[]+1)
            end
        end
    end
end

# Element-wise operations over two vectors
for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, name) in ((:vadd, "addition"),  (:vsub, "subtraction"),
                      (:vdiv, "division"), (:vmul, "multiplication"))
        f! = Symbol("$(f)!")

        @eval begin
            @doc """
            `$($f!)(result::Vector{$($T)}, X::Vector{$($T)}, Y::Vector{$($T)})`

            Implements element-wise **$($name)** over two **Vector{$($T)}** and overwrites
            the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
            """
            function ($f!)(result::Vector{$T}, X::Vector{$T}, Y::Vector{$T})
                ccall(($(string("vDSP_", f, suff), libacc)),  Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T},  Int64, Ptr{$T}, Int64,  UInt64),
                      Y, 1, X, 1, result, 1, length(result))
                return result
            end
        end

        @eval begin
            @doc """
            `$($f)(X::Vector{$($T)}, Y::Vector{$($T)})`

            Implements element-wise **$($name)** over two **Vector{$($T)}**. Allocates
            memory to store result. *Returns:* **Vector{$($T)}**
            """
            function ($f)(X::Vector{$T}, Y::Vector{$T})
                result = similar(X)
                ($f!)(result, X, Y)
                return result
            end
        end
    end
end
