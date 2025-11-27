## Array.jl ##

for (T, suff) in ((Float64, ""), (Float32, "f"))

    # 1-arg functions
    onearg_funcs = (:ceil,:floor,:sqrt,:rsqrt,:rec,
              :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
              :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
              :sinh,:cosh,:tanh,:asinh,:acosh,:atanh)
    # Combine this with the list of renamed 1-arg funcs
    onearg_funcs = (
        ((x,x) for x in onearg_funcs)...,
        (:trunc,:int),
        (:round,:nint),
        (:exponent,:logb),
        (:abs,:fabs)
    )

    for (f,fa) in onearg_funcs
        f! = Symbol("$(f)!")
        @eval begin
            # Allocating variant
            function ($f)(X::Array{$T})
                out = Array{$T}(undef, size(X))
                ($f!)(out, X)
            end

            # In-place mutating variant
            function ($f!)(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ref{Cint}),out,X,length(X))
                out
            end
        end
    end


    # 2 arg functions
    twoarg_funcs = (
        (:copysign,:copysign),
        (:rem,:fmod),
        (:div_float,:div),
        (:atan,:atan2)
    )

    for (f, fa) in twoarg_funcs
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ref{Cint}),out,X,Y,length(X))
                out
            end
        end
    end

    # for some bizarre reason, vvpow/vvpowf reverse the order of arguments.
    for (f, fa) in ((:pow,:pow),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, Y::Array{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, Y::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ref{Cint}),out,Y,X,length(X))
                out
            end
        end
    end

    # two-arg return
    for (f, fa) in ((:sincos,:sincos),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out1 = Array{$T}(undef, size(X))
                out2 = Array{$T}(undef, size(X))
                ($f!)(out1, out2, X)
            end
            function ($f!)(out1::Array{$T}, out2::Array{$T}, X::Array{$T})
                ccall(($(string("vv",f,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ptr{$T},Ref{Cint}),out1,out2,X,length(X))
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
                      (Ptr{Complex{$T}},Ptr{$T},Ref{Cint}),out,X,length(X))
                out
            end
        end
    end

    for (f, fa) in onearg_funcs
        f! = Symbol("$(f)!")
        @eval begin
            # Broadcasting override such that f.(X) turns into f(X)
            Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N}}}) where {Style, Axes, N} = ($f)(bc.args...)
            Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N}}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
        end
        if T == Float32
            @eval Base.broadcasted(::typeof($f), arg::Union{Array{F,N},Base.Broadcast.Broadcasted}) where {N,F<:Union{Float32,Float64}} = ($f)(maybecopy(arg))
        end
    end
    for (f, fa) in (twoarg_funcs...,(:pow,:pow))
        f! = Symbol("$(f)!")
        @eval begin
            # Broadcasting override such that f.(X) turns into f(X)
            Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N},Array{$T,N}}}) where {Style, Axes, N} = ($f)(bc.args...)
            Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T,N},Array{$T,N}}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
        end
        if T == Float32
            @eval Base.broadcasted(::typeof($f), arg1::Union{Array{F, N},Base.Broadcast.Broadcasted}, arg2::Union{Array{F, N},Base.Broadcast.Broadcasted}) where {N,F<:Union{Float32,Float64}} = ($f)(maybecopy(arg1), maybecopy(arg2))
        end
    end
end

# Functions over single vectors that return scalars/tuples
for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, fa) in ((:maximum, :maxv), (:minimum, :minv), (:mean, :meanv),
                    (:meanmag,  :meamgv), (:meansqr, :measqv), (:meanssqr, :mvessq),
                    (:sum, :sve), (:summag, :svemg), (:sumsqr, :svesq),
                    (:sumssqr, :svs))
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
            """ ->
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
            """ ->
            function ($f)(X::Vector{$T}, Y::Vector{$T})
                result = similar(X)
                ($f!)(result, X, Y)
                return result
            end
        end

        @eval begin
            # Broadcasting override such that f.(X) turns into f(X)
            Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N},Array{$T,N}}}) where {Style, Axes, N} = ($f)(bc.args...)
            Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T,N},Array{$T,N}}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
            Base.broadcasted(::typeof($f), arg1::Union{Array{$T, N},Base.Broadcast.Broadcasted}, arg2::Union{Array{$T, N},Base.Broadcast.Broadcasted}) where {N} = ($f)(maybecopy(arg1), maybecopy(arg2))
        end
    end
end

# Element-wise operations over a vector and a scalar
for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, name) in ((:vsadd, "addition"), (:vsdiv, "division"), (:vsmul, "multiplication"))
        f! = Symbol("$(f)!")

        @eval begin
            @doc """
            `$($f!)(result::Vector{$($T)}, X::Vector{$($T)}, c::$($T))`

            Implements vector-scalar **$($name)** over **Vector{$($T)}** and $($T) and overwrites
            the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
            """ ->
            function ($f!)(result::Vector{$T}, X::Vector{$T}, c::$T)
                ccall(($(string("vDSP_", f, suff), libacc)),  Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64,  UInt64),
                      X, 1, Ref(c), result, 1, length(result))
                return result
            end
        end

        @eval begin
            @doc """
            `$($f)(X::Vector{$($T)}, c::$($T))`

            Implements vector-scalar **$($name)** over **Vector{$($T)}** and $($T). Allocates
            memory to store result. *Returns:* **Vector{$($T)}**
            """ ->
            function ($f)(X::Vector{$T}, c::$T)
                result = similar(X)
                ($f!)(result, X, c)
                return result
            end
        end
    end
    f = :vssub
    f! = Symbol("$(f)!")

    @eval begin
        @doc """
        `$($f!)(result::Vector{$($T)}, X::Vector{$($T)}, c::$($T))`

        Implements vector-scalar **subtraction** over **Vector{$($T)}** and $($T) and overwrites
        the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
        """ ->
        function ($f!)(result::Vector{$T}, X::Vector{$T}, c::$T)
            ccall(($(string("vDSP_vsadd", suff), libacc)),  Cvoid,
                    (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64,  UInt64),
                    X, 1, Ref(-c), result, 1, length(result))
            return result
        end
    end

    @eval begin
        @doc """
        `$($f)(X::Vector{$($T)}, c::$($T))`

        Implements vector-scalar **subtraction** over **Vector{$($T)}** and $($T). Allocates
        memory to store result. *Returns:* **Vector{$($T)}**
        """ ->
        function ($f)(X::Vector{$T}, c::$T)
            result = similar(X)
            ($f!)(result, X, c)
            return result
        end
    end

    f = :svsub
    f! = Symbol("$(f)!")

    @eval begin
        @doc """
        `$($f!)(result::Vector{$($T)}, X::Vector{$($T)}, c::$($T))`

        Implements vector-scalar **subtraction** over $($T) and **Vector{$($T)}** and overwrites
        the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
        """ ->
        function ($f!)(result::Vector{$T}, X::Vector{$T}, c::$T)
            ccall(($(string("vDSP_vsadd", suff), libacc)),  Cvoid,
                    (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64,  UInt64),
                    -X, 1, Ref(c), result, 1, length(result))
            return result
        end
    end

    @eval begin
        @doc """
        `$($f)(X::Vector{$($T), c::$($T)})`

        Implements vector-scalar **subtraction** over $($T) and  **Vector{$($T)}**. Allocates
        memory to store result. *Returns:* **Vector{$($T)}**
        """ ->
        function ($f)(X::Vector{$T}, c::$T)
            result = similar(X)
            ($f!)(result, X, c)
            return result
        end
    end

    for f in (:vsadd, :vssub, :vsdiv, :vsmul)
        f! = Symbol("$(f)!")
        @eval begin
            # Broadcasting override such that f.(X) turns into f(X)
            Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N},$T}}) where {Style, Axes, N} = ($f)(bc.args...)
            Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T,N},$T}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
            Base.broadcasted(::typeof($f), arg1::Union{Array{$T, N},Base.Broadcast.Broadcasted}, arg2::$T) where {N} = ($f)(maybecopy(arg1), arg2)
        end
    end

    f = :svsub
    f! = Symbol("$(f)!")

    @eval begin
        # Broadcasting override such that f.(X) turns into f(X)
        Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N}, $T}}) where {Style, Axes, N} = ($f)(bc.args...)
        Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T,N}, $T}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
        Base.broadcasted(::typeof($f), arg1::Union{Array{$T, N},Base.Broadcast.Broadcasted}, arg2::$T) where {N} = ($f)(maybecopy(arg1), arg2)
    end
end

