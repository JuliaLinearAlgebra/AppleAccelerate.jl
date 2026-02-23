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

# ============================================================
# vDSP Unary vector operations
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, fa) in ((:vneg, :vneg), (:vnabs, :vnabs), (:vsq, :vsq),
                    (:vssq, :vssq), (:vfrac, :vfrac), (:vabs, :vabs))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::Vector{$T}, X::Vector{$T})
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      X, 1, result, 1, length(X))
                return result
            end
            function ($f)(X::Vector{$T})
                result = similar(X)
                ($f!)(result, X)
            end
        end
    end

    # vreverse: in-place only
    @eval begin
        function vreverse!(X::Vector{$T})
            ccall(($(string("vDSP_vrvrs", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, UInt64),
                  X, 1, length(X))
            return X
        end
        function vreverse(X::Vector{$T})
            Y = copy(X)
            vreverse!(Y)
        end
    end
end

# ============================================================
# vDSP Two-vector element-wise operations
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))

    for (f, fa) in ((:vmax, :vmax), (:vmin, :vmin),
                    (:vmaxmg, :vmaxmg), (:vminmg, :vminmg),
                    (:vdist, :vdist))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::Vector{$T}, X::Vector{$T}, Y::Vector{$T})
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      X, 1, Y, 1, result, 1, length(X))
                return result
            end
            function ($f)(X::Vector{$T}, Y::Vector{$T})
                result = similar(X)
                ($f!)(result, X, Y)
            end
        end
    end

    # vtmerg: tapered merge — C[n] = A[n] + (B[n] - A[n]) * n/(N-1)
    # Signature same as two-vector ops
    @eval begin
        function vtmerg!(result::Vector{$T}, X::Vector{$T}, Y::Vector{$T})
            ccall(($(string("vDSP_vtmerg", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  X, 1, Y, 1, result, 1, length(X))
            return result
        end
        function vtmerg(X::Vector{$T}, Y::Vector{$T})
            result = similar(X)
            vtmerg!(result, X, Y)
        end
    end
end

# ============================================================
# vDSP Scalar-vector divide: C[n] = A / B[n]
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function svdiv!(result::Vector{$T}, X::Vector{$T}, c::$T)
            ccall(($(string("vDSP_svdiv", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  Ref(c), X, 1, result, 1, length(X))
            return result
        end
        function svdiv(X::Vector{$T}, c::$T)
            result = similar(X)
            svdiv!(result, X, c)
        end
    end
end

# ============================================================
# vDSP Compound arithmetic: 3-vector ops (A, B, C → D)
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))

    # vam: (A+B)*C, vsbm: (A-B)*C
    for (f, fa) in ((:vam, :vam), (:vsbm, :vsbm))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, C::Vector{$T})
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      A, 1, B, 1, C, 1, result, 1, length(A))
                return result
            end
            function ($f)(A::Vector{$T}, B::Vector{$T}, C::Vector{$T})
                result = similar(A)
                ($f!)(result, A, B, C)
            end
        end
    end

    # venvlp: signal envelope — same 3-vector pattern
    @eval begin
        function venvlp!(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, C::Vector{$T})
            ccall(($(string("vDSP_venvlp", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, 1, B, 1, C, 1, result, 1, length(A))
            return result
        end
        function venvlp(A::Vector{$T}, B::Vector{$T}, C::Vector{$T})
            result = similar(A)
            venvlp!(result, A, B, C)
        end
    end

    # 4-vector ops (A, B, C, D → E)
    for (f, fa) in ((:vaam, :vaam), (:vsbsbm, :vsbsbm), (:vasbm, :vasbm))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, C::Vector{$T}, D::Vector{$T})
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      A, 1, B, 1, C, 1, D, 1, result, 1, length(A))
                return result
            end
            function ($f)(A::Vector{$T}, B::Vector{$T}, C::Vector{$T}, D::Vector{$T})
                result = similar(A)
                ($f!)(result, A, B, C, D)
            end
        end
    end

    # vpythg: sqrt((A-C)^2 + (B-D)^2) — same 4-vector pattern
    @eval begin
        function vpythg!(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, C::Vector{$T}, D::Vector{$T})
            ccall(($(string("vDSP_vpythg", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, 1, B, 1, C, 1, D, 1, result, 1, length(A))
            return result
        end
        function vpythg(A::Vector{$T}, B::Vector{$T}, C::Vector{$T}, D::Vector{$T})
            result = similar(A)
            vpythg!(result, A, B, C, D)
        end
    end

    # 2 vectors + 1 scalar → D: vasm: (A+B)*C_scalar, vsbsm: (A-B)*C_scalar
    for (f, fa) in ((:vasm, :vasm), (:vsbsm, :vsbsm))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, c::$T)
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                      A, 1, B, 1, Ref(c), result, 1, length(A))
                return result
            end
            function ($f)(A::Vector{$T}, B::Vector{$T}, c::$T)
                result = similar(A)
                ($f!)(result, A, B, c)
            end
        end
    end

    # vsma: A*B_scalar + C
    @eval begin
        function vsma!(result::Vector{$T}, A::Vector{$T}, b::$T, C::Vector{$T})
            ccall(($(string("vDSP_vsma", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(b), C, 1, result, 1, length(A))
            return result
        end
        function vsma(A::Vector{$T}, b::$T, C::Vector{$T})
            result = similar(A)
            vsma!(result, A, b, C)
        end
    end

    # vsmsa: A*B_scalar + C_scalar
    @eval begin
        function vsmsa!(result::Vector{$T}, A::Vector{$T}, b::$T, c::$T)
            ccall(($(string("vDSP_vsmsa", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(b), Ref(c), result, 1, length(A))
            return result
        end
        function vsmsa(A::Vector{$T}, b::$T, c::$T)
            result = similar(A)
            vsmsa!(result, A, b, c)
        end
    end

    # vaddsub: returns (A+B, A-B) as two vectors
    @eval begin
        function vaddsub!(add_result::Vector{$T}, sub_result::Vector{$T}, A::Vector{$T}, B::Vector{$T})
            ccall(($(string("vDSP_vaddsub", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, 1, B, 1, add_result, 1, sub_result, 1, length(A))
            return (add_result, sub_result)
        end
        function vaddsub(A::Vector{$T}, B::Vector{$T})
            add_result = similar(A)
            sub_result = similar(A)
            vaddsub!(add_result, sub_result, A, B)
        end
    end
end

# ============================================================
# vDSP Scalar reductions: dot product, distance squared
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function dot(X::Vector{$T}, Y::Vector{$T})
            val = Ref{$T}(0.0)
            ccall(($(string("vDSP_dotpr", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, UInt64),
                  X, 1, Y, 1, val, length(X))
            return val[]
        end
        function distancesq(X::Vector{$T}, Y::Vector{$T})
            val = Ref{$T}(0.0)
            ccall(($(string("vDSP_distancesq", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, UInt64),
                  X, 1, Y, 1, val, length(X))
            return val[]
        end
    end
end

# ============================================================
# vDSP Clipping & thresholding
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    # vclip: clip to [low, high]
    @eval begin
        function vclip!(result::Vector{$T}, X::Vector{$T}, low::$T, high::$T)
            ccall(($(string("vDSP_vclip", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(low), Ref(high), result, 1, length(X))
            return result
        end
        function vclip(X::Vector{$T}, low::$T, high::$T)
            result = similar(X)
            vclip!(result, X, low, high)
        end
    end

    # viclip: inverted clip (pass through values outside [low, high])
    @eval begin
        function viclip!(result::Vector{$T}, X::Vector{$T}, low::$T, high::$T)
            ccall(($(string("vDSP_viclip", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(low), Ref(high), result, 1, length(X))
            return result
        end
        function viclip(X::Vector{$T}, low::$T, high::$T)
            result = similar(X)
            viclip!(result, X, low, high)
        end
    end

    # vthr: threshold — X[n] >= threshold ? X[n] : threshold
    @eval begin
        function vthr!(result::Vector{$T}, X::Vector{$T}, threshold::$T)
            ccall(($(string("vDSP_vthr", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(threshold), result, 1, length(X))
            return result
        end
        function vthr(X::Vector{$T}, threshold::$T)
            result = similar(X)
            vthr!(result, X, threshold)
        end
    end

    # vthres: threshold — X[n] >= threshold ? X[n] : 0
    @eval begin
        function vthres!(result::Vector{$T}, X::Vector{$T}, threshold::$T)
            ccall(($(string("vDSP_vthres", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(threshold), result, 1, length(X))
            return result
        end
        function vthres(X::Vector{$T}, threshold::$T)
            result = similar(X)
            vthres!(result, X, threshold)
        end
    end

    # vcmprs: compress — keep X[n] where gate[n] != 0
    @eval begin
        function vcmprs!(result::Vector{$T}, X::Vector{$T}, gate::Vector{$T})
            ccall(($(string("vDSP_vcmprs", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  X, 1, gate, 1, result, 1, length(X))
            return result
        end
        function vcmprs(X::Vector{$T}, gate::Vector{$T})
            result = similar(X)
            vcmprs!(result, X, gate)
        end
    end
end

# ============================================================
# vDSP Type conversion
# ============================================================
function vdouble(X::Vector{Float32})
    result = Vector{Float64}(undef, length(X))
    ccall(("vDSP_vspdp", libacc), Cvoid,
          (Ptr{Float32}, Int64, Ptr{Float64}, Int64, UInt64),
          X, 1, result, 1, length(X))
    return result
end

function vsingle(X::Vector{Float64})
    result = Vector{Float32}(undef, length(X))
    ccall(("vDSP_vdpsp", libacc), Cvoid,
          (Ptr{Float64}, Int64, Ptr{Float32}, Int64, UInt64),
          X, 1, result, 1, length(X))
    return result
end

# ============================================================
# vDSP Ramp generation
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vramp(start::$T, step::$T, n::Integer)
            result = Vector{$T}(undef, n)
            ccall(($(string("vDSP_vramp", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  Ref(start), Ref(step), result, 1, n)
            return result
        end

        function vrampmul!(result::Vector{$T}, X::Vector{$T}, start::$T, step::$T)
            s = Ref{$T}(start)
            ccall(($(string("vDSP_vrampmul", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, s, Ref(step), result, 1, length(X))
            return result
        end
        function vrampmul(X::Vector{$T}, start::$T, step::$T)
            result = similar(X)
            vrampmul!(result, X, start, step)
        end
    end
end

# ============================================================
# vDSP Integration & running operations
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    # vrsum: running sum scaled by scale
    @eval begin
        function vrsum!(result::Vector{$T}, X::Vector{$T}, scale::$T)
            ccall(($(string("vDSP_vrsum", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(scale), result, 1, length(X))
            return result
        end
        function vrsum(X::Vector{$T}, scale::$T)
            result = similar(X)
            vrsum!(result, X, scale)
        end
    end

    # vsimps: Simpson's rule integration
    @eval begin
        function vsimps!(result::Vector{$T}, X::Vector{$T}, step::$T)
            ccall(($(string("vDSP_vsimps", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(step), result, 1, length(X))
            return result
        end
        function vsimps(X::Vector{$T}, step::$T)
            result = similar(X)
            vsimps!(result, X, step)
        end
    end

    # vtrapz: trapezoidal integration
    @eval begin
        function vtrapz!(result::Vector{$T}, X::Vector{$T}, step::$T)
            ccall(($(string("vDSP_vtrapz", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  X, 1, Ref(step), result, 1, length(X))
            return result
        end
        function vtrapz(X::Vector{$T}, step::$T)
            result = similar(X)
            vtrapz!(result, X, step)
        end
    end

    # vswsum: sliding window sum
    @eval begin
        function vswsum!(result::Vector{$T}, X::Vector{$T}, window::Integer)
            n_out = length(X) - window + 1
            ccall(($(string("vDSP_vswsum", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  X, 1, result, 1, n_out, window)
            return result
        end
        function vswsum(X::Vector{$T}, window::Integer)
            n_out = length(X) - window + 1
            result = Vector{$T}(undef, n_out)
            vswsum!(result, X, window)
        end
    end

    # vswmax: sliding window maximum
    @eval begin
        function vswmax!(result::Vector{$T}, X::Vector{$T}, window::Integer)
            n_out = length(X) - window + 1
            ccall(($(string("vDSP_vswmax", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  X, 1, result, 1, n_out, window)
            return result
        end
        function vswmax(X::Vector{$T}, window::Integer)
            n_out = length(X) - window + 1
            result = Vector{$T}(undef, n_out)
            vswmax!(result, X, window)
        end
    end
end

# ============================================================
# vDSP Interpolation
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    # vintb: interpolate D[n] = A[n] + t*(B[n]-A[n])
    @eval begin
        function vintb!(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, t::$T)
            ccall(($(string("vDSP_vintb", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, B, 1, Ref(t), result, 1, length(A))
            return result
        end
        function vintb(A::Vector{$T}, B::Vector{$T}, t::$T)
            result = similar(A)
            vintb!(result, A, B, t)
        end
    end

    # vlint: linear interpolation lookup
    @eval begin
        function vlint!(result::Vector{$T}, table::Vector{$T}, indices::Vector{$T})
            ccall(($(string("vDSP_vlint", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  table, indices, 1, result, 1, length(indices), length(table))
            return result
        end
        function vlint(table::Vector{$T}, indices::Vector{$T})
            result = Vector{$T}(undef, length(indices))
            vlint!(result, table, indices)
        end
    end

    # vqint: quadratic interpolation lookup
    @eval begin
        function vqint!(result::Vector{$T}, table::Vector{$T}, indices::Vector{$T})
            ccall(($(string("vDSP_vqint", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  table, indices, 1, result, 1, length(indices), length(table))
            return result
        end
        function vqint(table::Vector{$T}, indices::Vector{$T})
            result = Vector{$T}(undef, length(indices))
            vqint!(result, table, indices)
        end
    end
end

# ============================================================
# vDSP Polynomial evaluation
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vpoly!(result::Vector{$T}, coeffs::Vector{$T}, X::Vector{$T})
            ccall(($(string("vDSP_vpoly", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  coeffs, 1, X, 1, result, 1, length(X), length(coeffs) - 1)
            return result
        end
        function vpoly(coeffs::Vector{$T}, X::Vector{$T})
            result = similar(X)
            vpoly!(result, coeffs, X)
        end
    end
end

# ============================================================
# vDSP Normalization
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vnormalize!(result::Vector{$T}, X::Vector{$T})
            mean_out = Ref{$T}(0.0)
            stddev_out = Ref{$T}(0.0)
            ccall(($(string("vDSP_normalize", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, UInt64),
                  X, 1, result, 1, mean_out, stddev_out, length(X))
            return (result, mean_out[], stddev_out[])
        end
        function vnormalize(X::Vector{$T})
            result = similar(X)
            vnormalize!(result, X)
        end
    end
end

# ============================================================
# vDSP Zero crossings
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function nzcros(X::Vector{$T}, max_crossings::Integer=0)
            if max_crossings <= 0
                max_crossings = length(X)
            end
            indices = Vector{UInt64}(undef, max_crossings)
            count = Ref{UInt64}(0)
            ccall(($(string("vDSP_nzcros", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, UInt64, Ptr{UInt64}, Ptr{UInt64}, UInt64),
                  X, 1, max_crossings, indices, count, length(X))
            n = Int(count[])
            return (indices[1:n], n)
        end
    end
end

# ============================================================
# vDSP Decibel conversion
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vdbcon!(result::Vector{$T}, X::Vector{$T}, ref::$T, power::Bool=true)
            flag = power ? UInt32(0) : UInt32(1)
            ccall(($(string("vDSP_vdbcon", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64, UInt32),
                  X, 1, Ref(ref), result, 1, length(X), flag)
            return result
        end
        function vdbcon(X::Vector{$T}, ref::$T, power::Bool=true)
            result = similar(X)
            vdbcon!(result, X, ref, power)
        end
    end
end

