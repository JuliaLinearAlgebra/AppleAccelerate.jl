## array.jl ##

mksymtuple(x) = (x, x)

for (T, suff) in ((Float64, ""), (Float32, "f"))

    # 1-arg functions
    onearg_funcs = (:ceil,:floor,:sqrt,:rsqrt,:rec,
              :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
              :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
              :sinh,:cosh,:tanh,:asinh,:acosh,:atanh)
    # Combine this with the list of renamed 1-arg funcs
    onearg_funcs = (
        map(mksymtuple, onearg_funcs)...,
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

    # 1-arg cube root. Kept in its own loop (rather than added to the shared
    # `onearg_funcs` tuple) because `onearg_funcs` also drives the broadcasting
    # overrides below, which would install `Base.broadcasted(::typeof(cbrt), …)`
    # methods on top of `Base.cbrt` and change global broadcast behaviour.
    for (f, fa) in ((:cbrt,:cbrt),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T})
                out = Array{$T}(undef, size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::Array{$T}, X::Array{$T})
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ptr{$T},Ref{Cint}),out,X,length(X))
                out
            end
        end
    end

    # 2-arg functions where the C prototype takes the operands in the order
    # (z, y, x, n) and computes z = op(y, x). The Julia-facing call is
    # `f(X, Y) = op(X, Y)`, so we pass the C-args as (out, X, Y, n).
    for (f, fa) in ((:nextafter,:nextafter), (:remainder,:remainder))
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

    # vvpows(z, y, x, n): z[i] = pow(x[i], y) for a vector base x and a SCALAR
    # exponent y (passed by reference). Julia-facing: `pows(X, y) = X .^ y`.
    for (f, fa) in ((:pows,:pows),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::Array{$T}, y::$T)
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, y)
            end
            function ($f!)(out::Array{$T}, X::Array{$T}, y::$T)
                ccall(($(string("vv",fa,suff)),libacc),Cvoid,
                      (Ptr{$T},Ref{$T},Ptr{$T},Ref{Cint}),out,y,X,length(X))
                out
            end
        end
    end

    # two-arg return
    for (f, _) in ((:sincos,:sincos),)
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

    for (f, _) in onearg_funcs
        f! = Symbol("$(f)!")
        @eval begin
            # Broadcasting override such that f.(X) turns into f(X)
            Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N}}}) where {Style, Axes, N} = ($f)(bc.args...)
            Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N}}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
        end
        if T == Float32
            @eval begin
                Base.broadcasted(::typeof($f), arg::Array{<:Union{Float32,Float64}}) = ($f)(arg)
                Base.broadcasted(::typeof($f), arg::Base.Broadcast.Broadcasted) = ($f)(copy(arg))
            end
        end
    end
    for (f, _) in (twoarg_funcs...,(:pow,:pow))
        f! = Symbol("$(f)!")
        @eval begin
            # Broadcasting override such that f.(X) turns into f(X)
            Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N},Array{$T,N}}}) where {Style, Axes, N} = ($f)(bc.args...)
            Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T,N},Array{$T,N}}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
        end
        if T == Float32
            @eval begin
                Base.broadcasted(::typeof($f), arg1::Array{F}, arg2::Array{F}) where {F<:Union{Float32,Float64}} = ($f)(arg1, arg2)
                Base.broadcasted(::typeof($f), arg1::Array{<:Union{Float32,Float64}}, arg2::Base.Broadcast.Broadcasted) = ($f)(arg1, copy(arg2))
                Base.broadcasted(::typeof($f), arg1::Base.Broadcast.Broadcasted, arg2::Array{<:Union{Float32,Float64}}) = ($f)(copy(arg1), arg2)
                Base.broadcasted(::typeof($f), arg1::Base.Broadcast.Broadcasted, arg2::Base.Broadcast.Broadcasted) = ($f)(copy(arg1), copy(arg2))
            end
        end
    end
end

"""
    sincos(X::Array{T}) where T <: Union{Float32, Float64}

Compute the sine and cosine of each element simultaneously via vecLib
[`vvsincos`](https://developer.apple.com/documentation/accelerate/vvsincos(_:_:_:_:)).
Returns a tuple `(sin(X), cos(X))` of arrays. Faster than computing `sin` and `cos`
separately since both are produced in a single pass.

The mutating variant `sincos!(out_sin, out_cos, X)` stores results in preallocated arrays.
"""
sincos

"""
    cis(X::Array{T}) where T <: Union{Float32, Float64}

Compute `cos(x) + im*sin(x)` for each element via vecLib
[`vvcosisin`](https://developer.apple.com/documentation/accelerate/vvcosisin(_:_:_:)).
Returns a `Complex{T}` array. Equivalent to `exp.(im .* X)` but faster.

The mutating variant `cis!(out, X)` stores results in a preallocated complex array.
"""
cis

"""
    cbrt(X::Array{T}) where T <: Union{Float32, Float64}

Compute the cube root of each element via vecLib
[`vvcbrt`](https://developer.apple.com/documentation/accelerate/vvcbrt(_:_:_:)),
i.e. `out[i] = cbrt(X[i])`.

The mutating variant `cbrt!(out, X)` stores results in a preallocated array.
"""
cbrt

"""
    nextafter(X::Array{T}, Y::Array{T}) where T <: Union{Float32, Float64}

Compute the next machine-representable floating-point value after each `X[i]` in
the direction of `Y[i]`, via vecLib
[`vvnextafter`](https://developer.apple.com/documentation/accelerate/vvnextafter(_:_:_:_:)),
i.e. `out[i] = nextafter(X[i], Y[i])`.

The C prototype is `vvnextafter(z, y, x, n)` computing `z[i] = nextafter(y[i], x[i])`,
so the Julia call passes the arrays as `(out, X, Y)`.

The mutating variant `nextafter!(out, X, Y)` stores results in a preallocated array.
"""
nextafter

"""
    remainder(X::Array{T}, Y::Array{T}) where T <: Union{Float32, Float64}

Compute the IEEE-754 remainder of each `X[i]` divided by `Y[i]` via vecLib
[`vvremainder`](https://developer.apple.com/documentation/accelerate/vvremainder(_:_:_:_:)),
i.e. `out[i] = X[i] - k*Y[i]` where `k` is the integer nearest `X[i]/Y[i]`
(ties to even). The result satisfies `abs(out[i]) <= abs(Y[i])/2`. This matches
`rem(X[i], Y[i], RoundNearest)`.

The C prototype is `vvremainder(z, y, x, n)` computing `z[i] = y[i] - k*x[i]`,
so the Julia call passes the arrays as `(out, X, Y)`.

The mutating variant `remainder!(out, X, Y)` stores results in a preallocated array.
"""
remainder

"""
    pows(X::Array{T}, y::T) where T <: Union{Float32, Float64}

Raise each element of the vector base `X` to the scalar exponent `y` via vecLib
[`vvpows`](https://developer.apple.com/documentation/accelerate/vvpows(_:_:_:_:)),
i.e. `out[i] = X[i]^y`. Equivalent to `X .^ y`.

The C prototype is `void vvpows(double *z, const double *y, const double *x, const int *n)`
with `z[i] = pow(x[i], y)`, where `y` is a scalar exponent (passed by reference)
and `x` is the vector base. The Julia call passes `(out, X, y)`.

The mutating variant `pows!(out, X, y)` stores results in a preallocated array.
"""
pows

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
                index = Ref{UInt}(0)
                val = Ref{$T}(0.0)
                ccall(($(string("vDSP_", fa, suff), libacc)),  Cvoid,
                      (Ptr{$T}, Int64,  Ref{$T}, Ref{UInt}, UInt64),
                      X, 1, val, index, length(X))
                return (val[], Int(index[])+1)
            end
        end
    end
end

@doc """
    maximum(X::Vector{T}) where T <: Union{Float32, Float64}

Return the maximum value in `X` via vDSP. Equivalent to `Base.maximum(X)`.
Wraps [`vDSP_maxv`](https://developer.apple.com/documentation/accelerate/vdsp_maxv).
""" maximum

@doc """
    minimum(X::Vector{T}) where T <: Union{Float32, Float64}

Return the minimum value in `X` via vDSP. Equivalent to `Base.minimum(X)`.
Wraps [`vDSP_minv`](https://developer.apple.com/documentation/accelerate/vdsp_minv).
""" minimum

@doc """
    sum(X::Vector{T}) where T <: Union{Float32, Float64}

Return the sum of elements in `X` via vDSP. Equivalent to `Base.sum(X)`.
Wraps [`vDSP_sve`](https://developer.apple.com/documentation/accelerate/vdsp_sve).
""" sum

@doc """
    mean(X::Vector{T}) where T <: Union{Float32, Float64}

Return the arithmetic mean of elements in `X` via vDSP.
Wraps [`vDSP_meanv`](https://developer.apple.com/documentation/accelerate/vdsp_meanv).
""" mean

@doc """
    meanmag(X::Vector{T}) where T <: Union{Float32, Float64}

Return the mean of absolute values: `sum(abs.(X)) / length(X)`.
Wraps [`vDSP_meamgv`](https://developer.apple.com/documentation/accelerate/vdsp_meamgv).
""" meanmag

@doc """
    meansqr(X::Vector{T}) where T <: Union{Float32, Float64}

Return the mean of squares: `sum(X.^2) / length(X)`.
Wraps [`vDSP_measqv`](https://developer.apple.com/documentation/accelerate/vdsp_measqv).
""" meansqr

@doc """
    meanssqr(X::Vector{T}) where T <: Union{Float32, Float64}

Return the mean of signed squares: `sum(X .* abs.(X)) / length(X)`.
Wraps [`vDSP_mvessq`](https://developer.apple.com/documentation/accelerate/vdsp_mvessq).
""" meanssqr

@doc """
    summag(X::Vector{T}) where T <: Union{Float32, Float64}

Return the sum of absolute values: `sum(abs.(X))`.
Wraps [`vDSP_svemg`](https://developer.apple.com/documentation/accelerate/vdsp_svemg).
""" summag

@doc """
    sumsqr(X::Vector{T}) where T <: Union{Float32, Float64}

Return the sum of squares: `sum(X.^2)`.
Wraps [`vDSP_svesq`](https://developer.apple.com/documentation/accelerate/vdsp_svesq).
""" sumsqr

@doc """
    sumssqr(X::Vector{T}) where T <: Union{Float32, Float64}

Return the sum of signed squares: `sum(X .* abs.(X))`.
Wraps [`vDSP_svs`](https://developer.apple.com/documentation/accelerate/vdsp_svs).
""" sumssqr

@doc """
    findmax(X::Vector{T}) where T <: Union{Float32, Float64}

Return `(value, index)` of the maximum element in `X` via vDSP. Equivalent to `Base.findmax(X)`.
Wraps [`vDSP_maxvi`](https://developer.apple.com/documentation/accelerate/vdsp_maxvi).
""" findmax

@doc """
    findmin(X::Vector{T}) where T <: Union{Float32, Float64}

Return `(value, index)` of the minimum element in `X` via vDSP. Equivalent to `Base.findmin(X)`.
Wraps [`vDSP_minvi`](https://developer.apple.com/documentation/accelerate/vdsp_minvi).
""" findmin

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
            Base.broadcasted(::typeof($f), arg1::Array{$T}, arg2::Array{$T}) = ($f)(arg1, arg2)
            Base.broadcasted(::typeof($f), arg1::Array{$T}, arg2::Base.Broadcast.Broadcasted) = ($f)(arg1, copy(arg2))
            Base.broadcasted(::typeof($f), arg1::Base.Broadcast.Broadcasted, arg2::Array{$T}) = ($f)(copy(arg1), arg2)
        end
    end
end

# Disambiguating methods for (Broadcasted, Broadcasted) — resolves ambiguity between
# Float32 and Float64 broadcasted overrides when both arguments are Broadcasted objects.
for f in (:vadd, :vsub, :vmul, :vdiv)
    @eval Base.broadcasted(::typeof($f), arg1::Base.Broadcast.Broadcasted, arg2::Base.Broadcast.Broadcasted) = ($f)(copy(arg1), copy(arg2))
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
            # c - X == X * (-1) + c, computed in one pass via vDSP_vsmsa (D = A*B + C)
            ccall(($(string("vDSP_vsmsa", suff), libacc)),  Cvoid,
                    (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64,  UInt64),
                    X, 1, Ref(-one($T)), Ref(c), result, 1, length(result))
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
            Base.broadcasted(::typeof($f), arg1::Array{$T}, arg2::$T) = ($f)(arg1, arg2)
            Base.broadcasted(::typeof($f), arg1::Base.Broadcast.Broadcasted, arg2::$T) = ($f)(copy(arg1), arg2)
        end
    end

    f = :svsub
    f! = Symbol("$(f)!")

    @eval begin
        # Broadcasting override such that f.(X) turns into f(X)
        Base.copy(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T, N}, $T}}) where {Style, Axes, N} = ($f)(bc.args...)
        Base.copyto!(dest::Array{$T, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof($f), Tuple{Array{$T,N}, $T}}) where {Style, Axes, N} = ($f!)(dest, bc.args...)
        Base.broadcasted(::typeof($f), arg1::Array{$T}, arg2::$T) = ($f)(arg1, arg2)
        Base.broadcasted(::typeof($f), arg1::Base.Broadcast.Broadcasted, arg2::$T) = ($f)(copy(arg1), arg2)
    end
end

# ============================================================
# vDSP Unary vector operations
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc "Negate each element: `result[i] = -X[i]`. Wraps [`vDSP_vneg`](https://developer.apple.com/documentation/accelerate/vdsp_vneg)." vneg
@doc "Negative absolute value: `result[i] = -|X[i]|`. Wraps [`vDSP_vnabs`](https://developer.apple.com/documentation/accelerate/vdsp_vnabs)." vnabs
@doc "Square each element: `result[i] = X[i]^2`. Wraps [`vDSP_vsq`](https://developer.apple.com/documentation/accelerate/vdsp_vsq)." vsq
@doc "Signed square: `result[i] = X[i] * |X[i]|`. Wraps [`vDSP_vssq`](https://developer.apple.com/documentation/accelerate/vdsp_vssq)." vssq
@doc "Fractional part: `result[i] = X[i] - trunc(X[i])`. Wraps [`vDSP_vfrac`](https://developer.apple.com/documentation/accelerate/vdsp_vfrac)." vfrac
@doc "Absolute value: `result[i] = |X[i]|`. Wraps [`vDSP_vabs`](https://developer.apple.com/documentation/accelerate/vdsp_vabs)." vabs
@doc "Reverse `X` in-place. Wraps [`vDSP_vrvrs`](https://developer.apple.com/documentation/accelerate/vdsp_vrvrs)." vreverse!
@doc "Return a reversed copy of `X`. Wraps [`vDSP_vrvrs`](https://developer.apple.com/documentation/accelerate/vdsp_vrvrs)." vreverse

# ============================================================
# vDSP Two-vector element-wise operations
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc "Element-wise maximum: `result[i] = max(X[i], Y[i])`. Wraps [`vDSP_vmax`](https://developer.apple.com/documentation/accelerate/vdsp_vmax)." vmax
@doc "Element-wise minimum: `result[i] = min(X[i], Y[i])`. Wraps [`vDSP_vmin`](https://developer.apple.com/documentation/accelerate/vdsp_vmin)." vmin
@doc "Element-wise maximum magnitude: `result[i] = max(|X[i]|, |Y[i]|)`. Wraps [`vDSP_vmaxmg`](https://developer.apple.com/documentation/accelerate/vdsp_vmaxmg)." vmaxmg
@doc "Element-wise minimum magnitude: `result[i] = min(|X[i]|, |Y[i]|)`. Wraps [`vDSP_vminmg`](https://developer.apple.com/documentation/accelerate/vdsp_vminmg)." vminmg
@doc "Element-wise Euclidean distance: `result[i] = hypot(X[i], Y[i])`. Wraps [`vDSP_vdist`](https://developer.apple.com/documentation/accelerate/vdsp_vdist)." vdist
@doc "Tapered merge of two vectors. Wraps [`vDSP_vtmerg`](https://developer.apple.com/documentation/accelerate/vdsp_vtmerg)." vtmerg

# ============================================================
# vDSP Scalar-vector divide: C[n] = A / B[n]
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc "Scalar divided by vector: `result[i] = c / X[i]`. Wraps [`vDSP_svdiv`](https://developer.apple.com/documentation/accelerate/vdsp_svdiv)." svdiv

# ============================================================
# vDSP Compound arithmetic
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))

    # 3-vector ops: (A, B, C → D)
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

    # 4-vector ops: (A, B, C, D → E)
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

    # 2 vectors + 1 scalar → D
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

@doc "Vector add and multiply: `result[i] = (A[i] + B[i]) * C[i]`. Wraps [`vDSP_vam`](https://developer.apple.com/documentation/accelerate/vdsp_vam)." vam
@doc "Vector subtract and multiply: `result[i] = (A[i] - B[i]) * C[i]`. Wraps [`vDSP_vsbm`](https://developer.apple.com/documentation/accelerate/vdsp_vsbm)." vsbm
@doc "Signal envelope. Wraps [`vDSP_venvlp`](https://developer.apple.com/documentation/accelerate/vdsp_venvlp)." venvlp
@doc "Vector add, add, and multiply: `result[i] = (A[i] + B[i]) * (C[i] + D[i])`. Wraps [`vDSP_vaam`](https://developer.apple.com/documentation/accelerate/vdsp_vaam)." vaam
@doc "Vector subtract, subtract, and multiply: `result[i] = (A[i] - B[i]) * (C[i] - D[i])`. Wraps [`vDSP_vsbsbm`](https://developer.apple.com/documentation/accelerate/vdsp_vsbsbm)." vsbsbm
@doc "Vector add, subtract, and multiply: `result[i] = (A[i] + B[i]) * (C[i] - D[i])`. Wraps [`vDSP_vasbm`](https://developer.apple.com/documentation/accelerate/vdsp_vasbm)." vasbm
@doc "Pythagorean distance: `result[i] = sqrt((A[i]-C[i])^2 + (B[i]-D[i])^2)`. Wraps [`vDSP_vpythg`](https://developer.apple.com/documentation/accelerate/vdsp_vpythg)." vpythg
@doc "Vector add and scalar multiply: `result[i] = (A[i] + B[i]) * c`. Wraps [`vDSP_vasm`](https://developer.apple.com/documentation/accelerate/vdsp_vasm)." vasm
@doc "Vector subtract and scalar multiply: `result[i] = (A[i] - B[i]) * c`. Wraps [`vDSP_vsbsm`](https://developer.apple.com/documentation/accelerate/vdsp_vsbsm)." vsbsm
@doc "Vector scalar multiply and add: `result[i] = A[i] * b + C[i]`. Wraps [`vDSP_vsma`](https://developer.apple.com/documentation/accelerate/vdsp_vsma)." vsma
@doc "Vector scalar multiply and scalar add: `result[i] = A[i] * b + c`. Wraps [`vDSP_vsmsa`](https://developer.apple.com/documentation/accelerate/vdsp_vsmsa)." vsmsa
@doc "Simultaneous add and subtract: returns `(A .+ B, B .- A)`. Wraps [`vDSP_vaddsub`](https://developer.apple.com/documentation/accelerate/vdsp_vaddsub)." vaddsub

# --- Batch 1: additional compound arithmetic ---
for (T, suff) in ((Float32, ""), (Float64, "D"))

    # vma: D[n] = A[n]*B[n] + C[n]  (3-vector → 1-vector, same pattern as vam)
    # vmsb: D[n] = A[n]*B[n] - C[n]
    for (f, fa) in ((:vma, :vma), (:vmsb, :vmsb))
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

    # vmma: E[n] = A[n]*B[n] + C[n]*D[n]  (4-vector → 1-vector)
    # vmmsb: E[n] = A[n]*B[n] - C[n]*D[n]
    for (f, fa) in ((:vmma, :vmma), (:vmmsb, :vmmsb))
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

    # vmsa: D[n] = A[n]*B[n] + c  (2-vector + scalar → 1-vector)
    @eval begin
        function vmsa!(result::Vector{$T}, A::Vector{$T}, B::Vector{$T}, c::$T)
            ccall(($(string("vDSP_vmsa", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, B, 1, Ref(c), result, 1, length(A))
            return result
        end
        function vmsa(A::Vector{$T}, B::Vector{$T}, c::$T)
            result = similar(A)
            vmsa!(result, A, B, c)
        end
    end

    # vsmsb: D[n] = A[n]*b - C[n]  (vector*scalar - vector)
    @eval begin
        function vsmsb!(result::Vector{$T}, A::Vector{$T}, b::$T, C::Vector{$T})
            ccall(($(string("vDSP_vsmsb", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(b), C, 1, result, 1, length(A))
            return result
        end
        function vsmsb(A::Vector{$T}, b::$T, C::Vector{$T})
            result = similar(A)
            vsmsb!(result, A, b, C)
        end
    end

    # vsmsma: E[n] = A[n]*b + C[n]*d  (vector*scalar + vector*scalar)
    @eval begin
        function vsmsma!(result::Vector{$T}, A::Vector{$T}, b::$T, C::Vector{$T}, d::$T)
            ccall(($(string("vDSP_vsmsma", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(b), C, 1, Ref(d), result, 1, length(A))
            return result
        end
        function vsmsma(A::Vector{$T}, b::$T, C::Vector{$T}, d::$T)
            result = similar(A)
            vsmsma!(result, A, b, C, d)
        end
    end
end

@doc "Vector multiply and add: `result[i] = A[i]*B[i] + C[i]`. Wraps [`vDSP_vma`](https://developer.apple.com/documentation/accelerate/vdsp_vma)." vma
@doc "Vector multiply and subtract: `result[i] = A[i]*B[i] - C[i]`. Wraps [`vDSP_vmsb`](https://developer.apple.com/documentation/accelerate/vdsp_vmsb)." vmsb
@doc "Vector multiply, multiply, and add: `result[i] = A[i]*B[i] + C[i]*D[i]`. Wraps [`vDSP_vmma`](https://developer.apple.com/documentation/accelerate/vdsp_vmma)." vmma
@doc "Vector multiply, multiply, and subtract: `result[i] = A[i]*B[i] - C[i]*D[i]`. Wraps [`vDSP_vmmsb`](https://developer.apple.com/documentation/accelerate/vdsp_vmmsb)." vmmsb
@doc "Vector multiply and scalar add: `result[i] = A[i]*B[i] + c`. Wraps [`vDSP_vmsa`](https://developer.apple.com/documentation/accelerate/vdsp_vmsa)." vmsa
@doc "Vector scalar multiply and subtract: `result[i] = A[i]*b - C[i]`. Wraps [`vDSP_vsmsb`](https://developer.apple.com/documentation/accelerate/vdsp_vsmsb)." vsmsb
@doc "Vector scalar multiply and scalar multiply add: `result[i] = A[i]*b + C[i]*d`. Wraps [`vDSP_vsmsma`](https://developer.apple.com/documentation/accelerate/vdsp_vsmsma)." vsmsma

# ============================================================
# vDSP Scalar reductions: dot product, distance squared
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc "Dot product: `sum(X .* Y)`. Wraps [`vDSP_dotpr`](https://developer.apple.com/documentation/accelerate/vdsp_dotpr)." dot
@doc "Squared Euclidean distance: `sum((X .- Y).^2)`. Wraps [`vDSP_distancesq`](https://developer.apple.com/documentation/accelerate/vdsp_distancesq)." distancesq

# --- Batch 2: additional reductions ---
for (T, suff) in ((Float32, ""), (Float64, "D"))

    # rmsqv: root mean square — sqrt(sum(X.^2)/N)
    @eval begin
        function rmsqv(X::Vector{$T})
            val = Ref{$T}(0)
            ccall(($(string("vDSP_rmsqv", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, UInt64),
                  X, 1, val, length(X))
            return val[]
        end
    end

    # sve_svesq: simultaneous sum and sum-of-squares
    @eval begin
        function sve_svesq(X::Vector{$T})
            s = Ref{$T}(0)
            ssq = Ref{$T}(0)
            ccall(($(string("vDSP_sve_svesq", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, UInt64),
                  X, 1, s, ssq, length(X))
            return (s[], ssq[])
        end
    end

    # maxmgv / minmgv: max/min magnitude
    for (f, fa) in ((:maxmgv, :maxmgv), (:minmgv, :minmgv))
        @eval begin
            function ($f)(X::Vector{$T})
                val = Ref{$T}(0)
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, UInt64),
                      X, 1, val, length(X))
                return val[]
            end
        end
    end

    # maxmgvi / minmgvi: max/min magnitude with index
    for (f, fa) in ((:maxmgvi, :maxmgvi), (:minmgvi, :minmgvi))
        @eval begin
            function ($f)(X::Vector{$T})
                val = Ref{$T}(0)
                idx = Ref{UInt}(0)
                ccall(($(string("vDSP_", fa, suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Ptr{UInt}, UInt64),
                      X, 1, val, idx, length(X))
                return (val[], Int(idx[]) + 1)
            end
        end
    end
end

@doc "Root mean square: `sqrt(sum(X.^2) / length(X))`. Wraps [`vDSP_rmsqv`](https://developer.apple.com/documentation/accelerate/vdsp_rmsqv)." rmsqv
@doc "Simultaneous sum and sum-of-squares: returns `(sum(X), sum(X.^2))`. Wraps [`vDSP_sve_svesq`](https://developer.apple.com/documentation/accelerate/vdsp_sve_svesq)." sve_svesq
@doc "Maximum magnitude: `maximum(abs.(X))`. Wraps [`vDSP_maxmgv`](https://developer.apple.com/documentation/accelerate/vdsp_maxmgv)." maxmgv
@doc "Minimum magnitude: `minimum(abs.(X))`. Wraps [`vDSP_minmgv`](https://developer.apple.com/documentation/accelerate/vdsp_minmgv)." minmgv
@doc "Maximum magnitude with index: returns `(maximum(abs.(X)), index)`. Wraps [`vDSP_maxmgvi`](https://developer.apple.com/documentation/accelerate/vdsp_maxmgvi)." maxmgvi
@doc "Minimum magnitude with index: returns `(minimum(abs.(X)), index)`. Wraps [`vDSP_minmgvi`](https://developer.apple.com/documentation/accelerate/vdsp_minmgvi)." minmgvi

# ============================================================
# vDSP Clipping & thresholding
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
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

@doc "Clip values to `[low, high]`: `clamp.(X, low, high)`. Wraps [`vDSP_vclip`](https://developer.apple.com/documentation/accelerate/vdsp_vclip)." vclip
@doc "Inverted clip: pass through values outside `[low, high]`, zero inside. Wraps [`vDSP_viclip`](https://developer.apple.com/documentation/accelerate/vdsp_viclip)." viclip
@doc "Threshold: `result[i] = X[i] >= threshold ? X[i] : threshold`. Wraps [`vDSP_vthr`](https://developer.apple.com/documentation/accelerate/vdsp_vthr)." vthr
@doc "Threshold to zero: `result[i] = X[i] >= threshold ? X[i] : 0`. Wraps [`vDSP_vthres`](https://developer.apple.com/documentation/accelerate/vdsp_vthres)." vthres
@doc "Compress: gather elements of `X` where `gate` is nonzero. Wraps [`vDSP_vcmprs`](https://developer.apple.com/documentation/accelerate/vdsp_vcmprs)." vcmprs

# ============================================================
# vDSP Type conversion
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================

"""
    vdouble(X::Vector{Float32}) -> Vector{Float64}

Convert single-precision to double-precision. Wraps [`vDSP_vspdp`](https://developer.apple.com/documentation/accelerate/vdsp_vspdp).
"""
function vdouble(X::Vector{Float32})
    result = Vector{Float64}(undef, length(X))
    ccall(("vDSP_vspdp", libacc), Cvoid,
          (Ptr{Float32}, Int64, Ptr{Float64}, Int64, UInt64),
          X, 1, result, 1, length(X))
    return result
end

"""
    vsingle(X::Vector{Float64}) -> Vector{Float32}

Convert double-precision to single-precision. Wraps [`vDSP_vdpsp`](https://developer.apple.com/documentation/accelerate/vdsp_vdpsp).
"""
function vsingle(X::Vector{Float64})
    result = Vector{Float32}(undef, length(X))
    ccall(("vDSP_vdpsp", libacc), Cvoid,
          (Ptr{Float64}, Int64, Ptr{Float32}, Int64, UInt64),
          X, 1, result, 1, length(X))
    return result
end

# ============================================================
# vDSP Ramp generation
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc "Generate a ramp: `result[i] = start + i * step` for `i = 0, ..., n-1`. Wraps [`vDSP_vramp`](https://developer.apple.com/documentation/accelerate/vdsp_vramp)." vramp
@doc "Multiply vector by a generated ramp. Wraps [`vDSP_vrampmul`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmul)." vrampmul

for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vrampmul2!(O0::Vector{$T}, O1::Vector{$T}, I0::Vector{$T}, I1::Vector{$T}, start::$T, step::$T)
            n = length(I0)
            s = Ref{$T}(start)
            ccall(($(string("vDSP_vrampmul2", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  I0, I1, 1, s, Ref(step), O0, O1, 1, n)
            return (O0, O1)
        end
        function vrampmul2(I0::Vector{$T}, I1::Vector{$T}, start::$T, step::$T)
            O0 = similar(I0)
            O1 = similar(I1)
            vrampmul2!(O0, O1, I0, I1, start, step)
        end
    end
end

@doc "Stereo ramp multiply: multiply two vectors by the same ramp. Wraps [`vDSP_vrampmul2`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmul2)." vrampmul2

for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vavlin!(C::Vector{$T}, A::Vector{$T}, weight::$T)
            ccall(($(string("vDSP_vavlin", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(weight), C, 1, length(A))
            return C
        end
        function vavlin(A::Vector{$T}, C::Vector{$T}, weight::$T)
            result = copy(C)
            vavlin!(result, A, weight)
        end
    end
end

@doc "Vector linear average: `C[n] = (C[n] * weight + A[n]) / (weight + 1)`. Wraps [`vDSP_vavlin`](https://developer.apple.com/documentation/accelerate/vdsp_vavlin)." vavlin

# ============================================================
# vDSP Integration & running operations
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
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

@doc "Running sum scaled by `scale`. Wraps [`vDSP_vrsum`](https://developer.apple.com/documentation/accelerate/vdsp_vrsum)." vrsum
@doc "Simpson's rule integration with step size `step`. Wraps [`vDSP_vsimps`](https://developer.apple.com/documentation/accelerate/vdsp_vsimps)." vsimps
@doc "Trapezoidal integration with step size `step`. Wraps [`vDSP_vtrapz`](https://developer.apple.com/documentation/accelerate/vdsp_vtrapz)." vtrapz
@doc "Sliding window sum with window size `window`. Returns a vector of length `length(X) - window + 1`. Wraps [`vDSP_vswsum`](https://developer.apple.com/documentation/accelerate/vdsp_vswsum)." vswsum
@doc "Sliding window maximum with window size `window`. Returns a vector of length `length(X) - window + 1`. Wraps [`vDSP_vswmax`](https://developer.apple.com/documentation/accelerate/vdsp_vswmax)." vswmax

# ============================================================
# vDSP Interpolation
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
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

@doc "Vector interpolation: `result[i] = A[i] + t * (B[i] - A[i])`. Wraps [`vDSP_vintb`](https://developer.apple.com/documentation/accelerate/vdsp_vintb)." vintb
@doc "Linear interpolation from a lookup table using fractional `indices`. Wraps [`vDSP_vlint`](https://developer.apple.com/documentation/accelerate/vdsp_vlint)." vlint
@doc "Quadratic interpolation from a lookup table using fractional `indices`. Wraps [`vDSP_vqint`](https://developer.apple.com/documentation/accelerate/vdsp_vqint)." vqint

# ============================================================
# vDSP Polynomial evaluation
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc """
    vpoly(coeffs, X)

Evaluate polynomial at each point in `X`. Coefficients are highest degree first:
`[a_P, a_{P-1}, ..., a_1, a_0]`. Wraps [`vDSP_vpoly`](https://developer.apple.com/documentation/accelerate/vdsp_vpoly).
""" vpoly

# ============================================================
# vDSP Normalization
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc """
    vnormalize(X) -> (normalized, mean, stddev)

Normalize vector to zero mean and unit standard deviation: `(X .- mean) ./ stddev`.
Returns a tuple of `(normalized_vector, mean, stddev)`.
Wraps [`vDSP_normalize`](https://developer.apple.com/documentation/accelerate/vdsp_normalize).
""" vnormalize

# ============================================================
# vDSP Zero crossings
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc """
    nzcros(X, max_crossings=0) -> (indices, count)

Find zero crossings in `X`. Returns a tuple of `(crossing_indices, count)`.
If `max_crossings <= 0`, searches for up to `length(X)` crossings.
Wraps [`vDSP_nzcros`](https://developer.apple.com/documentation/accelerate/vdsp_nzcros).
""" nzcros

# ============================================================
# vDSP Decibel conversion
# See: https://developer.apple.com/documentation/accelerate/vdsp
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

@doc """
    vdbcon(X, ref, power=true)

Convert to decibels relative to `ref`. If `power=true`, computes `10*log10(X/ref)`;
if `power=false`, computes `20*log10(X/ref)`.
Wraps [`vDSP_vdbcon`](https://developer.apple.com/documentation/accelerate/vdsp_vdbcon).
""" vdbcon

# ============================================================
# Batch 3: Vector Utility — fill, swap, gather, generation,
#   clipping variants, sorting, table lookup, integer ops,
#   format conversion
# ============================================================

# --- Data fill / clear / swap ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vclr!(C::Vector{$T})
            ccall(($(string("vDSP_vclr", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, UInt64),
                  C, 1, length(C))
            return C
        end
        function vfill!(C::Vector{$T}, a::$T)
            ccall(($(string("vDSP_vfill", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  Ref(a), C, 1, length(C))
            return C
        end
        function vswap!(A::Vector{$T}, B::Vector{$T})
            ccall(($(string("vDSP_vswap", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, 1, B, 1, length(A))
            return (A, B)
        end
    end
end

@doc "Fill vector with zeros: `C[i] = 0`. Wraps [`vDSP_vclr`](https://developer.apple.com/documentation/accelerate/vdsp_vclr)." vclr!
@doc "Fill vector with scalar value: `C[i] = a`. Wraps [`vDSP_vfill`](https://developer.apple.com/documentation/accelerate/vdsp_vfill)." vfill!
@doc "Swap elements of two vectors in-place. Wraps [`vDSP_vswap`](https://developer.apple.com/documentation/accelerate/vdsp_vswap)." vswap!

# --- Gathering / indexing ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vgathr!(C::Vector{$T}, A::Vector{$T}, B::Vector{UInt})
            ccall(($(string("vDSP_vgathr", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{UInt}, Int64, Ptr{$T}, Int64, UInt64),
                  A, B, 1, C, 1, length(B))
            return C
        end
        function vgathr(A::Vector{$T}, B::Vector{UInt})
            C = Vector{$T}(undef, length(B))
            vgathr!(C, A, B)
        end
        function vindex!(C::Vector{$T}, A::Vector{$T}, B::Vector{$T})
            ccall(($(string("vDSP_vindex", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  A, B, 1, C, 1, length(B))
            return C
        end
        function vindex(A::Vector{$T}, B::Vector{$T})
            C = Vector{$T}(undef, length(B))
            vindex!(C, A, B)
        end
    end
end

@doc "Gather elements by index: `C[i] = A[B[i]]` where B contains 1-based UInt indices. Wraps [`vDSP_vgathr`](https://developer.apple.com/documentation/accelerate/vdsp_vgathr)." vgathr
@doc "Index with float indices: `C[i] = A[trunc(B[i])]` where B contains 0-based float indices. Wraps [`vDSP_vindex`](https://developer.apple.com/documentation/accelerate/vdsp_vindex)." vindex

# --- Generation ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vgen!(C::Vector{$T}, a::$T, b::$T)
            ccall(($(string("vDSP_vgen", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  Ref(a), Ref(b), C, 1, length(C))
            return C
        end
        function vgen(a::$T, b::$T, n::Integer)
            C = Vector{$T}(undef, n)
            vgen!(C, a, b)
        end
        function vgenp!(C::Vector{$T}, A::Vector{$T}, B::Vector{$T}, n::Integer)
            ccall(($(string("vDSP_vgenp", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  A, 1, B, 1, C, 1, n, length(A))
            return C
        end
        function vgenp(A::Vector{$T}, B::Vector{$T}, n::Integer)
            C = Vector{$T}(undef, n)
            vgenp!(C, A, B, n)
        end
    end
end

@doc "Generate linear ramp between two scalars: `C[i] = a + (b-a)*i/(N-1)`. Wraps [`vDSP_vgen`](https://developer.apple.com/documentation/accelerate/vdsp_vgen)." vgen
@doc "Piecewise linear interpolation from breakpoints. Wraps [`vDSP_vgenp`](https://developer.apple.com/documentation/accelerate/vdsp_vgenp)." vgenp

# --- Clipping / thresholding variants ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vclipc!(result::Vector{$T}, X::Vector{$T}, low::$T, high::$T)
            nlow = Ref{UInt64}(0)
            nhigh = Ref{UInt64}(0)
            ccall(($(string("vDSP_vclipc", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64, Ptr{UInt64}, Ptr{UInt64}),
                  X, 1, Ref(low), Ref(high), result, 1, length(X), nlow, nhigh)
            return (result, Int(nlow[]), Int(nhigh[]))
        end
        function vclipc(X::Vector{$T}, low::$T, high::$T)
            result = similar(X)
            vclipc!(result, X, low, high)
        end
        function vlim!(result::Vector{$T}, A::Vector{$T}, b::$T, c::$T)
            ccall(($(string("vDSP_vlim", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(b), Ref(c), result, 1, length(A))
            return result
        end
        function vlim(A::Vector{$T}, b::$T, c::$T)
            result = similar(A)
            vlim!(result, A, b, c)
        end
        function vthrsc!(result::Vector{$T}, A::Vector{$T}, b::$T, c::$T)
            ccall(($(string("vDSP_vthrsc", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(b), Ref(c), result, 1, length(A))
            return result
        end
        function vthrsc(A::Vector{$T}, b::$T, c::$T)
            result = similar(A)
            vthrsc!(result, A, b, c)
        end
    end
end

@doc "Clip with count: returns `(clipped, nlow, nhigh)`. Wraps [`vDSP_vclipc`](https://developer.apple.com/documentation/accelerate/vdsp_vclipc)." vclipc
@doc "Test limit: `result[i] = (b <= A[i]) ? c : -c`. Wraps [`vDSP_vlim`](https://developer.apple.com/documentation/accelerate/vdsp_vlim)." vlim
@doc "Threshold with signed constant. Wraps [`vDSP_vthrsc`](https://developer.apple.com/documentation/accelerate/vdsp_vthrsc)." vthrsc

# --- Sorting ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vsort!(X::Vector{$T}, ascending::Bool=true)
            order = ascending ? Cint(1) : Cint(-1)
            ccall(($(string("vDSP_vsort", suff)), libacc), Cvoid,
                  (Ptr{$T}, UInt64, Cint),
                  X, length(X), order)
            return X
        end
        function vsorti!(indices::Vector{UInt}, X::Vector{$T}, ascending::Bool=true)
            order = ascending ? Cint(1) : Cint(-1)
            ccall(($(string("vDSP_vsorti", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{UInt}, Ptr{Cvoid}, UInt64, Cint),
                  X, indices, C_NULL, length(X), order)
            return indices
        end
        function vsorti(X::Vector{$T}, ascending::Bool=true)
            indices = Vector{UInt}(undef, length(X))
            for i in 1:length(X)
                indices[i] = UInt(i - 1)
            end
            vsorti!(indices, X, ascending)
            return Int.(indices) .+ 1
        end
    end
end

@doc "Sort vector in-place. `ascending=true` for ascending order. Wraps [`vDSP_vsort`](https://developer.apple.com/documentation/accelerate/vdsp_vsort)." vsort!
@doc "Return sort permutation (1-based indices). Wraps [`vDSP_vsorti`](https://developer.apple.com/documentation/accelerate/vdsp_vsorti)." vsorti

@doc """
    vsorti!(indices::Vector{UInt}, X::Vector{T}, ascending::Bool=true)

In-place index sort. `indices` MUST be pre-filled with the 0-based identity
permutation `0:length(X)-1` before calling; `vDSP_vsorti` reorders those existing
indices rather than generating them. Passing an uninitialized (`undef`) or otherwise
non-identity buffer yields a silently wrong permutation. Use the allocating
[`vsorti`](@ref) if you want the initialization handled for you. Wraps
[`vDSP_vsorti`](https://developer.apple.com/documentation/accelerate/vdsp_vsorti).
""" vsorti!

# --- Table lookup ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vtabi!(D::Vector{$T}, A::Vector{$T}, s1::$T, s2::$T, C::Vector{$T})
            ccall(($(string("vDSP_vtabi", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Ptr{$T}, Ptr{$T}, UInt64, Ptr{$T}, Int64, UInt64),
                  A, 1, Ref(s1), Ref(s2), C, length(C), D, 1, length(A))
            return D
        end
        function vtabi(A::Vector{$T}, s1::$T, s2::$T, C::Vector{$T})
            D = Vector{$T}(undef, length(A))
            vtabi!(D, A, s1, s2, C)
        end
    end
end

@doc "Table lookup with interpolation: `D[i] = C[clamp(s1*A[i]+s2, 0, M-1)]`. Wraps [`vDSP_vtabi`](https://developer.apple.com/documentation/accelerate/vdsp_vtabi)." vtabi

# --- Integer operations (Int32) ---
function vaddi!(C::Vector{Int32}, A::Vector{Int32}, B::Vector{Int32})
    ccall(("vDSP_vaddi", libacc), Cvoid,
          (Ptr{Int32}, Int64, Ptr{Int32}, Int64, Ptr{Int32}, Int64, UInt64),
          A, 1, B, 1, C, 1, length(A))
    return C
end
function vaddi(A::Vector{Int32}, B::Vector{Int32})
    C = similar(A)
    vaddi!(C, A, B)
end

function vabsi!(C::Vector{Int32}, A::Vector{Int32})
    ccall(("vDSP_vabsi", libacc), Cvoid,
          (Ptr{Int32}, Int64, Ptr{Int32}, Int64, UInt64),
          A, 1, C, 1, length(A))
    return C
end
function vabsi(A::Vector{Int32})
    C = similar(A)
    vabsi!(C, A)
end

function vfilli!(C::Vector{Int32}, a::Int32)
    ccall(("vDSP_vfilli", libacc), Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Int64, UInt64),
          Ref(a), C, 1, length(C))
    return C
end

function veqvi!(C::Vector{Int32}, A::Vector{Int32}, B::Vector{Int32})
    ccall(("vDSP_veqvi", libacc), Cvoid,
          (Ptr{Int32}, Int64, Ptr{Int32}, Int64, Ptr{Int32}, Int64, UInt64),
          A, 1, B, 1, C, 1, length(A))
    return C
end
function veqvi(A::Vector{Int32}, B::Vector{Int32})
    C = similar(A)
    veqvi!(C, A, B)
end

@doc "Int32 vector addition: `C[i] = A[i] + B[i]`. Wraps [`vDSP_vaddi`](https://developer.apple.com/documentation/accelerate/vdsp_vaddi)." vaddi
@doc "Int32 absolute value: `C[i] = |A[i]|`. Wraps [`vDSP_vabsi`](https://developer.apple.com/documentation/accelerate/vdsp_vabsi)." vabsi
@doc "Fill Int32 vector with scalar value. Wraps [`vDSP_vfilli`](https://developer.apple.com/documentation/accelerate/vdsp_vfilli)." vfilli!
@doc "Int32 bitwise XNOR: `C[i] = ~(A[i] ^ B[i])`. Wraps [`vDSP_veqvi`](https://developer.apple.com/documentation/accelerate/vdsp_veqvi)." veqvi

# ============================================================
# Batch 4: Matrix Operations
# Note: vDSP uses row-major layout, Julia uses column-major.
# We swap dimensions so the user passes Julia matrices naturally.
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function mmul!(C::Matrix{$T}, A::Matrix{$T}, B::Matrix{$T})
            # Julia (col-major) A is m×p, B is p×n → C is m×n
            # vDSP sees transposed layouts, so we compute Bᵀ × Aᵀ = (AB)ᵀ
            m, p = size(A)
            p2, n = size(B)
            p == p2 || throw(DimensionMismatch("A columns ($p) ≠ B rows ($p2)"))
            size(C) == (m, n) || throw(DimensionMismatch("C must be $m×$n"))
            ccall(($(string("vDSP_mmul", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64, UInt64),
                  B, 1, A, 1, C, 1, UInt64(n), UInt64(m), UInt64(p))
            return C
        end
        function mmul(A::Matrix{$T}, B::Matrix{$T})
            m = size(A, 1)
            n = size(B, 2)
            C = Matrix{$T}(undef, m, n)
            mmul!(C, A, B)
        end
    end

    @eval begin
        function mtrans!(C::Matrix{$T}, A::Matrix{$T})
            m, n = size(A)
            size(C) == (n, m) || throw(DimensionMismatch("C must be $n×$m"))
            # vDSP_mtrans(__A, __IA, __C, __IC, __M, __N):
            #   __M = columns in A (rows in C), __N = rows in A (columns in C)
            # Julia col-major m×n viewed as row-major: n rows × m cols
            # We want C to be row-major m rows × n cols, so __M = m, __N = n
            ccall(($(string("vDSP_mtrans", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64, UInt64),
                  A, 1, C, 1, UInt64(m), UInt64(n))
            return C
        end
        function mtrans(A::Matrix{$T})
            m, n = size(A)
            C = Matrix{$T}(undef, n, m)
            mtrans!(C, A)
        end
    end

    @eval begin
        function mmov!(C::Matrix{$T}, A::Matrix{$T})
            m, n = size(A)
            ta = UInt64(m)  # column stride of source (= number of rows in col-major)
            tc = UInt64(size(C, 1))  # column stride of destination
            ccall(($(string("vDSP_mmov", suff)), libacc), Cvoid,
                  (Ptr{$T}, Ptr{$T}, UInt64, UInt64, UInt64, UInt64),
                  A, C, UInt64(m), UInt64(n), ta, tc)
            return C
        end
        function mmov(A::Matrix{$T})
            C = similar(A)
            mmov!(C, A)
        end
    end
end

@doc "Matrix multiply: `C = A * B`. Wraps [`vDSP_mmul`](https://developer.apple.com/documentation/accelerate/vdsp_mmul)." mmul
@doc "Matrix transpose: `C = Aᵀ`. Wraps [`vDSP_mtrans`](https://developer.apple.com/documentation/accelerate/vdsp_mtrans)." mtrans
@doc "Matrix copy (submatrix move). Wraps [`vDSP_mmov`](https://developer.apple.com/documentation/accelerate/vdsp_mmov)." mmov

# ============================================================
# Batch 6: Type Conversion (int ↔ float)
# ============================================================

# float → signed int (truncating)
for (T, suff) in ((Float32, ""), (Float64, "D"))
    for (intT, intname) in ((Int8, "8"), (Int16, "16"), (Int32, "32"))
        fname = Symbol("vfix$(intname)")
        fname! = Symbol("vfix$(intname)!")
        vdsp_name = string("vDSP_vfix", intname, suff)
        @eval begin
            function ($fname!)(C::Vector{$intT}, A::Vector{$T})
                ccall(($vdsp_name, libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$intT}, Int64, UInt64),
                      A, 1, C, 1, length(A))
                return C
            end
            function ($fname)(A::Vector{$T})
                C = Vector{$intT}(undef, length(A))
                ($fname!)(C, A)
            end
        end
    end
end

# float → unsigned int (truncating)
for (T, suff) in ((Float32, ""), (Float64, "D"))
    for (intT, intname) in ((UInt8, "8"), (UInt16, "16"), (UInt32, "32"))
        fname = Symbol("vfixu$(intname)")
        fname! = Symbol("vfixu$(intname)!")
        vdsp_name = string("vDSP_vfixu", intname, suff)
        @eval begin
            function ($fname!)(C::Vector{$intT}, A::Vector{$T})
                ccall(($vdsp_name, libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$intT}, Int64, UInt64),
                      A, 1, C, 1, length(A))
                return C
            end
            function ($fname)(A::Vector{$T})
                C = Vector{$intT}(undef, length(A))
                ($fname!)(C, A)
            end
        end
    end
end

# float → signed int (rounding)
for (T, suff) in ((Float32, ""), (Float64, "D"))
    for (intT, intname) in ((Int8, "8"), (Int16, "16"), (Int32, "32"))
        fname = Symbol("vfixr$(intname)")
        fname! = Symbol("vfixr$(intname)!")
        vdsp_name = string("vDSP_vfixr", intname, suff)
        @eval begin
            function ($fname!)(C::Vector{$intT}, A::Vector{$T})
                ccall(($vdsp_name, libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$intT}, Int64, UInt64),
                      A, 1, C, 1, length(A))
                return C
            end
            function ($fname)(A::Vector{$T})
                C = Vector{$intT}(undef, length(A))
                ($fname!)(C, A)
            end
        end
    end
end

# float → unsigned int (rounding)
for (T, suff) in ((Float32, ""), (Float64, "D"))
    for (intT, intname) in ((UInt8, "8"), (UInt16, "16"), (UInt32, "32"))
        fname = Symbol("vfixru$(intname)")
        fname! = Symbol("vfixru$(intname)!")
        vdsp_name = string("vDSP_vfixru", intname, suff)
        @eval begin
            function ($fname!)(C::Vector{$intT}, A::Vector{$T})
                ccall(($vdsp_name, libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$intT}, Int64, UInt64),
                      A, 1, C, 1, length(A))
                return C
            end
            function ($fname)(A::Vector{$T})
                C = Vector{$intT}(undef, length(A))
                ($fname!)(C, A)
            end
        end
    end
end

# signed int → float
for (T, suff) in ((Float32, ""), (Float64, "D"))
    for (intT, intname) in ((Int8, "8"), (Int16, "16"), (Int32, "32"))
        fname = Symbol("vflt$(intname)")
        fname! = Symbol("vflt$(intname)!")
        vdsp_name = string("vDSP_vflt", intname, suff)
        @eval begin
            function ($fname!)(C::Vector{$T}, A::Vector{$intT})
                ccall(($vdsp_name, libacc), Cvoid,
                      (Ptr{$intT}, Int64, Ptr{$T}, Int64, UInt64),
                      A, 1, C, 1, length(A))
                return C
            end
            function ($fname)(A::Vector{$intT}, ::Type{$T})
                C = Vector{$T}(undef, length(A))
                ($fname!)(C, A)
            end
        end
    end
end

# unsigned int → float
for (T, suff) in ((Float32, ""), (Float64, "D"))
    for (intT, intname) in ((UInt8, "8"), (UInt16, "16"), (UInt32, "32"))
        fname = Symbol("vfltu$(intname)")
        fname! = Symbol("vfltu$(intname)!")
        vdsp_name = string("vDSP_vfltu", intname, suff)
        @eval begin
            function ($fname!)(C::Vector{$T}, A::Vector{$intT})
                ccall(($vdsp_name, libacc), Cvoid,
                      (Ptr{$intT}, Int64, Ptr{$T}, Int64, UInt64),
                      A, 1, C, 1, length(A))
                return C
            end
            function ($fname)(A::Vector{$intT}, ::Type{$T})
                C = Vector{$T}(undef, length(A))
                ($fname!)(C, A)
            end
        end
    end
end

# ============================================================
# Batch 8: Image Convolution (2D)
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function f3x3!(C::Matrix{$T}, A::Matrix{$T}, F::Matrix{$T})
            nr, nc = size(A)
            size(F) == (3, 3) || throw(DimensionMismatch("Filter must be 3×3"))
            size(C) == (nr, nc) || throw(DimensionMismatch("C must match A size"))
            # vDSP interprets data in row-major order; Julia matrices are
            # column-major, so a (nr, nc) Julia matrix is laid out as an
            # (nc, nr) row-major image. Swap the dimensions so the convolution
            # operates on the intended matrix (and the result reads back
            # correctly as a column-major (nr, nc) matrix).
            ccall(($(string("vDSP_f3x3", suff)), libacc), Cvoid,
                  (Ptr{$T}, UInt64, UInt64, Ptr{$T}, Ptr{$T}),
                  A, UInt64(nc), UInt64(nr), F, C)
            return C
        end
        function f3x3(A::Matrix{$T}, F::Matrix{$T})
            C = similar(A)
            f3x3!(C, A, F)
        end
        function f5x5!(C::Matrix{$T}, A::Matrix{$T}, F::Matrix{$T})
            nr, nc = size(A)
            size(F) == (5, 5) || throw(DimensionMismatch("Filter must be 5×5"))
            size(C) == (nr, nc) || throw(DimensionMismatch("C must match A size"))
            # See note in `f3x3!`: swap dims for row-major vs column-major.
            ccall(($(string("vDSP_f5x5", suff)), libacc), Cvoid,
                  (Ptr{$T}, UInt64, UInt64, Ptr{$T}, Ptr{$T}),
                  A, UInt64(nc), UInt64(nr), F, C)
            return C
        end
        function f5x5(A::Matrix{$T}, F::Matrix{$T})
            C = similar(A)
            f5x5!(C, A, F)
        end
        function imgfir!(C::Matrix{$T}, A::Matrix{$T}, F::Matrix{$T})
            nr, nc = size(A)
            fr, fc = size(F)
            size(C) == (nr, nc) || throw(DimensionMismatch("C must match A size"))
            # See note in `f3x3!`: vDSP reads row-major while Julia is
            # column-major, so swap both the image and filter dimensions.
            ccall(($(string("vDSP_imgfir", suff)), libacc), Cvoid,
                  (Ptr{$T}, UInt64, UInt64, Ptr{$T}, Ptr{$T}, UInt64, UInt64),
                  A, UInt64(nc), UInt64(nr), F, C, UInt64(fc), UInt64(fr))
            return C
        end
        function imgfir(A::Matrix{$T}, F::Matrix{$T})
            C = similar(A)
            imgfir!(C, A, F)
        end
    end
end

@doc "2D convolution with a 3×3 filter. Border elements are set to zero. Wraps [`vDSP_f3x3`](https://developer.apple.com/documentation/accelerate/vdsp_f3x3)." f3x3
@doc "2D convolution with a 5×5 filter. Border elements are set to zero. Wraps [`vDSP_f5x5`](https://developer.apple.com/documentation/accelerate/vdsp_f5x5)." f5x5
@doc "General 2D image convolution with a P×Q filter. Border elements are set to zero. Wraps [`vDSP_imgfir`](https://developer.apple.com/documentation/accelerate/vdsp_imgfir)." imgfir
