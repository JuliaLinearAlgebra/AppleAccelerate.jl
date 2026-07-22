## array.jl ##

# ------------------------------------------------------------
# StridedArray support helpers
#
# vDSP functions take explicit stride arguments, so most wrappers below accept
# any StridedArray/StridedVector and pass `stride(X, 1)` through to vDSP.
# vecLib vv* functions (the transcendental math wrappers) take NO stride
# argument and require contiguous memory, as do a few vDSP entry points whose
# prototypes have no stride parameter for a given operand (e.g. the table in
# `vlint`) or whose stride semantics are unclear (e.g. sliding windows).
# ------------------------------------------------------------

@inline function _iscontiguous(A::StridedArray)
    expected = 1
    for i in 1:ndims(A)
        size(A, i) == 1 || stride(A, i) == expected || return false
        expected *= size(A, i)
    end
    return true
end
@inline _iscontiguous(A::DenseArray) = true

function _check_contiguous(As::StridedArray...)
    for A in As
        _iscontiguous(A) || throw(ArgumentError(
            "this operation requires contiguous (unit-stride) arrays, but got an " *
            "array with strides $(strides(A)); copy the data first (e.g. with `collect`)"))
    end
    return nothing
end

@inline function _check_unit_stride(A::StridedVector, name::Symbol)
    stride(A, 1) == 1 || throw(ArgumentError(
        "$name requires a unit-stride (contiguous) vector for this argument, but got " *
        "stride $(stride(A, 1)); copy the data first (e.g. with `collect`)"))
    return nothing
end

# vDSP index-returning reductions report the 0-based *memory offset* of the
# element (a multiple of the stride), so a positive stride is required to map
# it back to a 1-based vector index.
@inline function _check_positive_stride(A::StridedVector, name::Symbol)
    stride(A, 1) >= 1 || throw(ArgumentError(
        "$name requires a positive stride, but got stride $(stride(A, 1))"))
    return nothing
end

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
            function ($f)(X::StridedArray{$T})
                out = Array{$T}(undef, size(X))
                ($f!)(out, X)
            end

            # In-place mutating variant
            function ($f!)(out::StridedArray{$T}, X::StridedArray{$T})
                _check_contiguous(out, X)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, X, Ref{Cint}(length(X)))
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
            function ($f)(X::StridedArray{$T}, Y::StridedArray{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::StridedArray{$T}, X::StridedArray{$T}, Y::StridedArray{$T})
                _check_contiguous(out, X, Y)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, X, Y, Ref{Cint}(length(X)))
                out
            end
        end
    end

    # for some bizarre reason, vvpow/vvpowf reverse the order of arguments.
    for (f, fa) in ((:pow,:pow),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::StridedArray{$T}, Y::StridedArray{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::StridedArray{$T}, X::StridedArray{$T}, Y::StridedArray{$T})
                _check_contiguous(out, X, Y)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, Y, X, Ref{Cint}(length(X)))
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
            function ($f)(X::StridedArray{$T})
                out = Array{$T}(undef, size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::StridedArray{$T}, X::StridedArray{$T})
                _check_contiguous(out, X)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, X, Ref{Cint}(length(X)))
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
            function ($f)(X::StridedArray{$T}, Y::StridedArray{$T})
                size(X) == size(Y) || throw(DimensionMismatch("Arguments must have same shape"))
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, Y)
            end
            function ($f!)(out::StridedArray{$T}, X::StridedArray{$T}, Y::StridedArray{$T})
                _check_contiguous(out, X, Y)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, X, Y, Ref{Cint}(length(X)))
                out
            end
        end
    end

    # vvpows(z, y, x, n): z[i] = pow(x[i], y) for a vector base x and a SCALAR
    # exponent y (passed by reference). Julia-facing: `pows(X, y) = X .^ y`.
    for (f, fa) in ((:pows,:pows),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::StridedArray{$T}, y::$T)
                out = Array{$T}(undef, size(X))
                ($f!)(out, X, y)
            end
            function ($f!)(out::StridedArray{$T}, X::StridedArray{$T}, y::$T)
                _check_contiguous(out, X)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, Ref(y), X, Ref{Cint}(length(X)))
                out
            end
        end
    end

    # two-arg return
    for (f, _) in ((:sincos,:sincos),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::StridedArray{$T})
                out1 = Array{$T}(undef, size(X))
                out2 = Array{$T}(undef, size(X))
                ($f!)(out1, out2, X)
            end
            function ($f!)(out1::StridedArray{$T}, out2::StridedArray{$T}, X::StridedArray{$T})
                _check_contiguous(out1, out2, X)
                LibAccelerate.$(Symbol(string("vv",f,suff)))(out1, out2, X, Ref{Cint}(length(X)))
                out1, out2
            end
        end
    end

    # complex return
    for (f,fa) in ((:cis,:cosisin),)
        f! = Symbol("$(f)!")
        @eval begin
            function ($f)(X::StridedArray{$T})
                out = Array{Complex{$T}}(undef, size(X))
                ($f!)(out, X)
            end
            function ($f!)(out::StridedArray{Complex{$T}}, X::StridedArray{$T})
                _check_contiguous(out, X)
                LibAccelerate.$(Symbol(string("vv",fa,suff)))(out, X, Ref{Cint}(length(X)))
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
        # vDSP returns ±Inf for max/min of an empty vector; match Base and throw instead.
        # The mean family divides by length and returns NaN for an empty vector;
        # match Base.mean and throw. Sum-type reductions correctly return 0 for
        # an empty vector and are left unguarded.
        guard = if f in (:maximum, :minimum)
            :(isempty(X) && throw(ArgumentError($(string(f)) * " over an empty collection is not allowed")))
        elseif f in (:mean, :meanmag, :meansqr, :meanssqr)
            :(isempty(X) && throw(ArgumentError("mean of empty collection")))
        else
            :(nothing)
        end
        @eval begin
            function ($f)(X::StridedVector{$T})
                $guard
                val = Ref{$T}(0.0)
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(X,stride(X,1), val, length(X))
                return val[]
            end
        end
    end

    for (f, fa) in ((:findmax, :maxvi), (:findmin, :minvi))
        @eval begin
            function ($f)(X::StridedVector{$T})
                # vDSP returns (±Inf, 0) for an empty vector; match Base and throw instead.
                isempty(X) && throw(ArgumentError($(string(f)) * " over an empty collection is not allowed"))
                _check_positive_stride(X, $(QuoteNode(f)))
                index = Ref{UInt}(0)
                val = Ref{$T}(0.0)
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(X,stride(X,1), val, index, length(X))
                # vDSP reports the 0-based memory offset (a multiple of the stride)
                return (val[], div(Int(index[]), stride(X,1)) + 1)
            end
        end
    end
end

@doc """
    maximum(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the maximum value in `X` via vDSP. Equivalent to `Base.maximum(X)`.
Wraps [`vDSP_maxv`](https://developer.apple.com/documentation/accelerate/vdsp_maxv).
""" maximum

@doc """
    minimum(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the minimum value in `X` via vDSP. Equivalent to `Base.minimum(X)`.
Wraps [`vDSP_minv`](https://developer.apple.com/documentation/accelerate/vdsp_minv).
""" minimum

@doc """
    sum(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the sum of elements in `X` via vDSP. Equivalent to `Base.sum(X)`.
Wraps [`vDSP_sve`](https://developer.apple.com/documentation/accelerate/vdsp_sve).
""" sum

@doc """
    mean(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the arithmetic mean of elements in `X` via vDSP.
Wraps [`vDSP_meanv`](https://developer.apple.com/documentation/accelerate/vdsp_meanv).
""" mean

@doc """
    meanmag(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the mean of absolute values: `sum(abs.(X)) / length(X)`.
Wraps [`vDSP_meamgv`](https://developer.apple.com/documentation/accelerate/vdsp_meamgv).
""" meanmag

@doc """
    meansqr(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the mean of squares: `sum(X.^2) / length(X)`.
Wraps [`vDSP_measqv`](https://developer.apple.com/documentation/accelerate/vdsp_measqv).
""" meansqr

@doc """
    meanssqr(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the mean of signed squares: `sum(X .* abs.(X)) / length(X)`.
Wraps [`vDSP_mvessq`](https://developer.apple.com/documentation/accelerate/vdsp_mvessq).
""" meanssqr

@doc """
    summag(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the sum of absolute values: `sum(abs.(X))`.
Wraps [`vDSP_svemg`](https://developer.apple.com/documentation/accelerate/vdsp_svemg).
""" summag

@doc """
    sumsqr(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the sum of squares: `sum(X.^2)`.
Wraps [`vDSP_svesq`](https://developer.apple.com/documentation/accelerate/vdsp_svesq).
""" sumsqr

@doc """
    sumssqr(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return the sum of signed squares: `sum(X .* abs.(X))`.
Wraps [`vDSP_svs`](https://developer.apple.com/documentation/accelerate/vdsp_svs).
""" sumssqr

@doc """
    findmax(X::StridedVector{T}) where T <: Union{Float32, Float64}

Return `(value, index)` of the maximum element in `X` via vDSP. Equivalent to `Base.findmax(X)`.
Wraps [`vDSP_maxvi`](https://developer.apple.com/documentation/accelerate/vdsp_maxvi).
""" findmax

@doc """
    findmin(X::StridedVector{T}) where T <: Union{Float32, Float64}

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
            `$($f!)(result::StridedVector{$($T)}, X::StridedVector{$($T)}, Y::StridedVector{$($T)})`

            Implements element-wise **$($name)** over two **Vector{$($T)}** and overwrites
            the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
            """ ->
            function ($f!)(result::StridedVector{$T}, X::StridedVector{$T}, Y::StridedVector{$T})
                (length(X) == length(Y) == length(result)) ||
                    throw(DimensionMismatch("result, X and Y must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", f, suff)))(Y,stride(Y,1),X,stride(X,1),result,stride(result,1), length(result))
                return result
            end
        end

        @eval begin
            @doc """
            `$($f)(X::StridedVector{$($T)}, Y::StridedVector{$($T)})`

            Implements element-wise **$($name)** over two **Vector{$($T)}**. Allocates
            memory to store result. *Returns:* **Vector{$($T)}**
            """ ->
            function ($f)(X::StridedVector{$T}, Y::StridedVector{$T})
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
            `$($f!)(result::StridedVector{$($T)}, X::StridedVector{$($T)}, c::$($T))`

            Implements vector-scalar **$($name)** over **Vector{$($T)}** and $($T) and overwrites
            the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
            """ ->
            function ($f!)(result::StridedVector{$T}, X::StridedVector{$T}, c::$T)
                length(result) == length(X) ||
                    throw(DimensionMismatch("result and X must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", f, suff)))(X,stride(X,1), Ref(c),result,stride(result,1), length(result))
                return result
            end
        end

        @eval begin
            @doc """
            `$($f)(X::StridedVector{$($T)}, c::$($T))`

            Implements vector-scalar **$($name)** over **Vector{$($T)}** and $($T). Allocates
            memory to store result. *Returns:* **Vector{$($T)}**
            """ ->
            function ($f)(X::StridedVector{$T}, c::$T)
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
        `$($f!)(result::StridedVector{$($T)}, X::StridedVector{$($T)}, c::$($T))`

        Implements vector-scalar **subtraction** over **Vector{$($T)}** and $($T) and overwrites
        the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
        """ ->
        function ($f!)(result::StridedVector{$T}, X::StridedVector{$T}, c::$T)
            length(result) == length(X) ||
                throw(DimensionMismatch("result and X must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vsadd", suff)))(X,stride(X,1), Ref(-c),result,stride(result,1), length(result))
            return result
        end
    end

    @eval begin
        @doc """
        `$($f)(X::StridedVector{$($T)}, c::$($T))`

        Implements vector-scalar **subtraction** over **Vector{$($T)}** and $($T). Allocates
        memory to store result. *Returns:* **Vector{$($T)}**
        """ ->
        function ($f)(X::StridedVector{$T}, c::$T)
            result = similar(X)
            ($f!)(result, X, c)
            return result
        end
    end

    f = :svsub
    f! = Symbol("$(f)!")

    @eval begin
        @doc """
        `$($f!)(result::StridedVector{$($T)}, X::StridedVector{$($T)}, c::$($T))`

        Implements vector-scalar **subtraction** over $($T) and **Vector{$($T)}** and overwrites
        the result vector with computed value. *Returns:* **Vector{$($T)}** `result`
        """ ->
        function ($f!)(result::StridedVector{$T}, X::StridedVector{$T}, c::$T)
            length(result) == length(X) ||
                throw(DimensionMismatch("result and X must have the same length"))
            # c - X == X * (-1) + c, computed in one pass via vDSP_vsmsa (D = A*B + C)
            LibAccelerate.$(Symbol(string("vDSP_vsmsa", suff)))(X,stride(X,1), Ref(-one($T)), Ref(c),result,stride(result,1), length(result))
            return result
        end
    end

    @eval begin
        @doc """
        `$($f)(X::StridedVector{$($T), c::$($T)})`

        Implements vector-scalar **subtraction** over $($T) and  **Vector{$($T)}**. Allocates
        memory to store result. *Returns:* **Vector{$($T)}**
        """ ->
        function ($f)(X::StridedVector{$T}, c::$T)
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
            function ($f!)(result::StridedVector{$T}, X::StridedVector{$T})
                length(result) >= length(X) ||
                    throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(X,stride(X,1),result,stride(result,1),length(X))
                return result
            end
            function ($f)(X::StridedVector{$T})
                result = similar(X)
                ($f!)(result, X)
            end
        end
    end

    # vreverse: in-place only
    @eval begin
        function vreverse!(X::StridedVector{$T})
            LibAccelerate.$(Symbol(string("vDSP_vrvrs", suff)))(X,stride(X,1),length(X))
            return X
        end
        function vreverse(X::StridedVector{$T})
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
            function ($f!)(result::StridedVector{$T}, X::StridedVector{$T}, Y::StridedVector{$T})
                (length(X) == length(Y) == length(result)) ||
                    throw(DimensionMismatch("result, X and Y must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(X,stride(X,1),Y,stride(Y,1),result,stride(result,1),length(X))
                return result
            end
            function ($f)(X::StridedVector{$T}, Y::StridedVector{$T})
                result = similar(X)
                ($f!)(result, X, Y)
            end
        end
    end

    @eval begin
        function vtmerg!(result::StridedVector{$T}, X::StridedVector{$T}, Y::StridedVector{$T})
            (length(X) == length(Y) == length(result)) ||
                throw(DimensionMismatch("result, X and Y must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vtmerg", suff)))(X,stride(X,1),Y,stride(Y,1),result,stride(result,1),length(X))
            return result
        end
        function vtmerg(X::StridedVector{$T}, Y::StridedVector{$T})
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
        function svdiv!(result::StridedVector{$T}, X::StridedVector{$T}, c::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_svdiv", suff)))(Ref(c),X,stride(X,1),result,stride(result,1),length(X))
            return result
        end
        function svdiv(X::StridedVector{$T}, c::$T)
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
            function ($f!)(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T})
                (length(A) == length(B) == length(C) == length(result)) ||
                    throw(DimensionMismatch("result, A, B and C must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),result,stride(result,1),length(A))
                return result
            end
            function ($f)(A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T})
                result = similar(A)
                ($f!)(result, A, B, C)
            end
        end
    end

    @eval begin
        function venvlp!(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T})
            (length(A) == length(B) == length(C) == length(result)) ||
                throw(DimensionMismatch("result, A, B and C must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_venvlp", suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),result,stride(result,1),length(A))
            return result
        end
        function venvlp(A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T})
            result = similar(A)
            venvlp!(result, A, B, C)
        end
    end

    # 4-vector ops: (A, B, C, D → E)
    for (f, fa) in ((:vaam, :vaam), (:vsbsbm, :vsbsbm), (:vasbm, :vasbm))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T}, D::StridedVector{$T})
                (length(A) == length(B) == length(C) == length(D) == length(result)) ||
                    throw(DimensionMismatch("result, A, B, C and D must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),D,stride(D,1),result,stride(result,1),length(A))
                return result
            end
            function ($f)(A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T}, D::StridedVector{$T})
                result = similar(A)
                ($f!)(result, A, B, C, D)
            end
        end
    end

    @eval begin
        function vpythg!(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T}, D::StridedVector{$T})
            (length(A) == length(B) == length(C) == length(D) == length(result)) ||
                throw(DimensionMismatch("result, A, B, C and D must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vpythg", suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),D,stride(D,1),result,stride(result,1),length(A))
            return result
        end
        function vpythg(A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T}, D::StridedVector{$T})
            result = similar(A)
            vpythg!(result, A, B, C, D)
        end
    end

    # 2 vectors + 1 scalar → D
    for (f, fa) in ((:vasm, :vasm), (:vsbsm, :vsbsm))
        f! = Symbol("$(f)!")
        @eval begin
            function ($f!)(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, c::$T)
                (length(A) == length(B) == length(result)) ||
                    throw(DimensionMismatch("result, A and B must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(A,stride(A,1),B,stride(B,1),Ref(c),result,stride(result,1),length(A))
                return result
            end
            function ($f)(A::StridedVector{$T}, B::StridedVector{$T}, c::$T)
                result = similar(A)
                ($f!)(result, A, B, c)
            end
        end
    end

    @eval begin
        function vsma!(result::StridedVector{$T}, A::StridedVector{$T}, b::$T, C::StridedVector{$T})
            (length(A) == length(C) == length(result)) ||
                throw(DimensionMismatch("result, A and C must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vsma", suff)))(A,stride(A,1),Ref(b),C,stride(C,1),result,stride(result,1),length(A))
            return result
        end
        function vsma(A::StridedVector{$T}, b::$T, C::StridedVector{$T})
            result = similar(A)
            vsma!(result, A, b, C)
        end
    end

    @eval begin
        function vsmsa!(result::StridedVector{$T}, A::StridedVector{$T}, b::$T, c::$T)
            length(result) >= length(A) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(A) ($(length(A)))"))
            LibAccelerate.$(Symbol(string("vDSP_vsmsa", suff)))(A,stride(A,1),Ref(b),Ref(c),result,stride(result,1),length(A))
            return result
        end
        function vsmsa(A::StridedVector{$T}, b::$T, c::$T)
            result = similar(A)
            vsmsa!(result, A, b, c)
        end
    end

    @eval begin
        function vaddsub!(add_result::StridedVector{$T}, sub_result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T})
            (length(A) == length(B) == length(add_result) == length(sub_result)) ||
                throw(DimensionMismatch("add_result, sub_result, A and B must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vaddsub", suff)))(A,stride(A,1),B,stride(B,1),add_result,stride(add_result,1),sub_result,stride(sub_result,1),length(A))
            return (add_result, sub_result)
        end
        function vaddsub(A::StridedVector{$T}, B::StridedVector{$T})
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
            function ($f!)(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T})
                (length(A) == length(B) == length(C) == length(result)) ||
                    throw(DimensionMismatch("result, A, B and C must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),result,stride(result,1),length(A))
                return result
            end
            function ($f)(A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T})
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
            function ($f!)(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T}, D::StridedVector{$T})
                (length(A) == length(B) == length(C) == length(D) == length(result)) ||
                    throw(DimensionMismatch("result, A, B, C and D must have the same length"))
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),D,stride(D,1),result,stride(result,1),length(A))
                return result
            end
            function ($f)(A::StridedVector{$T}, B::StridedVector{$T}, C::StridedVector{$T}, D::StridedVector{$T})
                result = similar(A)
                ($f!)(result, A, B, C, D)
            end
        end
    end

    # vmsa: D[n] = A[n]*B[n] + c  (2-vector + scalar → 1-vector)
    @eval begin
        function vmsa!(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, c::$T)
            (length(A) == length(B) == length(result)) ||
                throw(DimensionMismatch("result, A and B must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vmsa", suff)))(A,stride(A,1),B,stride(B,1),Ref(c),result,stride(result,1),length(A))
            return result
        end
        function vmsa(A::StridedVector{$T}, B::StridedVector{$T}, c::$T)
            result = similar(A)
            vmsa!(result, A, B, c)
        end
    end

    # vsmsb: D[n] = A[n]*b - C[n]  (vector*scalar - vector)
    @eval begin
        function vsmsb!(result::StridedVector{$T}, A::StridedVector{$T}, b::$T, C::StridedVector{$T})
            (length(A) == length(C) == length(result)) ||
                throw(DimensionMismatch("result, A and C must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vsmsb", suff)))(A,stride(A,1),Ref(b),C,stride(C,1),result,stride(result,1),length(A))
            return result
        end
        function vsmsb(A::StridedVector{$T}, b::$T, C::StridedVector{$T})
            result = similar(A)
            vsmsb!(result, A, b, C)
        end
    end

    # vsmsma: E[n] = A[n]*b + C[n]*d  (vector*scalar + vector*scalar)
    @eval begin
        function vsmsma!(result::StridedVector{$T}, A::StridedVector{$T}, b::$T, C::StridedVector{$T}, d::$T)
            (length(A) == length(C) == length(result)) ||
                throw(DimensionMismatch("result, A and C must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vsmsma", suff)))(A,stride(A,1),Ref(b),C,stride(C,1),Ref(d),result,stride(result,1),length(A))
            return result
        end
        function vsmsma(A::StridedVector{$T}, b::$T, C::StridedVector{$T}, d::$T)
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
        function dot(X::StridedVector{$T}, Y::StridedVector{$T})
            length(X) == length(Y) ||
                throw(DimensionMismatch("X and Y must have the same length"))
            val = Ref{$T}(0.0)
            LibAccelerate.$(Symbol(string("vDSP_dotpr", suff)))(X,stride(X,1),Y,stride(Y,1),val,length(X))
            return val[]
        end
        function distancesq(X::StridedVector{$T}, Y::StridedVector{$T})
            length(X) == length(Y) ||
                throw(DimensionMismatch("X and Y must have the same length"))
            val = Ref{$T}(0.0)
            LibAccelerate.$(Symbol(string("vDSP_distancesq", suff)))(X,stride(X,1),Y,stride(Y,1),val,length(X))
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
        function rmsqv(X::StridedVector{$T})
            val = Ref{$T}(0)
            LibAccelerate.$(Symbol(string("vDSP_rmsqv", suff)))(X,stride(X,1),val,length(X))
            return val[]
        end
    end

    # sve_svesq: simultaneous sum and sum-of-squares
    @eval begin
        function sve_svesq(X::StridedVector{$T})
            s = Ref{$T}(0)
            ssq = Ref{$T}(0)
            LibAccelerate.$(Symbol(string("vDSP_sve_svesq", suff)))(X,stride(X,1),s,ssq,length(X))
            return (s[], ssq[])
        end
    end

    # maxmgv / minmgv: max/min magnitude
    for (f, fa) in ((:maxmgv, :maxmgv), (:minmgv, :minmgv))
        @eval begin
            function ($f)(X::StridedVector{$T})
                val = Ref{$T}(0)
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(X,stride(X,1),val,length(X))
                return val[]
            end
        end
    end

    # maxmgvi / minmgvi: max/min magnitude with index
    for (f, fa) in ((:maxmgvi, :maxmgvi), (:minmgvi, :minmgvi))
        @eval begin
            function ($f)(X::StridedVector{$T})
                _check_positive_stride(X, $(QuoteNode(f)))
                val = Ref{$T}(0)
                idx = Ref{UInt}(0)
                LibAccelerate.$(Symbol(string("vDSP_", fa, suff)))(X,stride(X,1),val,idx,length(X))
                # vDSP reports the 0-based memory offset (a multiple of the stride)
                return (val[], div(Int(idx[]), stride(X,1)) + 1)
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
        function vclip!(result::StridedVector{$T}, X::StridedVector{$T}, low::$T, high::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vclip", suff)))(X,stride(X,1),Ref(low),Ref(high),result,stride(result,1),length(X))
            return result
        end
        function vclip(X::StridedVector{$T}, low::$T, high::$T)
            result = similar(X)
            vclip!(result, X, low, high)
        end
        function viclip!(result::StridedVector{$T}, X::StridedVector{$T}, low::$T, high::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_viclip", suff)))(X,stride(X,1),Ref(low),Ref(high),result,stride(result,1),length(X))
            return result
        end
        function viclip(X::StridedVector{$T}, low::$T, high::$T)
            result = similar(X)
            viclip!(result, X, low, high)
        end
        function vthr!(result::StridedVector{$T}, X::StridedVector{$T}, threshold::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vthr", suff)))(X,stride(X,1),Ref(threshold),result,stride(result,1),length(X))
            return result
        end
        function vthr(X::StridedVector{$T}, threshold::$T)
            result = similar(X)
            vthr!(result, X, threshold)
        end
        function vthres!(result::StridedVector{$T}, X::StridedVector{$T}, threshold::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vthres", suff)))(X,stride(X,1),Ref(threshold),result,stride(result,1),length(X))
            return result
        end
        function vthres(X::StridedVector{$T}, threshold::$T)
            result = similar(X)
            vthres!(result, X, threshold)
        end
        function vcmprs!(result::StridedVector{$T}, X::StridedVector{$T}, gate::StridedVector{$T})
            length(gate) >= length(X) ||
                throw(DimensionMismatch("gate length ($(length(gate))) must be at least length(X) ($(length(X)))"))
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vcmprs", suff)))(X,stride(X,1),gate,stride(gate,1),result,stride(result,1),length(X))
            return result
        end
        function vcmprs(X::StridedVector{$T}, gate::StridedVector{$T})
            # vDSP_vcmprs packs the gated elements into the front of the buffer; the
            # meaningful result length is the number of nonzero gate entries, so
            # write into a full-length scratch buffer and trim to that count.
            result = similar(X)
            vcmprs!(result, X, gate)
            count = Base.count(!iszero, @view gate[1:length(X)])
            return result[1:count]
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
    vdouble(X::StridedVector{Float32}) -> Vector{Float64}

Convert single-precision to double-precision. Wraps [`vDSP_vspdp`](https://developer.apple.com/documentation/accelerate/vdsp_vspdp).
"""
function vdouble(X::StridedVector{Float32})
    result = Vector{Float64}(undef, length(X))
    LibAccelerate.vDSP_vspdp(X,stride(X,1),result,stride(result,1),length(X))
    return result
end

"""
    vsingle(X::StridedVector{Float64}) -> Vector{Float32}

Convert double-precision to single-precision. Wraps [`vDSP_vdpsp`](https://developer.apple.com/documentation/accelerate/vdsp_vdpsp).
"""
function vsingle(X::StridedVector{Float64})
    result = Vector{Float32}(undef, length(X))
    LibAccelerate.vDSP_vdpsp(X,stride(X,1),result,stride(result,1),length(X))
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
            LibAccelerate.$(Symbol(string("vDSP_vramp", suff)))(Ref(start),Ref(step),result,stride(result,1),n)
            return result
        end
        function vrampmul!(result::StridedVector{$T}, X::StridedVector{$T}, start::$T, step::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            s = Ref{$T}(start)
            LibAccelerate.$(Symbol(string("vDSP_vrampmul", suff)))(X,stride(X,1),s,Ref(step),result,stride(result,1),length(X))
            return result
        end
        function vrampmul(X::StridedVector{$T}, start::$T, step::$T)
            result = similar(X)
            vrampmul!(result, X, start, step)
        end
    end
end

@doc "Generate a ramp: `result[i] = start + i * step` for `i = 0, ..., n-1`. Wraps [`vDSP_vramp`](https://developer.apple.com/documentation/accelerate/vdsp_vramp)." vramp
@doc "Multiply vector by a generated ramp. Wraps [`vDSP_vrampmul`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmul)." vrampmul

for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vrampmul2!(O0::StridedVector{$T}, O1::StridedVector{$T}, I0::StridedVector{$T}, I1::StridedVector{$T}, start::$T, step::$T)
            n = length(I0)
            (length(I1) >= n && length(O0) >= n && length(O1) >= n) ||
                throw(DimensionMismatch("I1, O0 and O1 must each have length at least length(I0) ($n)"))
            stride(I0,1) == stride(I1,1) ||
                throw(ArgumentError("vrampmul2!: I0 and I1 must have the same stride (vDSP shares one input stride)"))
            stride(O0,1) == stride(O1,1) ||
                throw(ArgumentError("vrampmul2!: O0 and O1 must have the same stride (vDSP shares one output stride)"))
            s = Ref{$T}(start)
            LibAccelerate.$(Symbol(string("vDSP_vrampmul2", suff)))(I0,I1,stride(I0,1),s,Ref(step),O0,O1,stride(O0,1),n)
            return (O0, O1)
        end
        function vrampmul2(I0::StridedVector{$T}, I1::StridedVector{$T}, start::$T, step::$T)
            O0 = similar(I0)
            O1 = similar(I1)
            vrampmul2!(O0, O1, I0, I1, start, step)
        end
    end
end

@doc "Stereo ramp multiply: multiply two vectors by the same ramp. Wraps [`vDSP_vrampmul2`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmul2)." vrampmul2

for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vavlin!(C::StridedVector{$T}, A::StridedVector{$T}, weight::$T)
            length(A) == length(C) ||
                throw(DimensionMismatch("C and A must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vavlin", suff)))(A,stride(A,1),Ref(weight),C,stride(C,1),length(A))
            return C
        end
        # Operand order matches the mutating `vavlin!(C, A, weight)`: C is the
        # running accumulator, A is the new sample vector.
        function vavlin(C::StridedVector{$T}, A::StridedVector{$T}, weight::$T)
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
        function vrsum!(result::StridedVector{$T}, X::StridedVector{$T}, scale::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vrsum", suff)))(X,stride(X,1),Ref(scale),result,stride(result,1),length(X))
            return result
        end
        function vrsum(X::StridedVector{$T}, scale::$T)
            result = similar(X)
            vrsum!(result, X, scale)
        end
        function vsimps!(result::StridedVector{$T}, X::StridedVector{$T}, step::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vsimps", suff)))(X,stride(X,1),Ref(step),result,stride(result,1),length(X))
            return result
        end
        function vsimps(X::StridedVector{$T}, step::$T)
            result = similar(X)
            vsimps!(result, X, step)
        end
        function vtrapz!(result::StridedVector{$T}, X::StridedVector{$T}, step::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vtrapz", suff)))(X,stride(X,1),Ref(step),result,stride(result,1),length(X))
            return result
        end
        function vtrapz(X::StridedVector{$T}, step::$T)
            result = similar(X)
            vtrapz!(result, X, step)
        end
        function vswsum!(result::StridedVector{$T}, X::StridedVector{$T}, window::Integer)
            (1 <= window <= length(X)) ||
                throw(ArgumentError("window ($window) must satisfy 1 <= window <= length(X) ($(length(X)))"))
            n_out = length(X) - window + 1
            length(result) >= n_out ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) - window + 1 ($n_out)"))
            _check_unit_stride(X, :vswsum!)
            _check_unit_stride(result, :vswsum!)
            LibAccelerate.$(Symbol(string("vDSP_vswsum", suff)))(X,1,result,1,n_out,window)
            return result
        end
        function vswsum(X::StridedVector{$T}, window::Integer)
            (1 <= window <= length(X)) ||
                throw(ArgumentError("window ($window) must satisfy 1 <= window <= length(X) ($(length(X)))"))
            n_out = length(X) - window + 1
            result = Vector{$T}(undef, n_out)
            vswsum!(result, X, window)
        end
        function vswmax!(result::StridedVector{$T}, X::StridedVector{$T}, window::Integer)
            (1 <= window <= length(X)) ||
                throw(ArgumentError("window ($window) must satisfy 1 <= window <= length(X) ($(length(X)))"))
            n_out = length(X) - window + 1
            length(result) >= n_out ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) - window + 1 ($n_out)"))
            _check_unit_stride(X, :vswmax!)
            _check_unit_stride(result, :vswmax!)
            LibAccelerate.$(Symbol(string("vDSP_vswmax", suff)))(X,1,result,1,n_out,window)
            return result
        end
        function vswmax(X::StridedVector{$T}, window::Integer)
            (1 <= window <= length(X)) ||
                throw(ArgumentError("window ($window) must satisfy 1 <= window <= length(X) ($(length(X)))"))
            n_out = length(X) - window + 1
            result = Vector{$T}(undef, n_out)
            vswmax!(result, X, window)
        end
    end
end

@doc "Running sum scaled by `scale`. Wraps [`vDSP_vrsum`](https://developer.apple.com/documentation/accelerate/vdsp_vrsum)." vrsum
@doc "Simpson's rule integration with step size `step`. Wraps [`vDSP_vsimps`](https://developer.apple.com/documentation/accelerate/vdsp_vsimps)." vsimps
@doc "Trapezoidal integration with step size `step`. Wraps [`vDSP_vtrapz`](https://developer.apple.com/documentation/accelerate/vdsp_vtrapz)." vtrapz
@doc "Sliding window sum with window size `window`. Returns a vector of length `length(X) - window + 1`. `window` must satisfy `1 ≤ window ≤ length(X)`; for `vswsum!`, `result` must have length `≥ length(X) - window + 1`. Wraps [`vDSP_vswsum`](https://developer.apple.com/documentation/accelerate/vdsp_vswsum)." vswsum
@doc "Sliding window maximum with window size `window`. Returns a vector of length `length(X) - window + 1`. `window` must satisfy `1 ≤ window ≤ length(X)`; for `vswmax!`, `result` must have length `≥ length(X) - window + 1`. Wraps [`vDSP_vswmax`](https://developer.apple.com/documentation/accelerate/vdsp_vswmax)." vswmax

# ============================================================
# vDSP Interpolation
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vintb!(result::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, t::$T)
            (length(A) == length(B) == length(result)) ||
                throw(DimensionMismatch("result, A and B must have the same length"))
            LibAccelerate.$(Symbol(string("vDSP_vintb", suff)))(A,stride(A,1),B,stride(B,1),Ref(t),result,stride(result,1),length(A))
            return result
        end
        function vintb(A::StridedVector{$T}, B::StridedVector{$T}, t::$T)
            result = similar(A)
            vintb!(result, A, B, t)
        end
        function vlint!(result::StridedVector{$T}, table::StridedVector{$T}, indices::StridedVector{$T})
            length(result) >= length(indices) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(indices) ($(length(indices)))"))
            _check_unit_stride(table, :vlint!)
            LibAccelerate.$(Symbol(string("vDSP_vlint", suff)))(table,indices,stride(indices,1),result,stride(result,1),length(indices),length(table))
            return result
        end
        function vlint(table::StridedVector{$T}, indices::StridedVector{$T})
            result = Vector{$T}(undef, length(indices))
            vlint!(result, table, indices)
        end
        function vqint!(result::StridedVector{$T}, table::StridedVector{$T}, indices::StridedVector{$T})
            length(result) >= length(indices) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(indices) ($(length(indices)))"))
            _check_unit_stride(table, :vqint!)
            LibAccelerate.$(Symbol(string("vDSP_vqint", suff)))(table,indices,stride(indices,1),result,stride(result,1),length(indices),length(table))
            return result
        end
        function vqint(table::StridedVector{$T}, indices::StridedVector{$T})
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
        function vpoly!(result::StridedVector{$T}, coeffs::StridedVector{$T}, X::StridedVector{$T})
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vpoly", suff)))(coeffs,stride(coeffs,1),X,stride(X,1),result,stride(result,1),length(X),length(coeffs) - 1)
            return result
        end
        function vpoly(coeffs::StridedVector{$T}, X::StridedVector{$T})
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
        function vnormalize!(result::StridedVector{$T}, X::StridedVector{$T})
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            mean_out = Ref{$T}(0.0)
            stddev_out = Ref{$T}(0.0)
            LibAccelerate.$(Symbol(string("vDSP_normalize", suff)))(X,stride(X,1),result,stride(result,1),mean_out,stddev_out,length(X))
            return (result, mean_out[], stddev_out[])
        end
        function vnormalize(X::StridedVector{$T})
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
        function nzcros(X::StridedVector{$T}, max_crossings::Integer=0)
            if max_crossings <= 0
                max_crossings = length(X)
            end
            _check_unit_stride(X, :nzcros)
            indices = Vector{UInt64}(undef, max_crossings)
            count = Ref{UInt64}(0)
            LibAccelerate.$(Symbol(string("vDSP_nzcros", suff)))(X,1,max_crossings,indices,count,length(X))
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
        function vdbcon!(result::StridedVector{$T}, X::StridedVector{$T}, ref::$T, power::Bool=true)
            flag = power ? UInt32(0) : UInt32(1)
            LibAccelerate.$(Symbol(string("vDSP_vdbcon", suff)))(X,stride(X,1),Ref(ref),result,stride(result,1),length(X),flag)
            return result
        end
        function vdbcon(X::StridedVector{$T}, ref::$T, power::Bool=true)
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
        function vclr!(C::StridedVector{$T})
            LibAccelerate.$(Symbol(string("vDSP_vclr", suff)))(C,stride(C,1),length(C))
            return C
        end
        function vfill!(C::StridedVector{$T}, a::$T)
            LibAccelerate.$(Symbol(string("vDSP_vfill", suff)))(Ref(a),C,stride(C,1),length(C))
            return C
        end
        function vswap!(A::StridedVector{$T}, B::StridedVector{$T})
            length(B) >= length(A) ||
                throw(DimensionMismatch("B length ($(length(B))) must be at least length(A) ($(length(A)))"))
            LibAccelerate.$(Symbol(string("vDSP_vswap", suff)))(A,stride(A,1),B,stride(B,1),length(A))
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
        function vgathr!(C::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{UInt})
            length(C) >= length(B) ||
                throw(DimensionMismatch("C length ($(length(C))) must be at least length(B) ($(length(B)))"))
            n = length(A)
            @inbounds for i in eachindex(B)
                (1 <= B[i] <= n) ||
                    throw(BoundsError(A, Int(B[i])))
            end
            _check_unit_stride(A, :vgathr!)
            LibAccelerate.$(Symbol(string("vDSP_vgathr", suff)))(A,B,stride(B,1),C,stride(C,1),length(B))
            return C
        end
        function vgathr(A::StridedVector{$T}, B::StridedVector{UInt})
            C = Vector{$T}(undef, length(B))
            vgathr!(C, A, B)
        end
        function vindex!(C::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T})
            length(C) >= length(B) ||
                throw(DimensionMismatch("C length ($(length(C))) must be at least length(B) ($(length(B)))"))
            _check_unit_stride(A, :vindex!)
            LibAccelerate.$(Symbol(string("vDSP_vindex", suff)))(A,B,stride(B,1),C,stride(C,1),length(B))
            return C
        end
        function vindex(A::StridedVector{$T}, B::StridedVector{$T})
            C = Vector{$T}(undef, length(B))
            vindex!(C, A, B)
        end
    end
end

@doc "Gather elements by index: `C[i] = A[B[i]]` where B contains 1-based UInt indices in `1:length(A)` (a `BoundsError` is thrown otherwise). Wraps [`vDSP_vgathr`](https://developer.apple.com/documentation/accelerate/vdsp_vgathr)." vgathr
@doc "Index with float indices: `C[i] = A[trunc(B[i])]` where B contains 0-based float indices. Wraps [`vDSP_vindex`](https://developer.apple.com/documentation/accelerate/vdsp_vindex)." vindex

# --- Generation ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vgen!(C::StridedVector{$T}, a::$T, b::$T)
            LibAccelerate.$(Symbol(string("vDSP_vgen", suff)))(Ref(a),Ref(b),C,stride(C,1),length(C))
            return C
        end
        function vgen(a::$T, b::$T, n::Integer)
            C = Vector{$T}(undef, n)
            vgen!(C, a, b)
        end
        function vgenp!(C::StridedVector{$T}, A::StridedVector{$T}, B::StridedVector{$T}, n::Integer)
            length(C) >= n ||
                throw(DimensionMismatch("C length ($(length(C))) must be at least n ($n)"))
            length(B) >= length(A) ||
                throw(DimensionMismatch("B length ($(length(B))) must be at least length(A) ($(length(A)))"))
            LibAccelerate.$(Symbol(string("vDSP_vgenp", suff)))(A,stride(A,1),B,stride(B,1),C,stride(C,1),n,length(A))
            return C
        end
        function vgenp(A::StridedVector{$T}, B::StridedVector{$T}, n::Integer)
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
        function vclipc!(result::StridedVector{$T}, X::StridedVector{$T}, low::$T, high::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            nlow = Ref{UInt64}(0)
            nhigh = Ref{UInt64}(0)
            LibAccelerate.$(Symbol(string("vDSP_vclipc", suff)))(X,stride(X,1),Ref(low),Ref(high),result,stride(result,1),length(X),nlow,nhigh)
            return (result, Int(nlow[]), Int(nhigh[]))
        end
        function vclipc(X::StridedVector{$T}, low::$T, high::$T)
            result = similar(X)
            vclipc!(result, X, low, high)
        end
        function vlim!(result::StridedVector{$T}, A::StridedVector{$T}, b::$T, c::$T)
            length(result) >= length(A) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(A) ($(length(A)))"))
            LibAccelerate.$(Symbol(string("vDSP_vlim", suff)))(A,stride(A,1),Ref(b),Ref(c),result,stride(result,1),length(A))
            return result
        end
        function vlim(A::StridedVector{$T}, b::$T, c::$T)
            result = similar(A)
            vlim!(result, A, b, c)
        end
        function vthrsc!(result::StridedVector{$T}, A::StridedVector{$T}, b::$T, c::$T)
            length(result) >= length(A) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(A) ($(length(A)))"))
            LibAccelerate.$(Symbol(string("vDSP_vthrsc", suff)))(A,stride(A,1),Ref(b),Ref(c),result,stride(result,1),length(A))
            return result
        end
        function vthrsc(A::StridedVector{$T}, b::$T, c::$T)
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
        function vsort!(X::StridedVector{$T}, ascending::Bool=true)
            _check_unit_stride(X, :vsort!)
            order = ascending ? Cint(1) : Cint(-1)
            LibAccelerate.$(Symbol(string("vDSP_vsort", suff)))(X,length(X),order)
            return X
        end
        function vsorti!(indices::StridedVector{UInt}, X::StridedVector{$T}, ascending::Bool=true)
            _check_unit_stride(X, :vsorti!)
            _check_unit_stride(indices, :vsorti!)
            order = ascending ? Cint(1) : Cint(-1)
            LibAccelerate.$(Symbol(string("vDSP_vsorti", suff)))(X,indices,C_NULL,length(X),order)
            return indices
        end
        function vsorti(X::StridedVector{$T}, ascending::Bool=true)
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
    vsorti!(indices::StridedVector{UInt}, X::StridedVector{T}, ascending::Bool=true)

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
        function vtabi!(D::StridedVector{$T}, A::StridedVector{$T}, s1::$T, s2::$T, C::StridedVector{$T})
            length(D) >= length(A) ||
                throw(DimensionMismatch("D length ($(length(D))) must be at least length(A) ($(length(A)))"))
            _check_unit_stride(C, :vtabi!)
            LibAccelerate.$(Symbol(string("vDSP_vtabi", suff)))(A,stride(A,1),Ref(s1),Ref(s2),C,length(C),D,stride(D,1),length(A))
            return D
        end
        function vtabi(A::StridedVector{$T}, s1::$T, s2::$T, C::StridedVector{$T})
            D = Vector{$T}(undef, length(A))
            vtabi!(D, A, s1, s2, C)
        end
    end
end

@doc "Table lookup with interpolation: `D[i] = C[clamp(s1*A[i]+s2, 0, M-1)]`. Wraps [`vDSP_vtabi`](https://developer.apple.com/documentation/accelerate/vdsp_vtabi)." vtabi

# --- Integer operations (Int32) ---
function vaddi!(C::StridedVector{Int32}, A::StridedVector{Int32}, B::StridedVector{Int32})
    (length(A) == length(B) == length(C)) ||
        throw(DimensionMismatch("C, A and B must have the same length"))
    LibAccelerate.vDSP_vaddi(A,stride(A,1),B,stride(B,1),C,stride(C,1),length(A))
    return C
end
function vaddi(A::StridedVector{Int32}, B::StridedVector{Int32})
    C = similar(A)
    vaddi!(C, A, B)
end

function vabsi!(C::StridedVector{Int32}, A::StridedVector{Int32})
    length(C) >= length(A) ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
    LibAccelerate.vDSP_vabsi(A,stride(A,1),C,stride(C,1),length(A))
    return C
end
function vabsi(A::StridedVector{Int32})
    C = similar(A)
    vabsi!(C, A)
end

function vfilli!(C::StridedVector{Int32}, a::Int32)
    LibAccelerate.vDSP_vfilli(Ref(a),C,stride(C,1),length(C))
    return C
end

function veqvi!(C::StridedVector{Int32}, A::StridedVector{Int32}, B::StridedVector{Int32})
    (length(A) == length(B) == length(C)) ||
        throw(DimensionMismatch("C, A and B must have the same length"))
    LibAccelerate.vDSP_veqvi(A,stride(A,1),B,stride(B,1),C,stride(C,1),length(A))
    return C
end
function veqvi(A::StridedVector{Int32}, B::StridedVector{Int32})
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
            LibAccelerate.$(Symbol(string("vDSP_mmul", suff)))(B,1,A,1,C,1,UInt64(n),UInt64(m),UInt64(p))
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
            LibAccelerate.$(Symbol(string("vDSP_mtrans", suff)))(A,1,C,1,UInt64(m),UInt64(n))
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
            LibAccelerate.$(Symbol(string("vDSP_mmov", suff)))(A,C,UInt64(m),UInt64(n),ta,tc)
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
            function ($fname!)(C::StridedVector{$intT}, A::StridedVector{$T})
                length(C) >= length(A) ||
                    throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
                LibAccelerate.$(Symbol(vdsp_name))(A,stride(A,1),C,stride(C,1),length(A))
                return C
            end
            function ($fname)(A::StridedVector{$T})
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
            function ($fname!)(C::StridedVector{$intT}, A::StridedVector{$T})
                length(C) >= length(A) ||
                    throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
                LibAccelerate.$(Symbol(vdsp_name))(A,stride(A,1),C,stride(C,1),length(A))
                return C
            end
            function ($fname)(A::StridedVector{$T})
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
            function ($fname!)(C::StridedVector{$intT}, A::StridedVector{$T})
                length(C) >= length(A) ||
                    throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
                LibAccelerate.$(Symbol(vdsp_name))(A,stride(A,1),C,stride(C,1),length(A))
                return C
            end
            function ($fname)(A::StridedVector{$T})
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
            function ($fname!)(C::StridedVector{$intT}, A::StridedVector{$T})
                length(C) >= length(A) ||
                    throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
                LibAccelerate.$(Symbol(vdsp_name))(A,stride(A,1),C,stride(C,1),length(A))
                return C
            end
            function ($fname)(A::StridedVector{$T})
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
            function ($fname!)(C::StridedVector{$T}, A::StridedVector{$intT})
                length(C) >= length(A) ||
                    throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
                LibAccelerate.$(Symbol(vdsp_name))(A,stride(A,1),C,stride(C,1),length(A))
                return C
            end
            function ($fname)(A::StridedVector{$intT}, ::Type{$T})
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
            function ($fname!)(C::StridedVector{$T}, A::StridedVector{$intT})
                length(C) >= length(A) ||
                    throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
                LibAccelerate.$(Symbol(vdsp_name))(A,stride(A,1),C,stride(C,1),length(A))
                return C
            end
            function ($fname)(A::StridedVector{$intT}, ::Type{$T})
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
            LibAccelerate.$(Symbol(string("vDSP_f3x3", suff)))(A,UInt64(nc),UInt64(nr),F,C)
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
            LibAccelerate.$(Symbol(string("vDSP_f5x5", suff)))(A,UInt64(nc),UInt64(nr),F,C)
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
            LibAccelerate.$(Symbol(string("vDSP_imgfir", suff)))(A,UInt64(nc),UInt64(nr),F,C,UInt64(fc),UInt64(fr))
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

# ============================================================
# Batch 9: dotpr2, vrampmuladd, fixed-point (Q1.15 / Q8.24) kernels,
# 24-bit packed conversions, and misc integer vector ops
# See: https://developer.apple.com/documentation/accelerate/vdsp
# ============================================================

# --- vDSP_dotpr2: one X vector dotted against two Y vectors ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function dotpr2(A0::StridedVector{$T}, A1::StridedVector{$T}, B::StridedVector{$T})
            length(A0) == length(A1) == length(B) ||
                throw(DimensionMismatch("A0, A1, and B must all have the same length"))
            c0 = Ref{$T}(0)
            c1 = Ref{$T}(0)
            LibAccelerate.$(Symbol(string("vDSP_dotpr2", suff)))(A0,stride(A0,1),A1,stride(A1,1),B,stride(B,1),c0,c1,length(B))
            return (c0[], c1[])
        end
    end
end

@doc "Dual dot product: `(sum(A0 .* B), sum(A1 .* B))` — one vector `B` dotted against both `A0` and `A1`. Wraps [`vDSP_dotpr2`](https://developer.apple.com/documentation/accelerate/vdsp_dotpr2)." dotpr2

# --- vDSP_vrampmuladd / vrampmuladd2: ramp-multiply then accumulate ---
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vrampmuladd!(result::StridedVector{$T}, X::StridedVector{$T}, start::$T, step::$T)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            s = Ref{$T}(start)
            LibAccelerate.$(Symbol(string("vDSP_vrampmuladd", suff)))(X,stride(X,1),s,Ref(step),result,stride(result,1),length(X))
            return result
        end
        function vrampmuladd(X::StridedVector{$T}, start::$T, step::$T, C::StridedVector{$T})
            length(C) == length(X) ||
                throw(DimensionMismatch("C and X must have the same length"))
            result = copy(C)
            vrampmuladd!(result, X, start, step)
        end
        function vrampmuladd2!(O0::StridedVector{$T}, O1::StridedVector{$T}, I0::StridedVector{$T}, I1::StridedVector{$T}, start::$T, step::$T)
            n = length(I0)
            length(I1) == n ||
                throw(DimensionMismatch("I0 and I1 must have the same length"))
            length(O0) >= n && length(O1) >= n ||
                throw(DimensionMismatch("O0 and O1 must be at least length(I0)/length(I1)"))
            stride(I0,1) == stride(I1,1) ||
                throw(ArgumentError("vrampmuladd2!: I0 and I1 must have the same stride (vDSP shares one input stride)"))
            stride(O0,1) == stride(O1,1) ||
                throw(ArgumentError("vrampmuladd2!: O0 and O1 must have the same stride (vDSP shares one output stride)"))
            s = Ref{$T}(start)
            LibAccelerate.$(Symbol(string("vDSP_vrampmuladd2", suff)))(I0,I1,stride(I0,1),s,Ref(step),O0,O1,stride(O0,1),n)
            return (O0, O1)
        end
        function vrampmuladd2(I0::StridedVector{$T}, I1::StridedVector{$T}, start::$T, step::$T, C0::StridedVector{$T}, C1::StridedVector{$T})
            length(C0) == length(I0) && length(C1) == length(I1) ||
                throw(DimensionMismatch("C0/C1 must match I0/I1 in length"))
            O0 = copy(C0)
            O1 = copy(C1)
            vrampmuladd2!(O0, O1, I0, I1, start, step)
        end
    end
end

@doc """
Ramp-multiply then accumulate: `result[i] = C[i] + X[i] * (start + i*step)` for `i = 0, ..., length(X)-1`.
The mutating `vrampmuladd!(result, X, start, step)` reads/writes `result` in place as the accumulator
(i.e. `result[i] += X[i] * ramp[i]`). Wraps [`vDSP_vrampmuladd`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmuladd).
""" vrampmuladd
@doc """
Stereo ramp-multiply then accumulate: multiplies `I0`/`I1` by the same generated ramp and accumulates
into `C0`/`C1`. Wraps [`vDSP_vrampmuladd2`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmuladd2).
""" vrampmuladd2

# --- Fixed-point Q-format kernels (Q1.15 on Int16, Q8.24 on Int32) ---
#
# These operate directly on the raw integer storage: a Q1.15 value `v::Int16`
# represents the real number `v / 32768`, and a Q8.24 value `v::Int32`
# represents `v / 2^24`. Empirically verified (Jul 2026): the intrinsics
# compute exactly the same result as the floating-point analog applied to the
# decoded fractional values, then re-encoded at the *same* scale (no extra
# shift) — e.g. `dotpr_s1_15(A, B) == round(Int16, 32768 * dot(A/32768, B/32768))`
# for in-range results (out-of-range accumulation saturates/wraps per vDSP's
# fixed-point semantics; keep inputs small enough to avoid that in tests).
# One function per Q-format suffix (rather than a runtime `Val` flag) so the
# element type alone fully determines the format.
for (intT, qsuff) in ((Int16, "s1_15"), (Int32, "s8_24"))
    dotprname = Symbol("dotpr_", qsuff)
    dotpr2name = Symbol("dotpr2_", qsuff)
    vrampmulname! = Symbol("vrampmul_", qsuff, "!")
    vrampmulname = Symbol("vrampmul_", qsuff)
    vrampmul2name! = Symbol("vrampmul2_", qsuff, "!")
    vrampmul2name = Symbol("vrampmul2_", qsuff)
    vrampmuladdname! = Symbol("vrampmuladd_", qsuff, "!")
    vrampmuladdname = Symbol("vrampmuladd_", qsuff)
    vrampmuladd2name! = Symbol("vrampmuladd2_", qsuff, "!")
    vrampmuladd2name = Symbol("vrampmuladd2_", qsuff)
    @eval begin
        function ($dotprname)(A::StridedVector{$intT}, B::StridedVector{$intT})
            length(A) == length(B) ||
                throw(DimensionMismatch("A and B must have the same length"))
            c = Ref{$intT}(0)
            LibAccelerate.$(Symbol(string("vDSP_dotpr_", qsuff)))(A,stride(A,1),B,stride(B,1),c,length(A))
            return c[]
        end
        function ($dotpr2name)(A0::StridedVector{$intT}, A1::StridedVector{$intT}, B::StridedVector{$intT})
            length(A0) == length(A1) == length(B) ||
                throw(DimensionMismatch("A0, A1, and B must all have the same length"))
            c0 = Ref{$intT}(0)
            c1 = Ref{$intT}(0)
            LibAccelerate.$(Symbol(string("vDSP_dotpr2_", qsuff)))(A0,stride(A0,1),A1,stride(A1,1),B,stride(B,1),c0,c1,length(B))
            return (c0[], c1[])
        end
        function ($vrampmulname!)(result::StridedVector{$intT}, X::StridedVector{$intT}, start::$intT, step::$intT)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vrampmul_", qsuff)))(X,stride(X,1),Ref(start),Ref(step),result,stride(result,1),length(X))
            return result
        end
        function ($vrampmulname)(X::StridedVector{$intT}, start::$intT, step::$intT)
            result = similar(X)
            ($vrampmulname!)(result, X, start, step)
        end
        function ($vrampmul2name!)(O0::StridedVector{$intT}, O1::StridedVector{$intT}, I0::StridedVector{$intT}, I1::StridedVector{$intT}, start::$intT, step::$intT)
            n = length(I0)
            length(I1) == n ||
                throw(DimensionMismatch("I0 and I1 must have the same length"))
            length(O0) >= n && length(O1) >= n ||
                throw(DimensionMismatch("O0 and O1 must be at least length(I0)/length(I1)"))
            stride(I0,1) == stride(I1,1) ||
                throw(ArgumentError("$($(QuoteNode(vrampmul2name!))): I0 and I1 must have the same stride (vDSP shares one input stride)"))
            stride(O0,1) == stride(O1,1) ||
                throw(ArgumentError("$($(QuoteNode(vrampmul2name!))): O0 and O1 must have the same stride (vDSP shares one output stride)"))
            LibAccelerate.$(Symbol(string("vDSP_vrampmul2_", qsuff)))(I0,I1,stride(I0,1),Ref(start),Ref(step),O0,O1,stride(O0,1),n)
            return (O0, O1)
        end
        function ($vrampmul2name)(I0::StridedVector{$intT}, I1::StridedVector{$intT}, start::$intT, step::$intT)
            O0 = similar(I0)
            O1 = similar(I1)
            ($vrampmul2name!)(O0, O1, I0, I1, start, step)
        end
        function ($vrampmuladdname!)(result::StridedVector{$intT}, X::StridedVector{$intT}, start::$intT, step::$intT)
            length(result) >= length(X) ||
                throw(DimensionMismatch("result length ($(length(result))) must be at least length(X) ($(length(X)))"))
            LibAccelerate.$(Symbol(string("vDSP_vrampmuladd_", qsuff)))(X,stride(X,1),Ref(start),Ref(step),result,stride(result,1),length(X))
            return result
        end
        function ($vrampmuladdname)(X::StridedVector{$intT}, start::$intT, step::$intT, C::StridedVector{$intT})
            length(C) == length(X) ||
                throw(DimensionMismatch("C and X must have the same length"))
            result = copy(C)
            ($vrampmuladdname!)(result, X, start, step)
        end
        function ($vrampmuladd2name!)(O0::StridedVector{$intT}, O1::StridedVector{$intT}, I0::StridedVector{$intT}, I1::StridedVector{$intT}, start::$intT, step::$intT)
            n = length(I0)
            length(I1) == n ||
                throw(DimensionMismatch("I0 and I1 must have the same length"))
            length(O0) >= n && length(O1) >= n ||
                throw(DimensionMismatch("O0 and O1 must be at least length(I0)/length(I1)"))
            stride(I0,1) == stride(I1,1) ||
                throw(ArgumentError("$($(QuoteNode(vrampmuladd2name!))): I0 and I1 must have the same stride (vDSP shares one input stride)"))
            stride(O0,1) == stride(O1,1) ||
                throw(ArgumentError("$($(QuoteNode(vrampmuladd2name!))): O0 and O1 must have the same stride (vDSP shares one output stride)"))
            LibAccelerate.$(Symbol(string("vDSP_vrampmuladd2_", qsuff)))(I0,I1,stride(I0,1),Ref(start),Ref(step),O0,O1,stride(O0,1),n)
            return (O0, O1)
        end
        function ($vrampmuladd2name)(I0::StridedVector{$intT}, I1::StridedVector{$intT}, start::$intT, step::$intT, C0::StridedVector{$intT}, C1::StridedVector{$intT})
            length(C0) == length(I0) && length(C1) == length(I1) ||
                throw(DimensionMismatch("C0/C1 must match I0/I1 in length"))
            O0 = copy(C0)
            O1 = copy(C1)
            ($vrampmuladd2name!)(O0, O1, I0, I1, start, step)
        end
    end
    @eval @doc "Fixed-point ($($qsuff)) dot product: `sum(A .* B)` computed and stored in Q-format. Wraps [`vDSP_dotpr_$($qsuff)`](https://developer.apple.com/documentation/accelerate/vdsp_dotpr_$($qsuff))." $dotprname
    @eval @doc "Fixed-point ($($qsuff)) dual dot product: `B` dotted against both `A0` and `A1`. Wraps [`vDSP_dotpr2_$($qsuff)`](https://developer.apple.com/documentation/accelerate/vdsp_dotpr2_$($qsuff))." $dotpr2name
    @eval @doc "Fixed-point ($($qsuff)) ramp multiply. Wraps [`vDSP_vrampmul_$($qsuff)`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmul_$($qsuff))." $vrampmulname
    @eval @doc "Fixed-point ($($qsuff)) stereo ramp multiply. Wraps [`vDSP_vrampmul2_$($qsuff)`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmul2_$($qsuff))." $vrampmul2name
    @eval @doc "Fixed-point ($($qsuff)) ramp-multiply then accumulate. Wraps [`vDSP_vrampmuladd_$($qsuff)`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmuladd_$($qsuff))." $vrampmuladdname
    @eval @doc "Fixed-point ($($qsuff)) stereo ramp-multiply then accumulate. Wraps [`vDSP_vrampmuladd2_$($qsuff)`](https://developer.apple.com/documentation/accelerate/vdsp_vrampmuladd2_$($qsuff))." $vrampmuladd2name
end

# --- 24-bit packed integer conversions ---
#
# vDSP represents packed 24-bit integers as a 3-byte struct
# (`LibAccelerate.vDSP_int24`/`vDSP_uint24`, no padding). We surface these to
# users as ordinary `Int32`/`UInt32` vectors holding values in the 24-bit
# range, packing/unpacking to the 3-byte layout around each ccall.
@inline function _pack_int24(x::Integer)
    (-8388608 <= x <= 8388607) ||
        throw(ArgumentError("value $x out of range for a 24-bit signed integer (-8388608:8388607)"))
    u = reinterpret(UInt32, Int32(x)) & 0x00ffffff
    return LibAccelerate.vDSP_int24((UInt8(u & 0xff), UInt8((u >> 8) & 0xff), UInt8((u >> 16) & 0xff)))
end
@inline function _unpack_int24(v)
    b0, b1, b2 = v.bytes
    u = UInt32(b0) | (UInt32(b1) << 8) | (UInt32(b2) << 16)
    (u & 0x800000) != 0 && (u |= 0xff000000)
    return reinterpret(Int32, u)
end
@inline function _pack_uint24(x::Integer)
    (0 <= x <= 16777215) ||
        throw(ArgumentError("value $x out of range for a 24-bit unsigned integer (0:16777215)"))
    u = UInt32(x)
    return LibAccelerate.vDSP_uint24((UInt8(u & 0xff), UInt8((u >> 8) & 0xff), UInt8((u >> 16) & 0xff)))
end
@inline function _unpack_uint24(v)
    b0, b1, b2 = v.bytes
    return UInt32(b0) | (UInt32(b1) << 8) | (UInt32(b2) << 16)
end

function vflt24!(C::StridedVector{Float32}, A::AbstractVector{<:Integer})
    n = length(A)
    length(C) >= n ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($n)"))
    packed = Vector{LibAccelerate.vDSP_int24}(undef, n)
    @inbounds for i in 1:n
        packed[i] = _pack_int24(A[i])
    end
    LibAccelerate.vDSP_vflt24(packed, 1, C, stride(C,1), n)
    return C
end
vflt24(A::AbstractVector{<:Integer}) = vflt24!(Vector{Float32}(undef, length(A)), A)

function vfltu24!(C::StridedVector{Float32}, A::AbstractVector{<:Integer})
    n = length(A)
    length(C) >= n ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($n)"))
    packed = Vector{LibAccelerate.vDSP_uint24}(undef, n)
    @inbounds for i in 1:n
        packed[i] = _pack_uint24(A[i])
    end
    LibAccelerate.vDSP_vfltu24(packed, 1, C, stride(C,1), n)
    return C
end
vfltu24(A::AbstractVector{<:Integer}) = vfltu24!(Vector{Float32}(undef, length(A)), A)

function vfltsm24!(C::StridedVector{Float32}, A::AbstractVector{<:Integer}, b::Float32)
    n = length(A)
    length(C) >= n ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($n)"))
    packed = Vector{LibAccelerate.vDSP_int24}(undef, n)
    @inbounds for i in 1:n
        packed[i] = _pack_int24(A[i])
    end
    LibAccelerate.vDSP_vfltsm24(packed, 1, Ref(b), C, stride(C,1), n)
    return C
end
vfltsm24(A::AbstractVector{<:Integer}, b::Float32) = vfltsm24!(Vector{Float32}(undef, length(A)), A, b)

function vfltsmu24!(C::StridedVector{Float32}, A::AbstractVector{<:Integer}, b::Float32)
    n = length(A)
    length(C) >= n ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($n)"))
    packed = Vector{LibAccelerate.vDSP_uint24}(undef, n)
    @inbounds for i in 1:n
        packed[i] = _pack_uint24(A[i])
    end
    LibAccelerate.vDSP_vfltsmu24(packed, 1, Ref(b), C, stride(C,1), n)
    return C
end
vfltsmu24(A::AbstractVector{<:Integer}, b::Float32) = vfltsmu24!(Vector{Float32}(undef, length(A)), A, b)

function vsmfix24!(C::AbstractVector{<:Integer}, A::StridedVector{Float32}, b::Float32)
    n = length(A)
    length(C) >= n ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($n)"))
    packed = Vector{LibAccelerate.vDSP_int24}(undef, n)
    LibAccelerate.vDSP_vsmfix24(A, stride(A,1), Ref(b), packed, 1, n)
    @inbounds for i in 1:n
        C[i] = _unpack_int24(packed[i])
    end
    return C
end
vsmfix24(A::StridedVector{Float32}, b::Float32) = vsmfix24!(Vector{Int32}(undef, length(A)), A, b)

function vsmfixu24!(C::AbstractVector{<:Integer}, A::StridedVector{Float32}, b::Float32)
    n = length(A)
    length(C) >= n ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($n)"))
    packed = Vector{LibAccelerate.vDSP_uint24}(undef, n)
    LibAccelerate.vDSP_vsmfixu24(A, stride(A,1), Ref(b), packed, 1, n)
    @inbounds for i in 1:n
        C[i] = _unpack_uint24(packed[i])
    end
    return C
end
vsmfixu24(A::StridedVector{Float32}, b::Float32) = vsmfixu24!(Vector{UInt32}(undef, length(A)), A, b)

@doc "Convert packed 24-bit signed integers (given as `Int32` values in `-8388608:8388607`) to `Float32`. Wraps [`vDSP_vflt24`](https://developer.apple.com/documentation/accelerate/vdsp_vflt24)." vflt24
@doc "Convert packed 24-bit unsigned integers (given as `UInt32`/`Integer` values in `0:16777215`) to `Float32`. Wraps [`vDSP_vfltu24`](https://developer.apple.com/documentation/accelerate/vdsp_vfltu24)." vfltu24
@doc "Convert packed 24-bit signed integers to `Float32`, then scale by `b`: `C[i] = Float32(A[i]) * b`. Wraps [`vDSP_vfltsm24`](https://developer.apple.com/documentation/accelerate/vdsp_vfltsm24)." vfltsm24
@doc "Convert packed 24-bit unsigned integers to `Float32`, then scale by `b`: `C[i] = Float32(A[i]) * b`. Wraps [`vDSP_vfltsmu24`](https://developer.apple.com/documentation/accelerate/vdsp_vfltsmu24)." vfltsmu24
@doc "Scale `A` by `b` and truncate to packed 24-bit signed integers (returned as `Int32`): `C[i] = trunc(Int32, A[i] * b)`. Wraps [`vDSP_vsmfix24`](https://developer.apple.com/documentation/accelerate/vdsp_vsmfix24)." vsmfix24
@doc "Scale `A` by `b` and truncate to packed 24-bit unsigned integers (returned as `UInt32`): `C[i] = trunc(UInt32, A[i] * b)`. Wraps [`vDSP_vsmfixu24`](https://developer.apple.com/documentation/accelerate/vdsp_vsmfixu24)." vsmfixu24

# --- Integer vector ops ---
function vdivi!(C::StridedVector{Cint}, A::StridedVector{Cint}, B::StridedVector{Cint})
    length(A) == length(B) ||
        throw(DimensionMismatch("A and B must have the same length"))
    length(C) >= length(A) ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
    LibAccelerate.vDSP_vdivi(B, stride(B,1), A, stride(A,1), C, stride(C,1), length(A))
    return C
end
vdivi(A::StridedVector{Cint}, B::StridedVector{Cint}) = vdivi!(similar(A), A, B)

function vsaddi!(C::StridedVector{Cint}, A::StridedVector{Cint}, b::Cint)
    length(C) >= length(A) ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
    LibAccelerate.vDSP_vsaddi(A, stride(A,1), Ref(b), C, stride(C,1), length(A))
    return C
end
vsaddi(A::StridedVector{Cint}, b::Cint) = vsaddi!(similar(A), A, b)

function vsdivi!(C::StridedVector{Cint}, A::StridedVector{Cint}, b::Cint)
    length(C) >= length(A) ||
        throw(DimensionMismatch("C length ($(length(C))) must be at least length(A) ($(length(A)))"))
    LibAccelerate.vDSP_vsdivi(A, stride(A,1), Ref(b), C, stride(C,1), length(A))
    return C
end
vsdivi(A::StridedVector{Cint}, b::Cint) = vsdivi!(similar(A), A, b)

@doc "Integer vector divide: `C[i] = div(A[i], B[i])` (truncating). Wraps [`vDSP_vdivi`](https://developer.apple.com/documentation/accelerate/vdsp_vdivi)." vdivi
@doc "Integer scalar add: `C[i] = A[i] + b`. Wraps [`vDSP_vsaddi`](https://developer.apple.com/documentation/accelerate/vdsp_vsaddi)." vsaddi
@doc "Integer scalar divide: `C[i] = div(A[i], b)` (truncating). Wraps [`vDSP_vsdivi`](https://developer.apple.com/documentation/accelerate/vdsp_vsdivi)." vsdivi

# --- vDSP_vgathra: gather via an array of pointers ---
#
# `A` is the caller-owned pointer-source array (its full extent, including any
# entries skipped by `ptrstride`); `C` determines how many outputs (`N`) are
# gathered. vDSP reads pointer entries `1, 1+ptrstride, 1+2*ptrstride, ...` (one
# per output), so `A` must contain at least `(length(C)-1)*ptrstride + 1` entries
# or the call reads past the end of the (GC-owned) pointer buffer.
for (T, suff) in ((Float32, ""), (Float64, "D"))
    @eval begin
        function vgathra!(C::StridedVector{$T}, A::AbstractVector{<:StridedVector{$T}}, ptrstride::Integer=1)
            ptrstride >= 1 ||
                throw(ArgumentError("vgathra!: pointer-array stride must be positive"))
            n = length(C)
            needed = n == 0 ? 0 : (n - 1) * ptrstride + 1
            length(A) >= needed ||
                throw(DimensionMismatch("A must contain at least $needed pointer entries for $n outputs at stride $ptrstride, got $(length(A))"))
            GC.@preserve A begin
                ptrs = Vector{Ptr{$T}}(undef, length(A))
                @inbounds for i in 1:length(A)
                    _check_unit_stride(A[i], :vgathra!)
                    isempty(A[i]) && throw(BoundsError(A[i], 1))
                    ptrs[i] = pointer(A[i])
                end
                LibAccelerate.$(Symbol(string("vDSP_vgathra", suff)))(ptrs, ptrstride, C, stride(C,1), n)
            end
            return C
        end
        function vgathra(A::AbstractVector{<:StridedVector{$T}}, ptrstride::Integer=1)
            n = (ptrstride < 1 || isempty(A)) ? 0 : div(length(A) - 1, ptrstride) + 1
            C = Vector{$T}(undef, n)
            vgathra!(C, A, ptrstride)
        end
    end
end

@doc """
    vgathra(A::AbstractVector{<:StridedVector{T}}, ptrstride::Integer=1) -> Vector{T}

Gather the first element of each source vector in `A` into a contiguous result:
`C[i] = A[1 + (i-1)*ptrstride][1]`. `A` is a caller-owned array of unit-stride
vectors; internally an array of their pointers is built and passed to vDSP
(`GC.@preserve`d for the duration of the call). `ptrstride` selects every
`ptrstride`-th pointer from that array (must be positive); the allocating form
derives the output length from `length(A)` and `ptrstride`, while `vgathra!`
derives the requested output count `N` from `length(C)` and requires `A` to hold
at least `(N-1)*ptrstride + 1` entries. Wraps [`vDSP_vgathra`](https://developer.apple.com/documentation/accelerate/vdsp_vgathra)
/ [`vDSP_vgathraD`](https://developer.apple.com/documentation/accelerate/vdsp_vgathrad).
""" vgathra
