## Util.jl##

tupletypelength(a)=length(a.parameters)

@inline maybecopy(x::T) where {T <: Base.Broadcast.Broadcasted} = copy(x)
@inline maybecopy(x::T) where {T <: Array} = x

const OPS = Dict{Symbol,Tuple{Symbol, Symbol, Symbol}}(:+ => (:vadd, :vsadd, :vsadd),
:- => (:vsub, :vssub, :svsub),
:* => (:vmul, :vsmul, :vsmul),
:/ => (:vdiv, :vsdiv, :vsdiv),)

macro replaceBase(fs...)
    b = Expr(:block)
    for f in fs
        if f == :/
            fa = :div_float
        elseif f == :^
            fa = :pow
        elseif f == :+
            fa = :vadd
        elseif f == :-
            fa = :vsub
        elseif f == :*
            fa = :vmul
        else
            fa = f
        end
        arg_consumed = false
        if fa in (:ceil,:floor,:sqrt,:rsqrt,:rec,
                  :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
                  :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
                  :sinh,:cosh,:tanh,:asinh,:acosh,:atanh,
                  :trunc,:round,:exponent,:abs,:sincos,:cis)
            e = quote
                (Base.$f)(X::Array{T}) where {T <: Union{Float64,Float32}} = ($fa)(X)
                (Base.$f)(X::Union{Float64,Float32}) = ($fa)([X])[1]
                Base.broadcasted(::typeof(Base.$f), arg::Union{Array{F,N},Base.Broadcast.Broadcasted}) where {N,F<:Union{Float32,Float64}} = ($fa)(maybecopy(arg))
            end
            arg_consumed = true
        end
        if fa in (:copysign,:atan,:pow,:rem)
            e = quote
                (Base.$f)(X::Array{T},Y::Array{T}) where {T <: Union{Float32,Float64}} = ($fa)(X,Y)
                (Base.$f)(X::T,Y::T) where {T <: Union{Float32,Float64}} = ($fa)([X],[Y])[1]
                Base.broadcasted(::typeof(Base.$f), arg1::Union{Array{F, N},Base.Broadcast.Broadcasted}, arg2::Union{Array{F, N},Base.Broadcast.Broadcasted}) where {N,F<:Union{Float32,Float64}} = ($fa)(maybecopy(arg1), maybecopy(arg2))
            end
            arg_consumed = true
        end
        if f in (:+,:-,:*,:/)
            e = quote
                (Base.$f)(X::Array{T},Y::Array{T}) where {T <: Union{Float32,Float64}} = ($(OPS[f][1]))(X,Y)
                (Base.$f)(X::T,Y::T) where {T <: Union{Float32,Float64}} = ($(OPS[f][1]))([X],[Y])[1]
                Base.broadcasted(::typeof(Base.$f), arg1::Union{Array{F, N},Base.Broadcast.Broadcasted}, arg2::Union{Array{F, N},Base.Broadcast.Broadcasted}) where {N,F<:Union{Float32,Float64}} = ($(OPS[f][1]))(maybecopy(arg1), maybecopy(arg2))

                Base.broadcasted(::typeof(Base.$f), arg1::Union{Array{T, N},Base.Broadcast.Broadcasted}, arg2::T) where {N, T <: Union{Float32,Float64}} = ($(OPS[f][2]))(maybecopy(arg1), arg2)
                Base.broadcasted(::typeof(Base.$f), arg1::T, arg2::Union{Array{T, N},Base.Broadcast.Broadcasted}) where {N, T <: Union{Float32,Float64}} = ($(OPS[f][3]))(maybecopy(arg2), arg1)
        end
            arg_consumed = true
        end
        if !arg_consumed
            error("Function $f not defined by AppleAccelerate.jl")
        end
        push!(b.args,e)
    end
    b
end
