## Util.jl##

tupletypelength(a)=length(a.parameters)


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
            end
            arg_consumed = true
        end
        if fa in (:copysign,:atan,:pow,:rem,:div_float, :vadd, :vsub, :vmul)
            e = quote
                (Base.$f)(X::Array{T},Y::Array{T}) where {T <: Union{Float32,Float64}} = ($fa)(X,Y)
                (Base.$f)(X::T,Y::T) where {T <: Union{Float32,Float64}} = ($fa)([X],[Y])[1]
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
