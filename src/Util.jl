## Util.jl##

tupletypelength(a)=length(a.parameters)


macro replaceBase(fs...)
    b = Expr(:block)
    for f in fs
        if f == :./
            fa = :fdiv
        elseif f == :.^
            fa = :pow
        elseif f == :.+
            fa = :vadd
        elseif f == :.-
            fa = :vsub
        elseif f == :.*
            fa = :vmul
        else
            fa = f
        end
        if fa in (:ceil,:floor,:sqrt,:rsqrt,:rec,
                  :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
                  :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,:atan,
                  :sinh,:cosh,:tanh,:asinh,:acosh,:atanh,
                  :trunc,:round,:exponent,:abs,:sincos,:cis)
            e = quote
                (Base.$f)(X::Array{Float64}) = ($fa)(X)
                (Base.$f)(X::Array{Float32}) = ($fa)(X)
            end
        elseif fa in (:copysign,:atan2,:rem)
            e = quote
                (Base.$f)(X::Array{Float32}, Y::Array{Float32}) = ($fa)(X,Y)
                (Base.$f)(X::Array{Float64}, Y::Array{Float64}) = ($fa)(X,Y)
            end
        elseif fa in (:fdiv,:pow,:vadd,:vsub,:vmul)
            fb = Symbol(replace(string(f), r"^\.(.*)", s"\1"))
            e = quote
                Base.broadcast(::typeof(Base.$fb), X::Array{Float64}, Y::Array{Float64}) = ($fa)(X,Y)
                Base.broadcast(::typeof(Base.$fb), X::Array{Float32}, Y::Array{Float32}) = ($fa)(X,Y)
            end
        else
            error("Function $f not defined by AppleAccelerate.jl")
        end
        push!(b.args,e)
    end
    b
end
