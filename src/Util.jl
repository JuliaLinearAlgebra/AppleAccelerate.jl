## Util.jl##

tupletypelength(a)=length(a.parameters)


macro replaceBase(fs...)
    b = Expr(:block)
    for f in fs
        if f == :/
            fa = :fdiv
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
        fa! = Symbol("$(fa)!")
        if fa in (:ceil,:floor,:sqrt,:rsqrt,:rec,
                  :exp,:exp2,:expm1,:log,:log1p,:log2,:log10,
                  :sin,:sinpi,:cos,:cospi,:tan,:tanpi,:asin,:acos,
                  :sinh,:cosh,:tanh,:asinh,:acosh,:atanh,
                  :trunc,:round,:exponent,:abs,:sincos)

            e = quote
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copyto!)(dest::Array{Float64, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
                (Base.copyto!)(dest::Array{Float32, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
            end
        elseif fa in (:copysign,:pow,:rem,:fdiv, :vadd, :vsub, :vmul)
            e = quote
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}, Array{Float64, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}, Array{Float32, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copyto!)(dest::Array{Float64, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}, Array{Float64, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
                (Base.copyto!)(dest::Array{Float32, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}, Array{Float32, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
            end
        elseif fa == :atan
            e = quote
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copyto!)(dest::Array{Float64, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
                (Base.copyto!)(dest::Array{Float32, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}, Array{Float64, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}, Array{Float32, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copyto!)(dest::Array{Float64, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}, Array{Float64, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
                (Base.copyto!)(dest::Array{Float32, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}, Array{Float32, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
            end
        elseif fa == :cis
            e = quote
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copy)(bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}}}) where {Style, Axes, N} = ($fa)(bc.args...)
                (Base.copyto!)(dest::Array{Complex{Float64}, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float64, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
                (Base.copyto!)(dest::Array{Complex{Float32}, N}, bc::Base.Broadcast.Broadcasted{Style, Axes, typeof(Base.$f), Tuple{Array{Float32, N}}}) where {Style, Axes, N} = ($fa!)(dest, bc.args...)
            end
        else
            error("Function $f not defined by AppleAccelerate.jl")
        end
        push!(b.args,e)
    end
    b
end
