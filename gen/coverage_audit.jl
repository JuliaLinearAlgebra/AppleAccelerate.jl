# Exact raw-layer coverage audit.
#
#   julia --project=. gen/coverage_audit.jl [prefix]      # prefix defaults to "vDSP"
#
# Lists the `LibAccelerate` functions with the given name prefix that no method in
# the idiomatic layer can reach. Name-matching over the source text (grep, or the
# stem heuristic) badly over-reports here, because most wrappers build the raw name
# at macro-expansion time — `Symbol(string("vDSP_vfix", intname, suff))` — so the
# literal symbol never appears in `src/*.jl`. This script instead walks the *lowered
# IR* of every method defined in `AppleAccelerate` (where those names are already
# resolved) and collects every `LibAccelerate.foo` reference plus every direct
# `ccall` symbol. What survives is genuinely unreferenced.

using AppleAccelerate

const AA = AppleAccelerate
const L  = AA.LibAccelerate

const used = Set{Symbol}()

_ccall_name(s::QuoteNode) = _ccall_name(s.value)
_ccall_name(s::Expr) = s.head === :tuple ? _ccall_name(s.args[1]) : nothing
_ccall_name(s::Tuple) = _ccall_name(s[1])
_ccall_name(s::Union{Symbol,AbstractString}) = Symbol(s)
_ccall_name(_) = nothing

function walk(x)
    if x isa GlobalRef
        x.mod === L && push!(used, x.name)
    elseif x isa QuoteNode
        # `getproperty(LibAccelerate, :foo)` lowers to a QuoteNode'd symbol
        x.value isa Symbol && push!(used, x.value)
    elseif x isa Expr
        if x.head === :foreigncall
            n = _ccall_name(x.args[1])
            n === nothing || push!(used, n)
        end
        foreach(walk, x.args)
    end
end

_methods(f::Function) = methods(f)
# `methods(FFTSetup)` misses inner constructors of a parametric struct (they are
# methods of `Type{FFTSetup{T}}`, not of the `UnionAll`), so match on `Type{<:T}`.
_methods(T::Type) = [m.method for m in
    Base._methods_by_ftype(Tuple{Type{<:T},Vararg{Any}}, -1, Base.get_world_counter())]

const seen = Set{Module}()
function scan(m::Module)
    (m in seen || m === L) && return
    push!(seen, m)
    for n in names(m; all = true)
        isdefined(m, n) || continue
        v = getfield(m, n)
        if v isa Module
            parentmodule(v) === m && scan(v)
        elseif v isa Function || v isa Type
            for meth in _methods(v)
                meth.module in seen || continue
                ci = try Base.uncompressed_ir(meth) catch; nothing end
                ci === nothing || foreach(walk, ci.code)
            end
        end
    end
end

function main(prefix)
    scan(AA)
    raw = sort!([n for n in names(L; all = true)
                 if startswith(String(n), prefix) && isdefined(L, n) && getfield(L, n) isa Function])
    missing_ = [n for n in raw if !(n in used)]
    println(length(raw), " raw `", prefix, "*` functions; ",
            length(raw) - length(missing_), " reached by the idiomatic layer; ",
            length(missing_), " unreferenced:")
    foreach(n -> println("  ", n), missing_)
end

main(isempty(ARGS) ? "vDSP" : ARGS[1])
