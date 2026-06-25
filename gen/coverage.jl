# Binding-coverage report (informational only — never fails CI).
#
# For each in-scope Accelerate subframework, reports:
#   total    — C functions present in the generated raw layer (src/lib/LibAccelerate.jl)
#   wrapped  — of those, how many are referenced by name in the hand-written idiomatic
#              layer (src/*.jl outside lib/, and ext/*.jl)
#
# Usage:
#   julia --project=. gen/coverage.jl            # print table
#   julia --project=. gen/coverage.jl --write    # also (re)write docs/src/coverage.md
#
# Classification is by source header: each generated function name is attributed to the
# subframework whose header(s) declare it. Note: idiomatic modules that build C symbol names
# dynamically (e.g. `Symbol("vv", fa, suff)`) are undercounted by the static name scan; this
# is a floor on real coverage, called out in the generated page.

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SDK = readchomp(`xcrun --show-sdk-path`)
const VECLIB = joinpath(SDK,
    "System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers")

# Subframework → header files whose declarations belong to it. Order matters: the first group
# whose headers declare a symbol claims it (disambiguates the vU*/vS* overlap between
# vBasicOps and vBigNum by header membership).
header_groups() = [
    "vForce"     => [joinpath(VECLIB, "vForce.h")],
    "vBigNum"    => [joinpath(VECLIB, "vBigNum.h")],
    "vBasicOps"  => [joinpath(VECLIB, "vBasicOps.h")],
    "vfp"        => [joinpath(VECLIB, "vfp.h")],
    "vectorOps"  => [joinpath(VECLIB, "vectorOps.h")],
    "vDSP"       => filter(isfile, [joinpath(VECLIB, "vDSP.h"),
                                    joinpath(VECLIB, "vDSP_translate.h")]),
    "Sparse"     => filter(f -> isfile(f) && !endswith(f, "BLAS.h"),
                           readdir(joinpath(VECLIB, "Sparse"); join=true)),
    "BNNS"       => filter(isfile, readdir(joinpath(VECLIB, "BNNS"); join=true)),
    "Quadrature" => filter(isfile, readdir(joinpath(VECLIB, "Quadrature"); join=true)),
]

# All function names defined in the generated raw layer.
function generated_functions()
    text = read(joinpath(ROOT, "src", "lib", "LibAccelerate.jl"), String)
    names = String[]
    for m in eachmatch(r"^function\s+([A-Za-z_][A-Za-z0-9_]*)\("m, text)
        push!(names, m.captures[1])
    end
    return unique(names)
end

# Concatenated text of the hand-written idiomatic layer (everything that consumes the raw
# bindings): all of src/ except the generated file, plus ext/.
function idiomatic_text()
    bufs = String[]
    for dir in (joinpath(ROOT, "src"), joinpath(ROOT, "ext"))
        isdir(dir) || continue
        for (root, _, files) in walkdir(dir)
            occursin(joinpath("src", "lib"), root) && continue
            for f in files
                endswith(f, ".jl") || continue
                push!(bufs, read(joinpath(root, f), String))
            end
        end
    end
    return join(bufs, "\n")
end

# header text cache
function group_of(name, groups, cache)
    # The vBasicOps/vBigNum integer ops are macro-generated (token-pasted widths), so their
    # names never appear literally in the headers — classify them by name pattern. Widths
    # 128/256/512/1024 are vBigNum; narrower widths are vBasicOps.
    m = match(r"^v[US](\d+)", name)
    if m !== nothing
        return parse(Int, m.captures[1]) >= 128 ? "vBigNum" : "vBasicOps"
    end
    # Macro-generated, token-pasted symbols that don't appear literally in their headers.
    occursin(r"^_?Sparse|^sparse_|^_Dense", name) && return "Sparse"
    occursin(r"^_?BNNS|^bnns", name) && return "BNNS"
    # Everything else: attribute by header membership — the first group whose headers mention
    # the symbol name (as a whole word). Headers don't cross-mention each other's symbols, so
    # collisions are negligible.
    pat = Regex("\\b" * name * "\\b")
    for (grp, files) in groups
        for f in files
            txt = get!(cache, f) do
                isfile(f) ? read(f, String) : ""
            end
            occursin(pat, txt) && return grp
        end
    end
    return "other"
end

function report()
    funcs = generated_functions()
    groups = header_groups()
    cache = Dict{String,String}()
    idiom = idiomatic_text()

    total = Dict{String,Int}()
    wrapped = Dict{String,Int}()
    for g in first.(groups); total[g] = 0; wrapped[g] = 0; end
    total["other"] = 0; wrapped["other"] = 0

    for fn in funcs
        g = group_of(fn, groups, cache)
        total[g] += 1
        # exact-name reference in the idiomatic layer
        if occursin(Regex("\\b" * fn * "\\b"), idiom)
            wrapped[g] += 1
        end
    end

    order = vcat(first.(groups), "other")
    return order, total, wrapped
end

function format_table(order, total, wrapped)
    io = IOBuffer()
    println(io, "| Subframework | Wrapped | Total | Coverage |")
    println(io, "|--------------|--------:|------:|---------:|")
    tw = ttot = 0
    for g in order
        total[g] == 0 && continue
        pct = total[g] == 0 ? 0 : round(Int, 100 * wrapped[g] / total[g])
        println(io, "| $g | $(wrapped[g]) | $(total[g]) | $(pct)% |")
        tw += wrapped[g]; ttot += total[g]
    end
    pct = ttot == 0 ? 0 : round(Int, 100 * tw / ttot)
    println(io, "| **Total** | **$tw** | **$ttot** | **$(pct)%** |")
    return String(take!(io))
end

function main()
    order, total, wrapped = report()
    table = format_table(order, total, wrapped)
    println("\nLibAccelerate binding coverage (idiomatic wrappers / generated raw functions):\n")
    println(table)

    if "--write" in ARGS
        page = """
        # Binding Coverage

        How much of each Accelerate subframework has a hand-written idiomatic Julia API on top
        of the generated `LibAccelerate` raw bindings. **Total** counts the C functions present
        in the generated raw layer (which itself mirrors the in-scope headers); **Wrapped**
        counts those referenced by name in the idiomatic layer.

        !!! note
            Modules that construct C symbol names dynamically (e.g. `Symbol("vv", fa, suff)`)
            are undercounted by the static scan, so these numbers are a floor on real coverage.
            BLAS/LAPACK are intentionally excluded — they are forwarded via libblastrampoline,
            not ccall.

        $table
        """
        out = joinpath(ROOT, "docs", "src", "coverage.md")
        write(out, page)
        println("\nWrote $out")
    end
end

main()
