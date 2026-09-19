# Detect drift between the committed src/lib/LibAccelerate.jl and what the generator
# produces from the macOS SDK on this machine.
#
# Usage:
#   julia gen/check_drift.jl                       # regenerate into a temp tree, then compare
#   julia gen/check_drift.jl --generated <file>    # compare against an already-generated file
#   julia gen/check_drift.jl --committed <file>    # override the committed side
#   julia gen/check_drift.jl --summary <file>      # also write the Markdown report to <file>
#
# Exit status: 0 = no drift (informational differences may still be listed),
#              1 = drift (a shared definition changed), 2 = the generator itself failed.
#
# A byte-for-byte diff of the generated file is useless across machines (see README.md,
# "Reproducibility"): Clang.jl renumbers the `var"##Ctag#NNN"` anonymous types wholesale
# whenever the SDK adds or removes one, and the dead-symbol strip pass depends on which
# symbols the *running* macOS exports. So this compares the stable surface instead:
#
#   * functions      — by name; the `@ccall` signature (symbol, argument types, return type)
#   * named structs  — by name; field names and types
#   * named enums    — by name; base type, member names and values
#   * constants      — by name; value
#   * anonymous (`##Ctag#`) types — as a multiset of bodies, with the counter erased
#
# A name present on only one side is *informational* (a newer SDK adds API; a different
# macOS strips a different set of unexported symbols). A name present on both sides whose
# definition differs is *drift* and fails: the committed bindings no longer match the ABI
# the headers describe. An anonymous body that exists only in the committed file is also
# drift, since a changed anonymous type is indistinguishable from a removed one.
#
# Needs only Base; the regeneration step shells out to `julia --project=gen`.

const GEN_DIR = @__DIR__
const REPO = dirname(GEN_DIR)
const COMMITTED_DEFAULT = joinpath(REPO, "src", "lib", "LibAccelerate.jl")

const CTAG = r"var\"##Ctag#\d+\""
const ANON = "<anon>"

normalize(s::AbstractString) = replace(strip(s), CTAG => ANON, r"\s+" => " ")

struct Surface
    named::Dict{String,String}       # "kind name" => normalized definition
    anon::Dict{String,Int}           # normalized anonymous definition => count
end

"""
    parse_surface(path) -> Surface

Line-oriented parse of a Clang.jl-generated module. Relies on the generator's fixed
layout: top-level `function`/`struct`/`mutable struct`/`@enum` blocks open at column 0
and close with a bare `end` at column 0; constants are single `const` lines.
"""
function parse_surface(path::AbstractString)
    lines = readlines(path)
    named = Dict{String,String}()
    anon = Dict{String,Int}()
    blockhead = r"^(function|struct|mutable struct|@enum)\s+(.*)$"
    i, n = 1, length(lines)
    while i <= n
        line = lines[i]
        m = match(blockhead, line)
        if m !== nothing
            kind, head = m.captures[1], m.captures[2]
            body = String[]
            j = i + 1
            while j <= n && lines[j] != "end"
                isempty(strip(lines[j])) || push!(body, normalize(lines[j]))
                j += 1
            end
            record_block!(named, anon, String(kind), String(head), body)
            i = j + 1
        elseif startswith(line, "const ")
            c = match(r"^const\s+(\S+)\s*=\s*(.*)$", line)
            if c !== nothing
                name, val = c.captures[1], normalize(c.captures[2])
                if occursin(CTAG, name)
                    key = "const $ANON = $val"
                    anon[key] = get(anon, key, 0) + 1
                else
                    named["const $name"] = val
                end
            end
            i += 1
        else
            i += 1
        end
    end
    return Surface(named, anon)
end

function record_block!(named, anon, kind::String, head::String, body::Vector{String})
    if occursin(CTAG, head)
        # Anonymous type (or an accessor method on one): identity is the body alone.
        key = string(kind, " ", normalize(head), " { ", join(body, "; "), " }")
        anon[key] = get(anon, key, 0) + 1
        return
    end
    if kind == "function"
        fm = match(r"^([A-Za-z_][A-Za-z0-9_]*)\(", head)
        if fm !== nothing
            # Plain wrapper: the ABI is the @ccall line(s); Julia argument names are not.
            ccalls = filter(l -> occursin("@ccall", l), body)
            sig = isempty(ccalls) ? join(body, "; ") : join(ccalls, "; ")
            named["function $(fm.captures[1])"] = sig
        else
            # `Base.getproperty(x::Ptr{T}, f::Symbol)` & co. — union/bitfield accessors
            # whose body encodes field offsets. Keyed by the full method head.
            named["method $(normalize(head))"] = join(body, "; ")
        end
    elseif kind == "@enum"
        em = match(r"^(\S+?)(::\S+)?\s+begin\s*$", head)
        name = em === nothing ? normalize(head) : em.captures[1]
        base = em === nothing || em.captures[2] === nothing ? "" : em.captures[2]
        named["enum $name"] = string(base, " { ", join(body, "; "), " }")
    else
        named["$kind $(normalize(head))"] = join(body, "; ")
    end
    return
end

struct Drift
    changed::Vector{Tuple{String,String,String}}   # (key, committed, generated)
    anon_lost::Vector{String}                      # anonymous bodies only in committed
    only_committed::Vector{String}
    only_generated::Vector{String}
    anon_new::Vector{String}
end

function compare(committed::Surface, generated::Surface)
    changed = Tuple{String,String,String}[]
    only_c, only_g = String[], String[]
    for (k, v) in committed.named
        if haskey(generated.named, k)
            g = generated.named[k]
            g == v || push!(changed, (k, v, g))
        else
            push!(only_c, k)
        end
    end
    for k in keys(generated.named)
        haskey(committed.named, k) || push!(only_g, k)
    end
    anon_lost, anon_new = String[], String[]
    for (k, c) in committed.anon
        d = c - get(generated.anon, k, 0)
        d > 0 && append!(anon_lost, fill(k, d))
    end
    for (k, c) in generated.anon
        d = c - get(committed.anon, k, 0)
        d > 0 && append!(anon_new, fill(k, d))
    end
    return Drift(sort!(changed), sort!(anon_lost), sort!(only_c), sort!(only_g), sort!(anon_new))
end

isdrift(d::Drift) = !isempty(d.changed) || !isempty(d.anon_lost)

function report(io::IO, d::Drift; context::AbstractString = "")
    println(io, "# LibAccelerate generator drift report\n")
    isempty(context) || println(io, context, "\n")
    if isdrift(d)
        println(io, "**Result: DRIFT** — the committed bindings disagree with the headers.\n")
    else
        println(io, "**Result: no drift** in the shared surface.\n")
    end
    if !isempty(d.changed)
        println(io, "## Changed definitions ($(length(d.changed))) — failure\n")
        for (k, c, g) in d.changed
            println(io, "- `$k`\n  - committed: `$c`\n  - generated: `$g`")
        end
        println(io)
    end
    if !isempty(d.anon_lost)
        println(io, "## Anonymous types changed or removed ($(length(d.anon_lost))) — failure\n")
        foreach(k -> println(io, "- `$k`"), d.anon_lost)
        println(io)
    end
    for (title, items) in (
            ("Only in the generated output (new in this SDK / exported by this macOS)", d.only_generated),
            ("Only in the committed file (absent from this SDK / not exported by this macOS)", d.only_committed),
            ("Anonymous types only in the generated output", d.anon_new))
        isempty(items) && continue
        println(io, "## $title ($(length(items))) — informational\n")
        foreach(k -> println(io, "- `$k`"), items)
        println(io)
    end
    return
end

"""
    regenerate() -> path

Run the unmodified generator against a throwaway copy of `gen/` so the committed
`src/lib/LibAccelerate.jl` is never touched. `generate.jl` resolves every path relative
to its own directory, so a copied `gen/` next to an empty `src/lib/` is all it needs; the
real `gen/` project is still the active environment, so Clang.jl stays pinned.
"""
function regenerate()
    tmp = mktempdir()
    cp(GEN_DIR, joinpath(tmp, "gen"))
    mkpath(joinpath(tmp, "src", "lib"))
    julia = Base.julia_cmd()
    run(`$julia --project=$GEN_DIR -e 'using Pkg; Pkg.instantiate()'`)
    run(`$julia --project=$GEN_DIR $(joinpath(tmp, "gen", "generate.jl"))`)
    out = joinpath(tmp, "src", "lib", "LibAccelerate.jl")
    isfile(out) || error("generator finished but did not write $out")
    return out
end

function environment_context()
    sdk = try readchomp(`xcrun --show-sdk-version`) catch; "unknown" end
    os = try readchomp(`sw_vers -productVersion`) catch; "unknown" end
    return "Generated with macOS SDK $sdk on macOS $os (Julia $VERSION)."
end

function main(args)
    committed, generated, summary = COMMITTED_DEFAULT, nothing, nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--committed", "--generated", "--summary") && i < length(args)
            v = args[i + 1]
            a == "--committed" ? (committed = v) : a == "--generated" ? (generated = v) : (summary = v)
            i += 2
        else
            println(stderr, "usage: julia gen/check_drift.jl [--generated FILE] [--committed FILE] [--summary FILE]")
            return 2
        end
    end
    context = ""
    if generated === nothing
        try
            generated = regenerate()
        catch err
            msg = "# LibAccelerate generator drift report\n\n**Result: GENERATOR FAILED** — " *
                  "`gen/generate.jl` did not complete against this SDK.\n\n```\n" *
                  sprint(showerror, err) * "\n```\n\n" * environment_context() * "\n"
            print(msg)
            summary === nothing || write(summary, msg)
            return 2
        end
        context = environment_context()
    end
    d = compare(parse_surface(committed), parse_surface(generated))
    text = sprint(io -> report(io, d; context))
    print(text)
    summary === nothing || write(summary, text)
    return isdrift(d) ? 1 : 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
