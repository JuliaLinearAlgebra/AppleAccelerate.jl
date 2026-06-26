# Regenerate src/lib/LibAccelerate.jl from Apple's Accelerate C headers.
#
# Usage:
#   julia --project=gen gen/generate.jl
#
# Requires the macOS SDK headers (Xcode or the Command Line Tools). The headers
# are read only at generation time; the *committed* output has no such dependency,
# so end users need nothing but the runtime framework that ships with macOS.

using Clang.Generators
using Libdl

cd(@__DIR__)

# Resolve the SDK and vecLib header directory from the active toolchain rather
# than hard-coding a path — this tracks whatever Xcode/CLT the developer has.
const SDK = readchomp(`xcrun --show-sdk-path`)
const VECLIB = joinpath(
    SDK,
    "System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers",
)

options = load_options(joinpath(@__DIR__, "generator.toml"))

args = get_default_args()
push!(args, "-isysroot", SDK)
push!(args, "-I", VECLIB)
# Resolve `<CoreFoundation/CFAvailability.h>` and other framework headers that vDSP.h pulls
# in for its availability annotations.
push!(args, "-iframework", joinpath(SDK, "System", "Library", "Frameworks"))

# In-scope headers. Excluded by design:
#   - cblas*.h / blas_new.h / lapack*.h / clapack.h / fortran_blas.h
#     → BLAS/LAPACK are forwarded via libblastrampoline, not ccall.
#   - LinearAlgebra/  → C++ generics, not C-mappable.
#   - Sparse/BLAS.h   → C++ name-mangled dense×sparse multiply (hand-wrapped in sparse.jl).
# We include Sparse/Solve.h (the C solver API) directly rather than Sparse/Sparse.h so we
# don't pull in the C++ BLAS.h. Likewise the BNNS umbrella + graph headers pull in the
# struct/constant headers transitively.
headers = [
    joinpath(VECLIB, "vDSP.h"),
    joinpath(VECLIB, "vForce.h"),
    joinpath(VECLIB, "vBasicOps.h"),
    joinpath(VECLIB, "vfp.h"),
    joinpath(VECLIB, "vectorOps.h"),
    joinpath(VECLIB, "vBigNum.h"),
    joinpath(VECLIB, "Sparse", "Solve.h"),
    joinpath(VECLIB, "BNNS", "bnns.h"),
    # bnns_graph.h is routed through a shim that neutralizes trailing availability
    # attributes which otherwise break Clang.jl's anonymous-struct typedef resolution.
    joinpath(@__DIR__, "shims", "bnns_graph_shim.h"),
    joinpath(VECLIB, "Quadrature", "Quadrature.h"),
]

ctx = create_context(headers, args, options)
build!(ctx)

# Post-process: strip BLAS/LAPACK function wrappers that are transitively pulled in via
# Sparse/Types.h (which includes cblas.h for its index/order types). BLAS and LAPACK are
# forwarded through libblastrampoline, not ccall, so these wrappers are explicitly out of
# scope. We remove only the `function cblas_*`/`catlas_*`/`clapack_*` wrapper blocks; the
# CBLAS enum/constant type definitions are left in place because Sparse references them. This
# runs as part of the generator, so the committed output remains fully reproducible.
function strip_out_of_scope_blas!(path)
    lines = readlines(path)
    out = String[]
    blas = r"^function\s+(cblas_|catlas_|clapack_|appleblas_)"
    i, n, removed = 1, length(lines), 0
    while i <= n
        if occursin(blas, lines[i])
            while i <= n && lines[i] != "end"
                i += 1
            end
            i += 1                                   # skip the closing `end`
            i <= n && isempty(lines[i]) && (i += 1)  # swallow one trailing blank line
            removed += 1
        else
            push!(out, lines[i]); i += 1
        end
    end
    write(path, join(out, "\n") * "\n")
    return removed
end

# Post-process: strip generated wrappers whose `@ccall libacc.<sym>` target is not an
# exported symbol of the Accelerate framework. Clang.jl turns every C declaration it sees
# into a `function` wrapper, but several Accelerate "functions" have no exported symbol and
# would throw `could not load symbol` if ever called:
#   - The high-level Sparse/Dense Solve API (`SparseFactor`, `SparseSolve`, `SparseMultiply`,
#     `_DenseMatrixFromVector_*`, …) is defined entirely as `__attribute__((overloadable))`
#     inline functions/macros in Sparse/Solve.h. The *real* entry points are the
#     underscore-suffixed-by-type variants (`_SparseFactorSymmetric_Double`,
#     `_SparseSolveOpaque_Double`, …), which DO export and are used by src/sparse.jl.
#   - A few header-only BNNS graph setters (`BNNSGraphContextSetBatchSize`, …) are declared
#     but not exported in the shipping framework binary.
#   - `CF_ENUM` is a CoreFoundation macro that the generator misparsed as a function
#     (also covered by the generator.toml ignorelist; this is a belt-and-braces backstop).
# Rather than maintain a hand-curated denylist, we probe each wrapper's exact ccall symbol
# with `dlsym` against the live framework at generation time and drop the block iff the
# symbol does not resolve. This is deterministic on a given machine and self-maintaining:
# if Apple ever exports one of these, it simply stops being stripped. Wrappers whose symbol
# resolves (the vast majority, including all the real `_Sparse*` ones) are left untouched.
function strip_dead_symbol_wrappers!(path)
    h = dlopen("/System/Library/Frameworks/Accelerate.framework/Accelerate")
    lines = readlines(path)
    out = String[]
    funchead = r"^function\s+[A-Za-z_][A-Za-z0-9_]*\("
    ccallsym = r"@ccall\s+libacc\.([A-Za-z_][A-Za-z0-9_]*)\("
    i, n, removed = 1, length(lines), String[]
    while i <= n
        if occursin(funchead, lines[i])
            # Capture the whole `function … end` block, then decide whether to keep it.
            block = String[lines[i]]
            j = i + 1
            while j <= n && lines[j] != "end"
                push!(block, lines[j]); j += 1
            end
            j <= n && push!(block, lines[j])         # the closing `end`
            sym = nothing
            for b in block
                m = match(ccallsym, b)
                if m !== nothing
                    sym = m.captures[1]; break
                end
            end
            if sym !== nothing && dlsym_e(h, sym) == C_NULL
                push!(removed, sym)
                i = j + 1
                i <= n && isempty(lines[i]) && (i += 1)  # swallow one trailing blank line
            else
                append!(out, block)
                i = j + 1
            end
        else
            push!(out, lines[i]); i += 1
        end
    end
    write(path, join(out, "\n") * "\n")
    return removed
end

# Post-process: correct the double-precision complex typedefs. Apple's headers define the
# complex element types via anonymous `_Complex` typedefs (e.g. `typedef _Complex double
# __double_complex_t;`). Clang.jl mis-resolves these anonymous-`_Complex` typedefs and emits
# the double-precision ones as `ComplexF32` instead of `ComplexF64`, which would silently give
# half-width buffers to any double-complex symbol (e.g. `vvcosisin`, `SparseMatrix_Complex_Double`).
# The single-precision typedefs (`__float_complex_t`, `__SPARSE_float_complex`) are correctly
# `ComplexF32` and are left untouched. This is a targeted, idempotent rewrite of exactly the two
# offending `const … = ComplexF32` lines, so the committed output stays deterministic.
function fix_double_complex_typedefs!(path)
    lines = readlines(path)
    fixed = 0
    for name in ("__double_complex_t", "__SPARSE_double_complex")
        wrong = "const $name = ComplexF32"
        right = "const $name = ComplexF64"
        for i in eachindex(lines)
            if lines[i] == wrong
                lines[i] = right
                fixed += 1
            end
        end
    end
    write(path, join(lines, "\n") * "\n")
    return fixed
end

out_path = joinpath(@__DIR__, "..", "src", "lib", "LibAccelerate.jl")
removed = strip_out_of_scope_blas!(out_path)
@info "Stripped $removed out-of-scope BLAS/LAPACK wrappers (forwarded via libblastrampoline)"
dead = strip_dead_symbol_wrappers!(out_path)
@info "Stripped $(length(dead)) wrappers whose ccall symbol does not exist in the framework" dead
fixed = fix_double_complex_typedefs!(out_path)
@info "Corrected $fixed double-precision complex typedef(s) (Clang.jl anonymous-_Complex bug)"
