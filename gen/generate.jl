# Regenerate src/lib/LibAccelerate.jl from Apple's Accelerate C headers.
#
# Usage:
#   julia --project=gen gen/generate.jl
#
# Requires the macOS SDK headers (Xcode or the Command Line Tools). The headers
# are read only at generation time; the *committed* output has no such dependency,
# so end users need nothing but the runtime framework that ships with macOS.

using Clang.Generators

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
