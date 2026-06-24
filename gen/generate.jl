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

# Proof-of-concept scope: just the Quadrature subframework. Adding BNNS, Sparse,
# vForce, etc. is a matter of appending their headers here.
headers = [
    joinpath(VECLIB, "Quadrature", "Quadrature.h"),
]

ctx = create_context(headers, args, options)
build!(ctx)
