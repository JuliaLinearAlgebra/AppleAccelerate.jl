# Accelerate ships as a macOS *system framework*. Unlike the usual Clang.jl + JLL
# flow, there is no artifact to load — every Mac already has it. We ccall it by
# absolute path, mirroring the `libacc` const in the top-level AppleAccelerate module.
const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
