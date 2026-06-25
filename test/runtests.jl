using AppleAccelerate
using Test

if !Sys.isapple()
    @info("AppleAccelerate.jl will be tested only on macOS. Exiting.")
    exit(0)
end

# Extension-loading / package-contract tests run FIRST, in fresh subprocesses,
# so they can observe the bare `using AppleAccelerate` state (no LinearAlgebra /
# SparseArrays, no BLAS forwarding) before this process loads the weakdeps below.
include("extension_loading_tests.jl")

# Load the weakdep triggers so the package extensions (LinearAlgebra,
# SparseArrays, and the dual LinearAlgebra+SparseArrays extension) activate for
# the remaining test files, which use sparse/factorization/forwarding freely.
using LinearAlgebra
using SparseArrays
using AbstractFFTs
using DSP, FFTW, Random, Statistics
using Aqua

@testset "Aqua" begin
    Aqua.test_all(AppleAccelerate)
end

Random.seed!(7)
N = 1_000

include("lib_tests.jl")
include("array_tests.jl")
include("vmath_tests.jl")
include("complexarray_tests.jl")
include("dsp_tests.jl")
include("sparse_tests.jl")
include("quadrature_tests.jl")
include("bnns_tests.jl")
include("nnlib_ext_tests.jl")
# linalg_tests.jl runs LAST on purpose: it exercises BLAS/LAPACK forwarding
# (macOS >= 13.4) and cleanly `exit(0)`s on older systems. Keeping it last means
# that early exit can never skip the non-forwarding suites above.
include("linalg_tests.jl")
