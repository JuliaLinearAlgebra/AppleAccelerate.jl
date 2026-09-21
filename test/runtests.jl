using LinearAlgebra
using AppleAccelerate
using DSP, FFTW, Random, Statistics, Test
using Aqua

# The package exports nothing on any platform, and off macOS it is a no-op: it loads, but
# none of the subsystems (whose `include`s are gated on `Sys.isapple()`) are defined. One
# representative name per subsystem file; on macOS the same list is checked the other way
# round so a rename can't leave the non-Apple assertion vacuously true.
@testset "Namespace contract" begin
    @test names(AppleAccelerate) == [:AppleAccelerate]
    subsystem_names = [:LibAccelerate, :exp, :VMATH_COVERAGE, :SIMDMath, :vmags, :fft,
                       :AASparseMatrix, :integrate, :bnns_reduce, :scale_PlanarF]
    defined = names(AppleAccelerate; all = true, imported = false)
    for name in subsystem_names
        @test (name in defined) == Sys.isapple()
    end
end

if !Sys.isapple()
    @info("AppleAccelerate.jl will be tested only on macOS. Exiting.")
    exit(0)
end

@testset "Aqua" begin
    Aqua.test_all(AppleAccelerate)
end

Random.seed!(7)
N = 1_000

include("lib_tests.jl")
include("array_tests.jl")
include("vmath_tests.jl")
include("simdmath_tests.jl")
include("complexarray_tests.jl")
include("edgecase_tests.jl")
include("dsp_tests.jl")
include("fftplan_tests.jl")
include("sparse_tests.jl")
include("quadrature_tests.jl")
include("bnns_tests.jl")
include("vimage_tests.jl")
include("readme_tests.jl")
# linalg_tests.jl runs LAST on purpose: it exercises BLAS/LAPACK forwarding
# (macOS >= 13.4) and cleanly `exit(0)`s on older systems. Keeping it last means
# that early exit can never skip the non-forwarding suites above.
include("linalg_tests.jl")
