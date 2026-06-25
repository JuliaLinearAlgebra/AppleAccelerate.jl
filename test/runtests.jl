using LinearAlgebra
using AppleAccelerate
using AbstractFFTs
using DSP, FFTW, Random, Statistics, Test
using Aqua

if !Sys.isapple()
    @info("AppleAccelerate.jl will be tested only on macOS. Exiting.")
    exit(0)
end

@testset "Aqua" begin
    Aqua.test_all(AppleAccelerate)
end

Random.seed!(7)
N = 1_000

include("array_tests.jl")
include("complexarray_tests.jl")
include("dsp_tests.jl")
include("sparse_tests.jl")
include("linearalgebra_ext_tests.jl")
include("linalg_tests.jl")
include("quadrature_tests.jl")
