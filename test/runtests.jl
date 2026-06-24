using LinearAlgebra
using AppleAccelerate
using AbstractFFTs
using DSP, FFTW, Random, Statistics, Test
using Aqua

@testset "Aqua" begin
    Aqua.test_all(AppleAccelerate)
end

if !Sys.isapple()
    @info("AppleAccelerate.jl will be tested only on macOS. Exiting.")
    exit(0)
end

Random.seed!(7)
N = 1_000

include("array_tests.jl")
include("complexarray_tests.jl")
include("dsp_tests.jl")
include("sparse_tests.jl")
include("linalg_tests.jl")
