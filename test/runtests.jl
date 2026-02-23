using LinearAlgebra
using AppleAccelerate
using AbstractFFTs
using DSP, FFTW, Test, Random, Statistics
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

@testset "AppleAccelerate.jl" begin
    include("array_tests.jl")
    include("dsp_tests.jl")
end

include("sparse_tests.jl")
include("linalg_tests.jl")
