using LinearAlgebra

# Set up a debugging fallback function that prints out a stacktrace if the LinearAlgebra
# tests end up calling a function that we don't have forwarded.
accelerate_header = """
import LinearAlgebra
function debug_missing_function()
    println("Missing BLAS/LAPACK function!")
    display(stacktrace())
end
LinearAlgebra.BLAS.lbt_set_default_func(@cfunction(debug_missing_function, Cvoid, ()))
import AppleAccelerate
"""
eval(Meta.parseall(accelerate_header))

config = BLAS.get_config()
@info("Running with $(length(config.loaded_libs)) libraries loaded:")
display(config.loaded_libs)

using Test
@testset "Sanity Tests" begin
    @test LinearAlgebra.peakflops() > 0
    @test endswith(BLAS.lbt_find_backing_library("dgemm_", :lp64).libname, "Accelerate")
    @test endswith(BLAS.lbt_find_backing_library("dgemm_", :ilp64).libname, "Accelerate")

    # Accelerate has `_rook` symbols:
    @test endswith(BLAS.lbt_find_backing_library("dsytrf_rook_", :ilp64).libname, "Accelerate")
end

@testset "CBLAS dot test" begin
    a = ComplexF64[
        1 + 1im,
        2 - 2im,
        3 + 3im
    ]
    @test BLAS.dotc(a, a) ≈ ComplexF64(28)
    @test BLAS.dotu(a, a) ≈ ComplexF64(12im)

    a = Float32[1, 2, 3]
    @test BLAS.dot(a, a) ≈ 14f0
end

# Run all the LinearAlgebra stdlib tests, but with Accelerate.  We still
# use `Base.runtests()` to get multithreaded, distributed execution
# to cut down on CI times, and also to restart workers that trip over
# the testing RSS limit.  In order for distributed workers to use Accelerate,
# we'll modify the test source code so that it imports Accelerate:
@testset "Full LinearAlgebra test suite" begin; mktempdir() do dir
    cp(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia", "test"), dir; force=true, follow_symlinks=true)

    # Prepend `using AppleAccelerate` to `testdefs.jl`, so that all test workers load Acclerate
    testdefs_path = joinpath(dir, "testdefs.jl")
    chmod(testdefs_path, 0o644)
    testdefs_content = String(read(testdefs_path))
    open(testdefs_path, write=true) do io
        println(io, accelerate_header)
        println(io, testdefs_content)
    end

    run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(dir)/runtests.jl LinearAlgebra`)
end; end
