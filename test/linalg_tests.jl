if AppleAccelerate.get_macos_version() < v"13.4"
    @info("AppleAccelerate.jl needs macOS >= 13.4 for BLAS forwarding. Not testing forwarding capabilities.")
    exit(0)
end

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

@testset "Dense Linear Algebra" begin

@testset "Accelerate Forwarding Sanity Tests" begin
    ver = AppleAccelerate.get_macos_version()
    if ver >= v"13.4"
        @test LinearAlgebra.peakflops() > 0
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :lp64).libname, "Accelerate")
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :ilp64).libname, "Accelerate")

        # Accelerate has `_rook` symbols:
        @test endswith(BLAS.lbt_find_backing_library("dsytrf_rook_", :ilp64).libname, "Accelerate")
    end
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

@testset "BLAS threading tests" begin
    if AppleAccelerate.get_macos_version() >= v"15"
        @test AppleAccelerate.set_num_threads(1) == 1
        @test AppleAccelerate.get_num_threads() == 1
        @test AppleAccelerate.set_num_threads(4) > 1
        @test AppleAccelerate.get_num_threads() > 1
    else
        @test AppleAccelerate.get_num_threads() == 1
        @test AppleAccelerate.set_num_threads(1) == 1
    end
end

linalg_stdlib_test_path = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")

@testset verbose=false "LinearAlgebra.jl BLAS tests" begin
    joinpath(linalg_stdlib_test_path, "blas.jl") |> include
end

@testset verbose=false "LinearAlgebra.jl LAPACK tests" begin
    joinpath(linalg_stdlib_test_path, "lapack.jl") |> include
end

end # @testset "Dense Linear Algebra"
