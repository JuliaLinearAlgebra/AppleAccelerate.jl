# Tests for AppleAccelerate's LinearAlgebra integration (vDSP-backed matmul /
# transpose). See src/linearalgebra.jl for the architectural note on why this
# lives in the main module rather than a weakdep extension.
using LinearAlgebra

@testset "LinearAlgebra integration" begin

    @testset "accelerate_mul! / accelerate_mul ($T)" for T in (Float32, Float64)
        A = rand(T, 6, 4)
        B = rand(T, 4, 5)
        C = Matrix{T}(undef, 6, 5)
        ref = A * B
        @test AppleAccelerate.accelerate_mul!(C, A, B) === C
        @test C ≈ ref
        @test AppleAccelerate.accelerate_mul(A, B) ≈ ref

        # square case
        A2 = rand(T, 5, 5)
        B2 = rand(T, 5, 5)
        @test AppleAccelerate.accelerate_mul(A2, B2) ≈ A2 * B2

        # dimension mismatch propagates from the vDSP wrapper
        @test_throws DimensionMismatch AppleAccelerate.accelerate_mul!(Matrix{T}(undef, 6, 5), A, rand(T, 3, 5))
    end

    @testset "accelerate_transpose! / accelerate_transpose ($T)" for T in (Float32, Float64)
        A = rand(T, 6, 4)
        C = Matrix{T}(undef, 4, 6)
        @test AppleAccelerate.accelerate_transpose!(C, A) === C
        @test C == permutedims(A)
        @test AppleAccelerate.accelerate_transpose(A) == permutedims(A)

        # destination shape is validated
        @test_throws DimensionMismatch AppleAccelerate.accelerate_transpose!(Matrix{T}(undef, 6, 4), A)
    end

    @testset "consistency with LinearAlgebra/BLAS ($T)" for T in (Float32, Float64)
        # The ordinary LinearAlgebra path is BLAS-forwarded to Accelerate and must
        # agree with the direct vDSP path.
        A = rand(T, 8, 8)
        B = rand(T, 8, 8)
        @test AppleAccelerate.accelerate_mul(A, B) ≈ mul!(similar(A), A, B)
        @test AppleAccelerate.accelerate_transpose(A) ≈ permutedims(A)
    end

end
