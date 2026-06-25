using AppleAccelerate
using LinearAlgebra
using Test

@testset "BNNS core" begin
    @testset "matmul" begin
        for (m, k, n) in ((4, 3, 5), (1, 1, 1), (8, 8, 8), (10, 1, 7), (2, 6, 3))
            A = randn(Float32, m, k)
            B = randn(Float32, k, n)
            ref = A * B
            C = bnns_matmul(A, B)
            @test size(C) == (m, n)
            @test C ≈ ref atol=1f-4 rtol=1f-4

            # alpha scaling
            @test bnns_matmul(A, B; alpha = 2.5f0) ≈ 2.5f0 .* ref atol=1f-4 rtol=1f-4

            # in-place
            Cout = Matrix{Float32}(undef, m, n)
            @test bnns_matmul!(Cout, A, B) === Cout
            @test Cout ≈ ref atol=1f-4 rtol=1f-4
        end

        # dimension checks
        @test_throws DimensionMismatch bnns_matmul!(zeros(Float32, 2, 2),
                                                    randn(Float32, 2, 3), randn(Float32, 4, 2))
        @test_throws DimensionMismatch bnns_matmul!(zeros(Float32, 3, 3),
                                                    randn(Float32, 2, 3), randn(Float32, 3, 2))
    end

    @testset "activation" begin
        X = randn(Float32, 6, 4) .- 0.5f0
        refs = Dict(
            :identity => X,
            :relu     => max.(X, 0f0),
            :sigmoid  => 1 ./ (1 .+ exp.(-X)),
            :tanh     => tanh.(X),
            :abs      => abs.(X),
        )
        for (f, ref) in refs
            Y = bnns_activation(f, X)
            @test size(Y) == size(X)
            @test Y ≈ ref atol=1f-5 rtol=1f-5

            # in-place, including non-Vector destination round-trip
            Yout = similar(X)
            @test bnns_activation!(f, Yout, X) === Yout
            @test Yout ≈ ref atol=1f-5 rtol=1f-5
        end

        # vector input (the trivial-layout fast path)
        v = randn(Float32, 100) .- 0.5f0
        @test bnns_activation(:relu, v) ≈ max.(v, 0f0) atol=1f-6

        @test_throws ArgumentError bnns_activation(:nope, X)
        @test_throws DimensionMismatch bnns_activation!(:relu, zeros(Float32, 2), X)
    end

    @testset "BNNSArray descriptor" begin
        A = rand(Float32, 3, 4)
        d = BNNSArray(A)
        @test d.parent === A
        @test d.desc.size[1] == 3
        @test d.desc.size[2] == 4
        v = rand(Float32, 5)
        dv = BNNSArray(v)
        @test dv.desc.size[1] == 5
        @test_throws ArgumentError BNNSArray(rand(Float32, 2, 2, 2))
        @test_throws ArgumentError BNNSArray(rand(Float64, 2, 2))
    end
end
