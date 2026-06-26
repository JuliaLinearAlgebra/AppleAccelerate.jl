using AppleAccelerate
using AppleAccelerate: BNNSArray, bnns_matmul, bnns_matmul!, bnns_activation, bnns_activation!,
                      bnns_data_type
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

        # Non-Array (strided/transposed) inputs go through the `Array(A)` copy
        # path inside bnns_matmul!.
        let Al = randn(Float32, 4, 3), Bl = randn(Float32, 4, 5)
            @test bnns_matmul(transpose(Al), Bl) ≈ transpose(Al) * Bl atol = 1f-4 rtol = 1f-4
            adj = adjoint(randn(Float32, 3, 6))   # 6×3 adjoint wrapper
            Bk = randn(Float32, 3, 2)
            @test bnns_matmul(adj, Bk) ≈ adj * Bk atol = 1f-4 rtol = 1f-4
            @test bnns_matmul(adj, Bk) isa Matrix{Float32}
        end

        # View / transpose destination exercises the temp-buffer copy-back path
        # (Cc !== C) in bnns_matmul!.
        let A = randn(Float32, 4, 3), B = randn(Float32, 3, 5), ref = nothing
            ref = A * B
            Cbig = zeros(Float32, 8, 8)
            Cv = view(Cbig, 1:4, 1:5)
            @test bnns_matmul!(Cv, A, B) === Cv
            @test Cv ≈ ref atol = 1f-4 rtol = 1f-4
            @test all(==(0f0), Cbig[5:8, :])  # untouched region

            sq = randn(Float32, 3, 3)
            Ct = transpose(zeros(Float32, 3, 3))
            @test bnns_matmul!(Ct, sq, sq) === Ct
            @test Ct ≈ sq * sq atol = 1f-4 rtol = 1f-4
        end
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

        # In-place aliasing X === Y on a Vector (Yc === Xc, true in-place).
        let w = randn(Float32, 64) .- 0.5f0, wref = tanh.(copy(w))
            @test bnns_activation!(:tanh, w, w) === w
            @test w ≈ wref atol = 1f-5 rtol = 1f-5
        end

        # Transposed / strided input drives the `vec(Array(X))` copy path.
        let Xs = (randn(Float32, 5, 7) .- 0.5f0)
            Xt = transpose(Xs)              # 7×5 lazy transpose
            @test bnns_activation(:sigmoid, Xt) ≈ (1 ./ (1 .+ exp.(-Xt))) atol = 1f-5 rtol = 1f-5
            Yt = similar(Array(Xt))
            @test bnns_activation!(:relu, Yt, Xt) === Yt
            @test Yt ≈ max.(Array(Xt), 0f0) atol = 1f-6
        end

        @test_throws ArgumentError bnns_activation(:nope, X)
        @test_throws DimensionMismatch bnns_activation!(:relu, zeros(Float32, 2), X)
    end

    @testset "BNNSArray descriptor" begin
        A = rand(Float32, 3, 4)
        d = BNNSArray(A)
        @test d.parent === A
        @test d.desc.size[1] == 3
        @test d.desc.size[2] == 4
        # column-major strides reported explicitly: stride[2] == rows.
        @test d.desc.stride[1] == 1
        @test d.desc.stride[2] == 3
        v = rand(Float32, 5)
        dv = BNNSArray(v)
        @test dv.desc.size[1] == 5
        @test dv.desc.stride[1] == 1

        # Int32 element type goes through bnns_data_type(Int32).
        @test bnns_data_type(Float32) == bnns_data_type(Float32)
        @test bnns_data_type(Int32) != bnns_data_type(Float32)
        @test_throws ArgumentError bnns_data_type(Float64)
        di = BNNSArray(Int32[1 2 3; 4 5 6])
        @test di.desc.size[1] == 2
        @test di.desc.size[2] == 3
        @test di.parent isa Matrix{Int32}
        dvi = BNNSArray(Int32[1, 2, 3, 4])
        @test dvi.desc.size[1] == 4

        @test_throws ArgumentError BNNSArray(rand(Float32, 2, 2, 2))
        @test_throws ArgumentError BNNSArray(rand(Float64, 2, 2))
    end
end
