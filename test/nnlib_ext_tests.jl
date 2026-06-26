using AppleAccelerate
using NNlib
using Test

@testset "NNlib extension" begin
    # The extension must load once NNlib is present.
    @test Base.get_extension(AppleAccelerate, :AppleAccelerateNNlibExt) !== nothing

    @testset "batched_mul" begin
        for (m, k, n, b) in ((3, 4, 5, 2), (8, 8, 8, 4), (1, 3, 1, 3))
            A = randn(Float32, m, k, b)
            B = randn(Float32, k, n, b)
            ref = NNlib.batched_mul(A, B)
            C = similar(ref)
            NNlib.batched_mul!(C, A, B)
            @test C ≈ ref atol=1f-4 rtol=1f-4

            # alpha scaling through the 5-arg path (β == 0)
            C2 = similar(ref)
            NNlib.batched_mul!(C2, A, B, 2.0f0, 0.0f0)
            @test C2 ≈ 2.0f0 .* ref atol=1f-4 rtol=1f-4

            # β != 0 must fall back and still be correct
            C3 = ones(Float32, m, n, b)
            base = copy(C3)
            NNlib.batched_mul!(C3, A, B, 1.0f0, 0.5f0)
            @test C3 ≈ ref .+ 0.5f0 .* base atol=1f-4 rtol=1f-4
        end

        # Unequal batch dims (broadcast batch) must take the @invoke fallback,
        # not the BNNS fast path (which requires equal size(.,3)).
        let A = randn(Float32, 3, 4, 5), B = randn(Float32, 4, 2, 1)
            ref = NNlib.batched_mul(A, B)
            C = similar(ref)
            @test NNlib.batched_mul!(C, A, B) === C
            @test C ≈ ref atol = 1f-4 rtol = 1f-4
            @test size(C, 3) == 5
        end

        # Transposed batch slices: BatchedTranspose isn't Array{Float32,3}, so
        # this dispatches to the generic NNlib method (sanity cross-check).
        let A = randn(Float32, 4, 3, 2), B = randn(Float32, 4, 5, 2)
            ref = NNlib.batched_mul(NNlib.batched_transpose(A), B)
            @test NNlib.batched_mul(NNlib.batched_transpose(A), B) ≈ ref atol = 1f-4 rtol = 1f-4
        end
    end

    @testset "activations" begin
        X = randn(Float32, 7, 5) .- 0.5f0
        @test AppleAccelerate.bnns_act(relu, X)    ≈ relu.(X)    atol=1f-5
        @test AppleAccelerate.bnns_act(sigmoid, X) ≈ sigmoid.(X) atol=1f-5
        @test AppleAccelerate.bnns_act(tanh, X)    ≈ tanh.(X)    atol=1f-5
        @test AppleAccelerate.bnns_act(abs, X)     ≈ abs.(X)     atol=1f-6
        # σ / tanh_fast / identity alias entries in the _NNLIB_ACT map.
        @test AppleAccelerate.bnns_act(NNlib.σ, X)         ≈ NNlib.σ.(X)         atol=1f-5
        @test AppleAccelerate.bnns_act(NNlib.tanh_fast, X) ≈ NNlib.tanh_fast.(X) atol=1f-5
        @test AppleAccelerate.bnns_act(identity, X)        == X
        @test AppleAccelerate.bnns_act(identity, X) !== X  # returns a fresh array
        # unsupported activation falls back to broadcasting f
        @test AppleAccelerate.bnns_act(NNlib.softplus, X) ≈ NNlib.softplus.(X) atol=1f-5
    end
end
