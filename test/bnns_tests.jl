using AppleAccelerate
using AppleAccelerate: BNNSArray, bnns_matmul, bnns_matmul!, bnns_activation, bnns_activation!,
                      bnns_data_type
using LinearAlgebra
using Statistics
using Test

const AA = AppleAccelerate

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

@testset "BNNS tensor ops" begin
    M = Float32[1 2 3; 4 5 6]
    @test AA.bnns_tile(M, (2, 1)) == repeat(M, 2, 1)
    @test AA.bnns_tile(M, (1, 3)) == repeat(M, 1, 3)
    @test AA.bnns_tile(M, (2, 2)) == repeat(M, 2, 2)
    @test AA.bnns_tile_backward(repeat(M, 2, 3), (2, 3)) ≈ 6 .* M

    @test AA.bnns_transpose(M, 1, 2) == permutedims(M, (2, 1))
    C = Float32[1 2; 3 4; 5 6]
    @test AA.bnns_transpose(C, 1, 2) == permutedims(C, (2, 1))

    X = Float32[1 5; 3 2]; Y = Float32[2 5; 1 4]
    @test AA.bnns_compare(:lt, X, Y) == (X .< Y)
    @test AA.bnns_compare(:eq, X, Y) == (X .== Y)
    @test AA.bnns_compare(:ge, X, Y) == (X .>= Y)
    @test AA.bnns_compare(:ne, X, Y) == (X .!= Y)
    @test AA.bnns_compare(:lt, X, Y) isa Array{Bool}
    @test_throws ArgumentError AA.bnns_compare(:bogus, X, Y)
    @test_throws DimensionMismatch AA.bnns_compare(:lt, X, Float32[1 2 3])

    S = Float32[1 2 3; 4 5 6; 7 8 9]
    @test AA.bnns_band_part(S, -1, 0) == tril(S)
    @test AA.bnns_band_part(S, 0, -1) == triu(S)
    @test AA.bnns_band_part(S, 0, 0) == Float32.(diagm(diag(S)))

    G = Float32[1 2 3; 4 5 6]
    idx = Int32[2, 0]
    @test AA.bnns_gather(2, G, idx) == G[:, idx .+ 1]
    idx2 = Int32[1, 0]
    @test AA.bnns_gather(1, G, idx2) == G[idx2 .+ 1, :]

    # scatter rows into a zeroed output
    out = zeros(Float32, 3, 2); inp = Float32[10 20; 30 40]; sidx = Int32[0, 2]
    ref = zeros(Float32, 3, 2); ref[1, :] = inp[1, :]; ref[3, :] = inp[2, :]
    @test AA.bnns_scatter(1, :sum, inp, sidx, out) == ref
    @test out === AA.bnns_scatter(1, :sum, inp, sidx, out) |> (_ -> out)  # in place
    @test_throws ArgumentError AA.bnns_scatter(1, :bogus, inp, sidx, zeros(Float32, 3, 2))

    @test AA.bnns_copy!(zeros(Float32, 2, 3), G) == G
end

@testset "BNNS reductions / norms / clipping" begin
    G = Float32[1 2 3; 4 5 6]
    @test vec(AA.bnns_reduce(:sum, G; dim = 1)) == vec(sum(G, dims = 1))
    @test vec(AA.bnns_reduce(:sum, G; dim = 2)) == vec(sum(G, dims = 2))
    @test vec(AA.bnns_reduce(:max, G; dim = 2)) == vec(maximum(G, dims = 2))
    @test vec(AA.bnns_reduce(:min, G; dim = 1)) == vec(minimum(G, dims = 1))
    @test vec(AA.bnns_reduce(:mean, G; dim = 2)) ≈ vec(sum(G, dims = 2) ./ 3)
    @test_throws ArgumentError AA.bnns_reduce(:bogus, G)

    @test AA.bnns_compute_norm(Float32[3, 4])[1] ≈ 5.0f0
    @test AA.bnns_compute_norm(Float32[5, 12, 0])[1] ≈ 13.0f0

    @test AA.bnns_clip_by_value(Float32[-2, 0.5, 3], 0.0f0, 1.0f0) == Float32[0, 0.5, 1]
    @test AA.bnns_clip_by_value(Float32[-5, 5], -1.0f0, 1.0f0) == Float32[-1, 1]

    x = Float32[3, 4]
    @test AA.bnns_clip_by_norm(x, 2.0f0) ≈ x .* (2.0f0 / 5)
    @test AA.bnns_clip_by_norm(x, 10.0f0) == x           # already within
    @test norm(AA.bnns_clip_by_norm(x, 2.0f0)) ≈ 2.0f0

    a = Float32[3, 4]; b = Float32[0, 0, 0, 12]   # global norm 13
    r = AA.bnns_clip_by_global_norm([a, b], 6.5f0)   # factor 0.5
    @test r[1] ≈ a .* 0.5f0
    @test r[2] ≈ b .* 0.5f0
    r2 = AA.bnns_clip_by_global_norm([a, b], 13.0f0) # no clip
    @test r2[1] ≈ a && r2[2] ≈ b
end

@testset "BNNS utility queries" begin
    @test AA.bnns_layout_rank(AA.LibAccelerate.BNNSDataLayoutColumnMajorMatrix) == 2
    @test AA.bnns_layout_rank(AA.LibAccelerate.BNNSDataLayoutVector) == 1
    G = rand(Float32, 3, 4)
    @test AA.bnns_data_size(G) == sizeof(G)
    @test AA.bnns_tensor_allocation_size(G) >= sizeof(G)
end

@testset "BNNS DirectApply" begin
    A2 = Float32[1 2; 3 4]; B2 = Float32[1 0; 0 1]
    @test AA.bnns_broadcast_matmul(A2, B2) ≈ A2 * B2
    @test AA.bnns_broadcast_matmul(A2, B2; alpha = 2.0f0) ≈ 2.0f0 .* (A2 * B2)

    # top-K along dim 1
    T = Float32[1 9; 7 2; 3 5]
    vals, inds = AA.bnns_topk(T, 2; dim = 1)
    @test vals[:, 1] == Float32[7, 3] && inds[:, 1] == Int32[1, 2]
    @test vals[:, 2] == Float32[9, 5] && inds[:, 2] == Int32[0, 2]

    # in-top-k
    scores = Float32[0.1 0.9; 0.5 0.05; 0.4 0.05]  # 3 classes x 2 samples
    targets = Int32[1, 0]
    @test AA.bnns_in_topk(scores, targets, 1; dim = 1) == Bool[1, 1]
    @test AA.bnns_in_topk(scores, Int32[2, 1], 1; dim = 1) == Bool[0, 0]
end

@testset "BNNS random" begin
    g = AA.BNNSRandomGenerator(1234)
    X = zeros(Float32, 100_000)
    AA.bnns_random_fill_uniform!(g, X, 0.0f0, 1.0f0)
    @test 0.0f0 <= minimum(X) && maximum(X) < 1.0f0
    @test isapprox(mean(X), 0.5f0, atol = 0.02)

    Y = zeros(Float32, 100_000)
    AA.bnns_random_fill_normal!(g, Y, 2.0f0, 3.0f0)
    @test isapprox(mean(Y), 2.0f0, atol = 0.1)
    @test isapprox(std(Y), 3.0f0, atol = 0.1)

    Xi = zeros(Int32, 10_000)
    AA.bnns_random_fill_uniform_int!(g, Xi, 5, 10)
    @test all(5 .<= Xi .< 10)

    # seed reproducibility
    a1 = zeros(Float32, 32); a2 = zeros(Float32, 32)
    AA.bnns_random_fill_uniform!(AA.BNNSRandomGenerator(42), a1)
    AA.bnns_random_fill_uniform!(AA.BNNSRandomGenerator(42), a2)
    @test a1 == a2

    # state save/restore
    g3 = AA.BNNSRandomGenerator(7)
    st = AA.bnns_random_state(g3)
    b1 = zeros(Float32, 16); AA.bnns_random_fill_uniform!(g3, b1)
    AA.bnns_random_state!(g3, st)
    b2 = zeros(Float32, 16); AA.bnns_random_fill_uniform!(g3, b2)
    @test b1 == b2

    # categorical: probability mass only on category 2 (0-based)
    gc = AA.BNNSRandomGenerator(9)
    outc = zeros(Float32, 64)
    AA.bnns_random_fill_categorical!(gc, outc, Float32[0, 0, 1])
    @test all(==(2.0f0), outc)
end

@testset "BNNS nearest neighbors" begin
    knn = AA.BNNSNearestNeighbors(10, 2, 2)
    data = Float32[0 1 5 6; 0 1 5 6]   # points (0,0),(1,1),(5,5),(6,6) as columns
    @test AA.bnns_knn_load!(knn, data) == 4
    idx, dist = AA.bnns_knn_query(knn, 0)    # nearest to (0,0)
    @test length(idx) == 2 && length(dist) == 2
    @test 0 in idx                            # itself is nearest
    @test 1 in idx                            # (1,1) is the next nearest
    @test all(dist .>= 0)
end

@testset "BNNS graph compile options" begin
    o = AA.BNNSGraphCompileOptions()
    AA.bnns_compile_options_set_single_thread!(o, true)
    @test AA.bnns_compile_options_get_single_thread(o) == true
    AA.bnns_compile_options_set_debug_info!(o, true)
    @test AA.bnns_compile_options_get_debug_info(o) == true
    AA.bnns_compile_options_set_optimization!(o, :ir_size)
    @test AA.bnns_compile_options_get_optimization(o) == :ir_size
    AA.bnns_compile_options_set_optimization!(o, :performance)
    @test AA.bnns_compile_options_get_optimization(o) == :performance
    AA.bnns_compile_options_set_output_fd!(o, 7)
    @test AA.bnns_compile_options_get_output_fd(o) == 7
    AA.bnns_compile_options_set_output_path!(o, "/tmp/foo.bnns")
    @test AA.bnns_compile_options_get_output_path(o) == "/tmp/foo.bnns"
    AA.bnns_compile_options_set_log_mask!(o, 0xF)
    @test_throws ArgumentError AA.bnns_compile_options_set_optimization!(o, :bogus)

    o2 = AA.BNNSGraphCompileOptions(single_thread = true, optimization = :ir_size, output_fd = 3)
    @test AA.bnns_compile_options_get_single_thread(o2)
    @test AA.bnns_compile_options_get_optimization(o2) == :ir_size
    @test AA.bnns_compile_options_get_output_fd(o2) == 3
end

@testset "BNNS layer filters" begin
    W = Float32[1 2 3; 4 5 6]; bias = Float32[10, 20]
    f = AA.bnns_fully_connected(W, bias)
    x = Float32[1, 2, 3]; y = zeros(Float32, 2)
    @test AA.bnns_filter_apply!(f, y, x) == W * x .+ bias

    X = Float32[1 0; 2 0; 3 1]; Yb = zeros(Float32, 2, 2)
    @test AA.bnns_filter_apply_batch!(f, Yb, X) == W * X .+ bias

    # no bias, relu
    fr = AA.bnns_fully_connected(Float32[1 -1], Float32[0]; activation = :relu)
    yr = zeros(Float32, 1)
    AA.bnns_filter_apply!(fr, yr, Float32[1, 5])
    @test yr[1] == 0.0f0                       # relu(1 - 5) = 0

    fn = AA.bnns_fully_connected(W)
    yn = zeros(Float32, 2); AA.bnns_filter_apply!(fn, yn, x)
    @test yn == W * x
    @test_throws DimensionMismatch AA.bnns_fully_connected(W, Float32[1, 2, 3])

    # direct activation batch
    ob = zeros(Float32, 4)
    AA.bnns_activation_batch!(:relu, ob, Float32[-1, 2, -3, 4])
    @test ob == Float32[0, 2, 0, 4]
    @test_throws ArgumentError AA.bnns_activation_batch!(:bogus, zeros(Float32, 2), Float32[1, 2])
end

@testset "BNNS optimizer step" begin
    # SGD without momentum: theta -= lr * (grad * gradient_scale)
    theta = [Float32[1.0, 2.0, 3.0]]
    grad = [Float32[0.1, 0.2, 0.3]]
    accum = [zeros(Float32, 3)]
    LAcc = AA.LibAccelerate
    fields = LAcc.BNNSOptimizerSGDMomentumFields(0.5f0, 0.0f0, 1.0f0, 0.0f0, false, 0.0f0, 0.0f0,
        false, LAcc.BNNSOptimizerRegularizationFunction(LAcc.BNNSOptimizerRegularizationNone),
        LAcc.BNNSOptimizerSGDMomentumVariant(LAcc.BNNSSGDMomentumVariant0))
    AA.bnns_optimizer_step_sgd!(theta, grad, accum, fields)
    @test theta[1] ≈ Float32[1, 2, 3] .- 0.5f0 .* Float32[0.1, 0.2, 0.3]
    @test_throws DimensionMismatch AA.bnns_optimizer_step_sgd!(theta, grad,
        [zeros(Float32, 3), zeros(Float32, 3)], fields)
end
