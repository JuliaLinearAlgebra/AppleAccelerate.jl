using AppleAccelerate
using AppleAccelerate: BNNSArray, bnns_data_type
using Statistics
using Test

const AA = AppleAccelerate

@testset "BNNS core" begin
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
    @test AA.bnns_transpose(M, 1, 2) == permutedims(M, (2, 1))
    C = Float32[1 2; 3 4; 5 6]
    @test AA.bnns_transpose(C, 1, 2) == permutedims(C, (2, 1))

    G = Float32[1 2 3; 4 5 6]
    @test AA.bnns_copy!(zeros(Float32, 2, 3), G) == G
end

@testset "BNNS reductions" begin
    G = Float32[1 2 3; 4 5 6]
    @test vec(AA.bnns_reduce(:sum, G; dim = 1)) == vec(sum(G, dims = 1))
    @test vec(AA.bnns_reduce(:sum, G; dim = 2)) == vec(sum(G, dims = 2))
    @test vec(AA.bnns_reduce(:max, G; dim = 2)) == vec(maximum(G, dims = 2))
    @test vec(AA.bnns_reduce(:min, G; dim = 1)) == vec(minimum(G, dims = 1))
    @test vec(AA.bnns_reduce(:mean, G; dim = 2)) ≈ vec(sum(G, dims = 2) ./ 3)
    @test_throws ArgumentError AA.bnns_reduce(:bogus, G)
end

@testset "BNNS utility queries" begin
    @test AA.bnns_layout_rank(AA.LibAccelerate.BNNSDataLayoutColumnMajorMatrix) == 2
    @test AA.bnns_layout_rank(AA.LibAccelerate.BNNSDataLayoutVector) == 1
    G = rand(Float32, 3, 4)
    @test AA.bnns_data_size(G) == sizeof(G)
    @test AA.bnns_tensor_allocation_size(G) >= sizeof(G)
end

@testset "BNNS DirectApply" begin
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
