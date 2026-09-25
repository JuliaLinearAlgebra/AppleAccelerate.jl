using AppleAccelerate
using AppleAccelerate: BNNSArray, bnns_data_type
using Random: shuffle
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

        # Float16 descriptors carry the half-precision data type code.
        dh = BNNSArray(rand(Float16, 3, 2))
        @test dh.desc.data_type == bnns_data_type(Float16) != bnns_data_type(Float32)

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

    # Every describable element type, on a non-square 3-D array so an axis mix-up fails.
    for T in (Float16, Float32, Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64, Bool)
        A = T === Bool ? rand(Bool, 2, 3, 4) : T.(rand(1:100, 2, 3, 4))
        @test AA.bnns_transpose(A, 1, 3) == permutedims(A, (3, 2, 1))
        @test AA.bnns_transpose(A, 2, 3) == permutedims(A, (1, 3, 2))
        @test AA.bnns_copy!(similar(A), A) == A
    end
    @test_throws ArgumentError AA.bnns_transpose(G, 1, 3)
    @test_throws ArgumentError AA.bnns_transpose(rand(Float64, 2, 2), 1, 2)

    # Converting copies: the Float16 <-> Float32 bridge, ints -> Float32, Float32 -> Int32.
    F = Float32[1.5 -2.25 3; 4 5.125 -6]
    @test AA.bnns_copy!(zeros(Float16, 2, 3), F) == Float16.(F)
    @test AA.bnns_copy!(zeros(Float32, 2, 3), Float16.(F)) == F
    for T in (Int8, Int16, Int32, UInt8)
        @test AA.bnns_copy!(zeros(Float32, 2, 3), T.(G)) == G
    end
    @test AA.bnns_copy!(zeros(Int32, 2, 3), G) == Int32.(G)
    # Pairs BNNS does not implement must be rejected before the call: Float32 -> Int8
    # aborts the process inside BNNSCopy on Intel macOS 15 instead of returning a status.
    @test_throws ArgumentError AA.bnns_copy!(zeros(Int8, 2, 3), G)
    @test_throws ArgumentError AA.bnns_copy!(zeros(Float64, 2, 3), G)
    @test_throws ArgumentError AA.bnns_copy!(zeros(Float32, 2, 3), Int64.(G))
    # BNNSCopy does not broadcast; a shape mismatch must not reach it.
    @test_throws DimensionMismatch AA.bnns_copy!(zeros(Float32, 2, 3), Float32[1 2 3])
    @test_throws DimensionMismatch AA.bnns_copy!(zeros(Float32, 2, 3), zeros(Float32, 3, 2))
end

@testset "BNNS reductions" begin
    G = Float32[1 2 3; 4 5 6]
    @test vec(AA.bnns_reduce(:sum, G; dim = 1)) == vec(sum(G, dims = 1))
    @test vec(AA.bnns_reduce(:sum, G; dim = 2)) == vec(sum(G, dims = 2))
    @test vec(AA.bnns_reduce(:max, G; dim = 2)) == vec(maximum(G, dims = 2))
    @test vec(AA.bnns_reduce(:min, G; dim = 1)) == vec(minimum(G, dims = 1))
    @test vec(AA.bnns_reduce(:mean, G; dim = 2)) ≈ vec(sum(G, dims = 2) ./ 3)
    @test_throws ArgumentError AA.bnns_reduce(:bogus, G)
    @test_throws ArgumentError AA.bnns_reduce(:sum, G; dim = 3)

    refs = (sum = x -> sum(x, dims = 2), mean = x -> mean(x, dims = 2),
            max = x -> maximum(x, dims = 2), min = x -> minimum(x, dims = 2),
            sumsquare = x -> sum(abs2, x, dims = 2), l1 = x -> sum(abs, x, dims = 2),
            l2 = x -> sqrt.(sum(abs2, x, dims = 2)), product = x -> prod(x, dims = 2),
            logsumexp = x -> log.(sum(exp.(x), dims = 2)))
    @testset "reduce $T" for T in (Float16, Float32)
        A = T.(randn(3, 5, 2))
        for (f, ref) in pairs(refs)
            R = AA.bnns_reduce(f, A; dim = 2)
            @test R isa Array{T,3} && size(R) == (3, 1, 2)
            @test Float64.(R) ≈ ref(Float64.(A)) rtol = (T === Float16 ? 2e-2 : 1e-5)
        end
    end
    @testset "reduce Int32" begin
        A = Int32.(rand(-4:4, 3, 5, 2))
        for f in (:sum, :max, :min, :sumsquare, :l1, :product)
            R = AA.bnns_reduce(f, A; dim = 2)
            @test R isa Array{Int32,3} && R == refs[f](Int64.(A))
        end
        for f in (:mean, :l2, :logsumexp)     # inexact in integers
            @test_throws ArgumentError AA.bnns_reduce(f, A)
        end
    end
    # Other integer widths return status 0 with wrong values, so they are not offered.
    @test_throws MethodError AA.bnns_reduce(:sum, Int8[1 2; 3 4])
    @test_throws MethodError AA.bnns_reduce(:sum, rand(Float64, 2, 2))
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

    # top-K across element types, against a sort-based reference. Distinct values
    # per slice so the index comparison is unambiguous.
    @testset "topk $E" for E in (Float16, Float32, Int8, Int16, Int32, UInt8, UInt16)
        A = E.(reduce(hcat, [shuffle(1:9) for _ in 1:4]))          # 9 × 4
        for (dim, B) in ((1, A), (2, permutedims(A)))
            vals, inds = AA.bnns_topk(B, 3; dim)
            @test vals isa Matrix{E} && inds isa Matrix{Int32}
            @test vals == mapslices(c -> sort(c, rev = true)[1:3], B; dims = dim)
            @test inds == mapslices(c -> sortperm(c, rev = true)[1:3] .- 1, B; dims = dim)
        end
    end
    @test_throws MethodError AA.bnns_topk(Int64[1 2; 3 4], 1)   # BNNS rejects wide ints
    @test_throws ArgumentError AA.bnns_topk(T, 4; dim = 1)
    @test_throws ArgumentError AA.bnns_topk(T, 1; dim = 3)

    @testset "in_topk $E" for E in (Float16, Float32)
        s = E.(reshape(shuffle(1:200), 5, 40) ./ 200)
        t = Int32.(rand(0:4, 40))
        want = [t[j] in (sortperm(s[:, j], rev = true)[1:2] .- 1) for j in 1:40]
        @test AA.bnns_in_topk(s, t, 2; dim = 1) == want
    end
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

    @testset "uniform_int $E" for E in (Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64)
        Z = zeros(E, 5_000)
        AA.bnns_random_fill_uniform_int!(g, Z, 3, 9)
        @test sort(unique(Z)) == E.(3:8)             # half-open, every value hit
        @test_throws ArgumentError AA.bnns_random_fill_uniform_int!(g, Z, 9, 3)
    end
    @test extrema(AA.bnns_random_fill_uniform_int!(g, zeros(Int8, 5_000), -5, 5)) == (-5, 4)
    @test_throws ArgumentError AA.bnns_random_fill_uniform_int!(g, zeros(Int8, 4), 0, 300)
    @test_throws ArgumentError AA.bnns_random_fill_uniform_int!(g, zeros(UInt8, 4), -1, 3)

    # Half precision fills.
    H = zeros(Float16, 100_000)
    AA.bnns_random_fill_uniform!(g, H, 2, 4)
    @test all(2 .<= H .<= 4) && isapprox(mean(Float64.(H)), 3.0, atol = 0.02)
    AA.bnns_random_fill_normal!(g, H, 1, 2)
    @test isapprox(mean(Float64.(H)), 1.0, atol = 0.05) && isapprox(std(Float64.(H)), 2.0, atol = 0.05)
    Hc = zeros(Float16, 20_000)
    AA.bnns_random_fill_categorical!(g, Hc, Float16[0.1, 0.7, 0.2])
    @test [count(==(k), Hc) for k in 0:2] ./ 20_000 ≈ [0.1, 0.7, 0.2] atol = 0.02
    # mixed precisions mis-sample silently in BNNS, so they are not dispatchable
    @test_throws MethodError AA.bnns_random_fill_categorical!(g, Hc, Float32[0.5, 0.5])

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

# --- BNNS Graph: compile -> introspect -> run ---------------------------------
#
# `BNNSGraphCompileFromFile` takes a compiled Core ML model: a `.mlmodelc`
# directory whose `model.mil` is a *textual* MIL program. The programs below are
# written by hand (constants inline, so no weight blob is needed), which lets the
# whole pipeline be tested without Xcode's `coremlcompiler` or coremltools.
function _write_mlmodelc(body::AbstractString)
    mc = joinpath(mktempdir(), "model.mlmodelc")
    mkpath(mc)
    write(joinpath(mc, "model.mil"), """
    program(1.3)
    [buildInfo = dict<string, string>({{"coremlc-component-MIL", "handwritten"}})]
    {
    $body
    }
    """)
    return mc
end

# The BNNS Graph API needs macOS 15+.
if AA.get_macos_version() >= v"15"
@testset "BNNS graph inference" begin
    # z = relu(x * W + b): non-square everywhere, so a layout mix-up cannot pass.
    Wm = Float32[1 -2 3 0.5; -1 0.25 2 -3; 0.5 1 -1 2]          # 3 × 4
    bm = Float32[0.5, -1, 0.25, 2]
    # MIL tensor literals are nested, row by row.
    lit(v::AbstractVector) = "[" * join(v, ", ") * "]"
    lit(m::AbstractMatrix) = "[" * join((lit(r) for r in eachrow(m)), ", ") * "]"
    dense = _write_mlmodelc("""
        func main<ios16>(tensor<fp32, [2, 3]> x) {
                tensor<fp32, [3, 4]> W = const()[name = string("W"), val = tensor<fp32, [3, 4]>($(lit(Wm)))];
                tensor<fp32, [4]> b = const()[name = string("b"), val = tensor<fp32, [4]>($(lit(bm)))];
                bool no = const()[name = string("no"), val = bool(false)];
                tensor<fp32, [2, 4]> xw = matmul(transpose_x = no, transpose_y = no, x = x, y = W)[name = string("xw")];
                tensor<fp32, [2, 4]> s = add(x = xw, y = b)[name = string("s")];
                tensor<fp32, [2, 4]> z = relu(x = s)[name = string("z")];
            } -> (z);
    """)
    g = AA.BNNSGraph(dense)
    @test AA.bnns_graph_function_names(g) == ["main"]
    @test AA.bnns_graph_input_names(g) == ["x"]
    @test AA.bnns_graph_output_names(g) == ["z"]
    @test AA.bnns_graph_argument_names(g) == ["z", "x"]         # outputs first
    intents = AA.bnns_graph_argument_intents(g)
    @test intents isa Vector{Symbol}
    @test intents == [:out, :in]
    @test AA.bnns_graph_argument_position(g, "x") == 1
    @test (AA.bnns_graph_input_count(g), AA.bnns_graph_output_count(g), AA.bnns_graph_argument_count(g)) == (1, 1, 2)

    c = AA.BNNSGraphContext(g)
    info = AA.bnns_graph_arguments(c)
    @test [a.name for a in info] == ["z", "x"]
    @test info[1].intent === :out && info[2].intent === :in
    @test all(a -> a.eltype === Float32, info)
    @test info[1].shape == (2, 4) && info[1].size == (4, 2)     # Julia size = reversed MIL shape
    @test info[2].shape == (2, 3) && info[2].size == (3, 2)
    @test [a.name for a in AA.bnns_graph_arguments(g)] == ["z", "x"]
    @test occursin("shape=(2, 3)", sprint(show, info[2]))

    X = Float32[1 2 3; -4 5 -6]                                  # MIL-order input
    ref = max.(X * Wm .+ bm', 0)

    # zero-copy reversed-dimension convention
    out = AA.bnns_graph_run(c, "x" => permutedims(X))
    @test collect(keys(out)) == ["z"] && out["z"] isa Matrix{Float32}
    @test size(out["z"]) == (4, 2)
    @test permutedims(out["z"]) ≈ ref
    # model index order, inputs given as a NamedTuple / Dict
    @test AA.bnns_graph_run(c, (x = X,); mil_order = true)["z"] ≈ ref
    @test AA.bnns_graph_run(c, Dict("x" => X); mil_order = true)["z"] ≈ ref

    # in-place, reusing a page-aligned workspace
    ws = AA.bnns_graph_workspace(c)
    @test length(ws) >= AA.bnns_graph_context_workspace_size(c)
    @test UInt(pointer(ws)) % 16384 == 0
    Z = zeros(Float32, 4, 2)
    for _ in 1:3
        fill!(Z, NaN)
        @test AA.bnns_graph_run!(c, "z" => Z, "x" => permutedims(X); workspace = ws) == ("z" => Z)
        @test permutedims(Z) ≈ ref
    end
    # a second context on the same graph is independent
    @test AA.bnns_graph_run(AA.BNNSGraphContext(g), "x" => permutedims(X))["z"] == Z

    # argument validation happens before BNNS sees a pointer
    Xr = permutedims(X)
    @test_throws DimensionMismatch AA.bnns_graph_run(c, "x" => X)                       # un-reversed size
    @test_throws DimensionMismatch AA.bnns_graph_run!(c, "z" => zeros(Float32, 2, 4), "x" => Xr)
    @test_throws ArgumentError AA.bnns_graph_run(c, "x" => Float64.(Xr))                # eltype
    @test_throws ArgumentError AA.bnns_graph_run(c, "x" => view(Xr, :, :))              # not a dense Array
    @test_throws ArgumentError AA.bnns_graph_run(c, "nope" => Xr)
    @test_throws ArgumentError AA.bnns_graph_run(c)                                      # missing input
    @test_throws ArgumentError AA.bnns_graph_run!(c, "x" => Xr, "z" => Z)               # intents swapped
    @test_throws ArgumentError AA.bnns_graph_run!(c, ("z" => Z, "z" => Z), "x" => Xr)   # duplicate
    misaligned = unsafe_wrap(Vector{UInt8}, pointer(ws) + 1, length(ws) - 1)             # aliases `ws`
    @test_throws ArgumentError AA.bnns_graph_run!(c, "z" => Z, "x" => Xr; workspace = misaligned)
    @test_throws DimensionMismatch AA.bnns_graph_execute!(c, AA.LibAccelerate.bnns_graph_argument_t[])

    # compile failures throw instead of returning a null graph
    @test_throws ArgumentError AA.BNNSGraph(joinpath(mktempdir(), "missing.mlmodelc"))
    @test_throws ErrorException AA.BNNSGraph(_write_mlmodelc("    func main<ios16>(tensor<fp32, [2]> x) { nonsense } -> (x);"))

    @testset "Float16 graph" begin
        half = _write_mlmodelc("""
            func main<ios16>(tensor<fp16, [2, 3]> x) {
                    tensor<fp16, [3]> b = const()[name = string("b"), val = tensor<fp16, [3]>([1.0, -2.0, 0.5])];
                    tensor<fp16, [2, 3]> s = add(x = x, y = b)[name = string("s")];
                    tensor<fp16, [2, 3]> z = relu(x = s)[name = string("z")];
                } -> (z);
        """)
        ch = AA.BNNSGraphContext(AA.BNNSGraph(half))
        @test all(a -> a.eltype === Float16, AA.bnns_graph_arguments(ch))
        Xh = Float16[1 2 3; -4 0 6]
        @test AA.bnns_graph_run(ch, "x" => Xh; mil_order = true)["z"] == max.(Xh .+ Float16[1 -2 0.5], 0)
        # Float32 data in and out of a half-precision graph via the converting copy
        X32 = Float32[0.5 1.5 -2; 3 -0.25 8]
        zh = AA.bnns_graph_run(ch, "x" => AA.bnns_copy!(zeros(Float16, 2, 3), X32); mil_order = true)["z"]
        @test AA.bnns_copy!(zeros(Float32, 2, 3), zh) == max.(X32 .+ Float32[1 -2 0.5], 0)
        @test_throws ArgumentError AA.bnns_graph_run(ch, "x" => X32; mil_order = true)
    end

    @testset "Int32 / Bool graph, several outputs" begin
        ints = _write_mlmodelc("""
            func main<ios16>(tensor<int32, [4]> a, tensor<int32, [4]> b) {
                    tensor<int32, [4]> s = add(x = a, y = b)[name = string("s")];
                    tensor<bool, [4]> gt = greater(x = a, y = b)[name = string("gt")];
                } -> (s, gt);
        """)
        ci = AA.BNNSGraphContext(AA.BNNSGraph(ints))
        ai = AA.bnns_graph_arguments(ci)
        @test [(a.name, a.intent, a.eltype) for a in ai] ==
              [("s", :out, Int32), ("gt", :out, Bool), ("a", :in, Int32), ("b", :in, Int32)]
        a = Int32[1, 5, 3, 0]; b = Int32[2, 2, 3, -1]
        o = AA.bnns_graph_run(ci, "a" => a, "b" => b)
        @test o["s"] == a .+ b && o["s"] isa Vector{Int32}
        @test o["gt"] == (a .> b) && o["gt"] isa Vector{Bool}
    end

    @testset "dynamic batch / shapes" begin
        dyn = _write_mlmodelc("""
            func main<ios16>(tensor<fp32, [?, 3]> x) {
                    tensor<fp32, [?, 3]> z = relu(x = x)[name = string("z")];
                } -> (z);
        """)
        cd_ = AA.BNNSGraphContext(AA.BNNSGraph(dyn))
        @test AA.bnns_graph_arguments(cd_)[2].shape == (0, 3)        # 0 = not bound yet
        @test_throws ErrorException AA.bnns_graph_run(cd_, "x" => zeros(Float32, 3, 5))
        AA.bnns_graph_context_set_batch_size!(cd_, 5)
        @test AA.bnns_graph_arguments(cd_)[1].size == (3, 5)
        Xd = randn(Float32, 3, 5)
        @test AA.bnns_graph_run(cd_, "x" => Xd)["z"] == max.(Xd, 0)
        @test AA.bnns_graph_context_set_dynamic_shapes!(cd_, "x" => (3, 7)) == ["z" => (3, 7), "x" => (3, 7)]
        Xd = randn(Float32, 3, 7)
        @test AA.bnns_graph_run(cd_, "x" => Xd; workspace = AA.bnns_graph_workspace(cd_))["z"] == max.(Xd, 0)
        @test_throws DimensionMismatch AA.bnns_graph_run(cd_, "x" => randn(Float32, 3, 5))
        @test_throws ArgumentError AA.bnns_graph_context_set_dynamic_shapes!(cd_, "z" => (3, 7))
        @test_throws ArgumentError AA.bnns_graph_context_set_dynamic_shapes!(cd_, "q" => (3, 7))
        @test_throws ArgumentError AA.bnns_graph_context_set_batch_size!(cd_, 0)
    end

    # Graph memory is released only after its contexts: churn finalizers in both orders.
    for _ in 1:20
        gg = AA.BNNSGraph(dense); cc = AA.BNNSGraphContext(gg)
        AA.bnns_graph_run(cc, "x" => permutedims(X))
    end
    GC.gc(); GC.gc()
    @test AA.bnns_graph_run(c, "x" => permutedims(X))["z"] == Z
end
end
