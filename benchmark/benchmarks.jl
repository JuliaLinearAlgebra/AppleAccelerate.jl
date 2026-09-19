# Regression benchmark suite for AppleAccelerate.jl.
#
# Defines `SUITE::BenchmarkGroup`, the BenchmarkTools / PkgBenchmark convention. It reuses
# the workloads of test/bench/ (the scripts behind docs/src/benchmarks.md) but measures
# only the AppleAccelerate side: the question here is "did this revision get slower than
# the last release?", not "how does Accelerate compare with OpenBLAS/FFTW/SuiteSparse?".
# A subset of sizes keeps a full run to a few minutes.
#
#   julia --project=benchmark benchmark/compare.jl v0.8.0 HEAD    # compare two revisions
#   julia --project=benchmark -e 'include("benchmark/benchmarks.jl"); run(SUITE; verbose=true)'
#
# Everything runs single-threaded, like test/bench/, so results are comparable run to run.

using AppleAccelerate, BenchmarkTools, LinearAlgebra, Random, SparseArrays

AppleAccelerate.set_num_threads(1)
BLAS.set_num_threads(1)
Random.seed!(7)
# Per-benchmark time budget. ~90 benchmarks, so the 2 s default is a few minutes a revision.
BenchmarkTools.DEFAULT_PARAMETERS.seconds = parse(Float64, get(ENV, "BENCH_SECONDS", "2"))

const AA = AppleAccelerate
const SUITE = BenchmarkGroup()

# compare.jl runs *this* file against older revisions too; skip a benchmark whose API a
# revision predates instead of failing the whole run.
has(names::Symbol...) = all(n -> isdefined(AA, n), names)

# --- Array operations (test/bench/bench_array.jl) ---------------------------------------
let g = addgroup!(SUITE, "array")
    for T in (Float64, Float32), N in (10_000, 1_000_000)
        X = rand(T, N) .+ T(0.5)
        Y = rand(T, N) .+ T(0.5)
        for f in (:exp, :log, :sin, :cos, :sqrt, :abs, :sum, :maximum, :minimum)
            has(f) || continue
            fn = getfield(AA, f)
            g[string(f), string(T), N] = @benchmarkable $fn($X)
        end
        for f in (:vadd, :vmul, :dot)
            has(f) || continue
            fn = getfield(AA, f)
            g[string(f), string(T), N] = @benchmarkable $fn($X, $Y)
        end
    end
end

# --- Dense linear algebra through LBT (test/bench/bench_dense.jl) -----------------------
let g = addgroup!(SUITE, "dense")
    for T in (Float64, Float32), N in (256, 1024)
        A = randn(T, N, N); B = randn(T, N, N); C = similar(A); b = randn(T, N)
        S = Hermitian(A'A + N * I)
        g["gemm", string(T), N] = @benchmarkable mul!($C, $A, $B)
        g["lu", string(T), N] = @benchmarkable lu($A)
        g["cholesky", string(T), N] = @benchmarkable cholesky($S)
        g["solve", string(T), N] = @benchmarkable $A \ $b
    end
    let A = randn(512, 512)
        g["qr", "Float64", 512] = @benchmarkable qr($A)
        g["svd", "Float64", 512] = @benchmarkable svd($A)
    end
end

# --- FFT (test/bench/bench_fft.jl) ------------------------------------------------------
if has(:fft, :plan_fft)
    let g = addgroup!(SUITE, "fft")
        for T in (ComplexF64, ComplexF32), n in (1024, 65536)
            x = randn(T, n)
            setup = AA.plan_fft(x)
            g["planned 1D", string(T), n] = @benchmarkable AA.fft($x, $setup)
            g["no-plan 1D", string(T), n] = @benchmarkable AA.fft($x)   # setup-cache path
        end
        for T in (ComplexF64, ComplexF32)
            x = randn(T, 256, 256)
            setup = AA.plan_fft(x)
            g["planned 2D", string(T), 256] = @benchmarkable AA.fft($x, $setup)
        end
        if has(:rfft, :plan_rfft)
            for T in (Float64, Float32), n in (1024, 65536)
                x = randn(T, n)
                setup = AA.plan_rfft(x)
                g["planned real 1D", string(T), n] = @benchmarkable AA.rfft($x, $setup)
            end
        end
    end
end

# --- Sparse (test/bench/bench_sparse.jl) ------------------------------------------------
if has(:AASparseMatrix, :cholesky, :solve)
    let g = addgroup!(SUITE, "sparse")
        for T in (Float64, Float32), N in (1000, 5000)
            R = sprandn(T, N, N, 5 / N)
            S = SparseMatrixCSC{T,Int64}(R * R' + T(N) * I)
            A = AA.AASparseMatrix(S)
            x = randn(T, N)
            g["spmv", string(T), N] = @benchmarkable $A * $x
            g["cholesky+solve", string(T), N] =
                @benchmarkable AA.solve(AA.cholesky($A), $x)
        end
    end
end
