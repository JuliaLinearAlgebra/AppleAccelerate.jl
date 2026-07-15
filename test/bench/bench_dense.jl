# Dense Linear Algebra Benchmarks: OpenBLAS vs Apple Accelerate
#
# This script benchmarks OpenBLAS FIRST (before loading AppleAccelerate),
# then loads AppleAccelerate and re-runs the same benchmarks.
# Benchmarks run single-threaded by default (set BENCH_THREADS to sweep thread
# counts) for a fair, reproducible comparison. See issue #132.
#
# Benchmark count: ~42 problems per backend (×2 = 84 timing calls)
#   GEMM:            2 types × 6 sizes = 12
#   Factorizations:  4 ops × 5 sizes   = 20
#   Linear solve:    2 types × 5 sizes = 10

using LinearAlgebra, Printf

# Thread count for the BLAS backends. Defaults to 1 for a reproducible
# kernel-quality comparison, but can be overridden to measure how OpenBLAS
# scales across cores vs. Accelerate's (flat) SME co-processor throughput:
#   BENCH_THREADS=4 julia --project=test/bench test/bench/bench_dense.jl
# See issue #132.
const BENCH_THREADS = parse(Int, get(ENV, "BENCH_THREADS", "1"))
BLAS.set_num_threads(BENCH_THREADS)
@info "Dense benchmark BLAS threads" BENCH_THREADS

# Timing helper: 1 warmup call + N timed runs, return minimum
function bench_min(f; runs=5)
    f()  # warmup
    t_best = Inf
    for _ in 1:runs
        t_best = min(t_best, @elapsed f())
    end
    t_best
end

function gemm_gflops(N, elapsed_s)
    2.0 * N^3 / elapsed_s / 1e9
end

function bench_gemm(T, sizes)
    results = []
    for N in sizes
        A = randn(T, N, N)
        B = randn(T, N, N)
        C = similar(A)
        t = bench_min(() -> mul!(C, A, B))
        gf = gemm_gflops(N, t)
        push!(results, (N=N, time=t, gflops=gf))
        @printf("  %s N=%d: %.1f GFLOPS (%.1fμs)\n", T, N, gf, t*1e6)
    end
    results
end

function bench_factorizations(sizes)
    results = Dict{String, Vector}()
    for (name, make_matrix, f) in [
            ("LU",       N -> randn(Float64, N, N),                                    lu),
            ("QR",       N -> randn(Float64, N, N),                                    qr),
            ("Cholesky", N -> let M = randn(Float64, N, N); M'*M + N*I end,           cholesky),
            ("SVD",      N -> randn(Float64, N, N),                                    svd)]
        timings = []
        for N in sizes
            A = make_matrix(N)
            t = bench_min(() -> f(A))
            push!(timings, (N=N, time=t))
            @printf("  %s N=%d: %.1fμs\n", name, N, t*1e6)
        end
        results[name] = timings
    end
    results
end

function bench_solve(T, sizes)
    results = []
    for N in sizes
        A = randn(T, N, N) + T(N) * I
        b = randn(T, N)
        t = bench_min(() -> A \ b)
        push!(results, (N=N, time=t))
        @printf("  %s N=%d: %.1fμs\n", T, N, t*1e6)
    end
    results
end

gemm_sizes = (64, 256, 512, 1024, 2048, 4096)
fact_sizes = (64, 256, 512, 1024, 2048)
solve_sizes = (64, 256, 512, 1024, 2048)

# ============================================================
# Phase 1: OpenBLAS benchmarks (before loading AppleAccelerate)
# ============================================================
println("="^70)
println("Phase 1: OpenBLAS Benchmarks")
println("Current BLAS: ", BLAS.get_config())
println("="^70)

println("\nBenchmarking GEMM (OpenBLAS)...")
openblas_gemm_f64 = bench_gemm(Float64, gemm_sizes)
openblas_gemm_f32 = bench_gemm(Float32, gemm_sizes)

println("Benchmarking Factorizations (OpenBLAS)...")
openblas_facts = bench_factorizations(fact_sizes)

println("Benchmarking Linear Solve (OpenBLAS)...")
openblas_solve_f64 = bench_solve(Float64, solve_sizes)
openblas_solve_f32 = bench_solve(Float32, solve_sizes)

# ============================================================
# Phase 2: Load Accelerate and re-benchmark
# ============================================================
println("\n" * "="^70)
println("Phase 2: Loading AppleAccelerate...")
using AppleAccelerate
AppleAccelerate.set_num_threads(BENCH_THREADS)
BLAS.set_num_threads(BENCH_THREADS)
println("Current BLAS: ", BLAS.get_config())
println("="^70)

println("\nBenchmarking GEMM (Accelerate)...")
accel_gemm_f64 = bench_gemm(Float64, gemm_sizes)
accel_gemm_f32 = bench_gemm(Float32, gemm_sizes)

println("Benchmarking Factorizations (Accelerate)...")
accel_facts = bench_factorizations(fact_sizes)

println("Benchmarking Linear Solve (Accelerate)...")
accel_solve_f64 = bench_solve(Float64, solve_sizes)
accel_solve_f32 = bench_solve(Float32, solve_sizes)

# ============================================================
# Print comparison tables
# ============================================================
println("\n" * "="^70)
println("DENSE LINEAR ALGEBRA BENCHMARK RESULTS")
println("="^70)

# GEMM table
println("\n### GEMM (mul!) — GFLOPS (higher is better)\n")
println("| Type | N | OpenBLAS GFLOPS | Accelerate GFLOPS | Speedup |")
println("|------|---|-----------------|-------------------|---------|")
for (T_str, ob, ac) in [("Float64", openblas_gemm_f64, accel_gemm_f64),
                          ("Float32", openblas_gemm_f32, accel_gemm_f32)]
    for (o, a) in zip(ob, ac)
        ratio = a.gflops / o.gflops
        winner = ratio >= 1.0 ? "$(round(ratio, digits=2))×" : "$(round(1/ratio, digits=2))× slower"
        @printf("| %-7s | %4d | %15.1f | %17.1f | %s |\n",
                T_str, o.N, o.gflops, a.gflops, winner)
    end
end

# Factorizations table
println("\n### Factorizations (Float64) — time in μs (lower is better)\n")
println("| Operation | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |")
println("|-----------|---|---------------|-----------------|---------|")
for name in ("LU", "QR", "Cholesky", "SVD")
    ob = openblas_facts[name]
    ac = accel_facts[name]
    for (o, a) in zip(ob, ac)
        ratio = o.time / a.time
        winner = ratio >= 1.0 ? "$(round(ratio, digits=2))×" : "$(round(1/ratio, digits=2))× slower"
        @printf("| %-9s | %4d | %13.1f | %15.1f | %s |\n",
                name, o.N, o.time*1e6, a.time*1e6, winner)
    end
end

# Linear solve table
println("\n### Linear Solve (A\\b) — time in μs (lower is better)\n")
println("| Type | N | OpenBLAS (μs) | Accelerate (μs) | Speedup |")
println("|------|---|---------------|-----------------|---------|")
for (T_str, ob, ac) in [("Float64", openblas_solve_f64, accel_solve_f64),
                          ("Float32", openblas_solve_f32, accel_solve_f32)]
    for (o, a) in zip(ob, ac)
        ratio = o.time / a.time
        winner = ratio >= 1.0 ? "$(round(ratio, digits=2))×" : "$(round(1/ratio, digits=2))× slower"
        @printf("| %-7s | %4d | %13.1f | %15.1f | %s |\n",
                T_str, o.N, o.time*1e6, a.time*1e6, winner)
    end
end
println()
