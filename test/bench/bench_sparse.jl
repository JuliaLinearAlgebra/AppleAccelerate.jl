# Sparse Solver Benchmarks: Apple Sparse Solvers vs SuiteSparse
#
# This script benchmarks SuiteSparse FIRST (before loading AppleAccelerate),
# then loads AppleAccelerate and re-runs the same benchmarks with Apple's
# sparse solvers. This ensures SuiteSparse uses OpenBLAS internally,
# not Accelerate-forwarded BLAS.
# All benchmarks run single-threaded for a fair, reproducible comparison.
#
# Benchmark count: ~28 problems
#   SpMV:     2 types × 4 sizes = 8
#   QR:       2 types × 5 sizes = 10
#   Cholesky: 2 types × 5 sizes = 10

using SparseArrays, LinearAlgebra, Printf

# Force single-threaded for fair comparison
BLAS.set_num_threads(1)

# Timing helper: 1 warmup call + N timed runs, return minimum
function bench_min(f; runs=5)
    f()  # warmup
    t_best = Inf
    for _ in 1:runs
        t_best = min(t_best, @elapsed f())
    end
    t_best
end

# Storage for results from both phases
const SUITE_SPMV = Dict{Tuple{Type,Int}, Float64}()
const SUITE_QR   = Dict{Tuple{Type,Int}, Float64}()
const SUITE_CHOL = Dict{Tuple{Type,Int}, Float64}()

# Pre-generate matrices so both phases use the same data
const spmv_sizes = (1000, 5000, 10000, 50000)
const solve_sizes = (100, 500, 1000, 2000, 5000)
const spmv_types = (Float64, Float32)

const SPMV_MATRICES = Dict{Tuple{Type,Int}, SparseMatrixCSC}()
const SPMV_VECTORS  = Dict{Tuple{Type,Int}, Vector}()
const QR_MATRICES   = Dict{Tuple{Type,Int}, SparseMatrixCSC}()
const QR_VECTORS    = Dict{Tuple{Type,Int}, Vector}()
const CHOL_MATRICES = Dict{Tuple{Type,Int}, SparseMatrixCSC}()
const CHOL_VECTORS  = Dict{Tuple{Type,Int}, Vector}()

println("Generating test matrices...")
for T in spmv_types
    for N in spmv_sizes
        SPMV_MATRICES[(T,N)] = sprandn(T, N, N, 0.01)
        SPMV_VECTORS[(T,N)]  = randn(T, N)
    end
    for N in solve_sizes
        QR_MATRICES[(T,N)]   = sprandn(T, N, N, 0.01) + T(20) * I
        QR_VECTORS[(T,N)]    = randn(T, N)
        A_raw = sprandn(T, N, N, 0.01)
        CHOL_MATRICES[(T,N)] = A_raw' * A_raw + T(N) * I
        CHOL_VECTORS[(T,N)]  = randn(T, N)
    end
end

# ============================================================
# Phase 1: SuiteSparse benchmarks (before loading AppleAccelerate)
# ============================================================
println("="^70)
println("Phase 1: SuiteSparse Benchmarks")
println("Current BLAS: ", BLAS.get_config())
println("="^70)

println("\nBenchmarking SpMV (SuiteSparse)...")
for T in spmv_types
    for N in spmv_sizes
        A_jl = SPMV_MATRICES[(T,N)]
        x = SPMV_VECTORS[(T,N)]
        t = bench_min(() -> A_jl * x)
        SUITE_SPMV[(T,N)] = t
        @printf("  %s N=%d: %.1fμs\n", T, N, t*1e6)
    end
end

println("Benchmarking QR Solve (SuiteSparse)...")
for T in spmv_types
    for N in solve_sizes
        A_jl = QR_MATRICES[(T,N)]
        b = QR_VECTORS[(T,N)]
        t = bench_min(() -> A_jl \ b)
        SUITE_QR[(T,N)] = t
        @printf("  %s N=%d: %.1fμs\n", T, N, t*1e6)
    end
end

println("Benchmarking Cholesky Solve (SuiteSparse)...")
for T in spmv_types
    for N in solve_sizes
        A_spd = CHOL_MATRICES[(T,N)]
        b = CHOL_VECTORS[(T,N)]
        t = bench_min(() -> begin
            F = cholesky(A_spd)
            F \ b
        end)
        SUITE_CHOL[(T,N)] = t
        @printf("  %s N=%d: %.1fμs\n", T, N, t*1e6)
    end
end

# ============================================================
# Phase 2: Load AppleAccelerate and benchmark Apple sparse solvers
# ============================================================
println("\n" * "="^70)
println("Phase 2: Loading AppleAccelerate...")
using AppleAccelerate
import AppleAccelerate: AASparseMatrix, AAFactorization, solve, factor!,
    SparseFactorizationQR, SparseFactorizationCholesky
AppleAccelerate.set_num_threads(1)
BLAS.set_num_threads(1)
println("Current BLAS: ", BLAS.get_config())
println("="^70)

const SPARSE_RESULTS = Dict{String, Vector{NamedTuple{(:type, :size, :apple_us, :suite_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}}}()

function print_sparse_table(label)
    results = SPARSE_RESULTS[label]
    println("\n### $label\n")
    println("| Type | Size | Apple (μs) | SuiteSparse (μs) | Speedup |")
    println("|------|------|-----------|-------------------|---------|")
    for r in results
        winner = r.speedup >= 1.0 ? "$(round(r.speedup, digits=2))×" : "$(round(1/r.speedup, digits=2))× slower"
        @printf("| %-7s | %6s | %10.1f | %17.1f | %s |\n",
                r.type, r.size, r.apple_us, r.suite_us, winner)
    end
end

# --- Sparse Matrix-Vector Multiply ---
println("\nBenchmarking SpMV (Apple)...")
let results = NamedTuple{(:type, :size, :apple_us, :suite_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}[]
    for T in spmv_types
        for N in spmv_sizes
            A_jl = SPMV_MATRICES[(T,N)]
            A_aa = AASparseMatrix(A_jl)
            x = SPMV_VECTORS[(T,N)]
            t_apple = bench_min(() -> A_aa * x)
            t_suite = SUITE_SPMV[(T,N)]
            ratio = t_suite / t_apple
            push!(results, (type=string(T), size=string(N), apple_us=t_apple*1e6, suite_us=t_suite*1e6, speedup=ratio))
            @printf("  %s N=%d: Apple=%.1fμs SuiteSparse=%.1fμs (%.2f×)\n", T, N, t_apple*1e6, t_suite*1e6, ratio)
        end
    end
    SPARSE_RESULTS["Sparse Matrix-Vector Multiply (density=0.01)"] = results
end

# --- QR Factorize + Solve ---
println("Benchmarking QR Factorize + Solve (Apple)...")
let results = NamedTuple{(:type, :size, :apple_us, :suite_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}[]
    for T in spmv_types
        for N in solve_sizes
            A_jl = QR_MATRICES[(T,N)]
            b = QR_VECTORS[(T,N)]
            t_apple = bench_min(() -> begin
                f = AAFactorization(A_jl)
                factor!(f, SparseFactorizationQR)
                solve(f, b)
            end)
            t_suite = SUITE_QR[(T,N)]
            ratio = t_suite / t_apple
            push!(results, (type=string(T), size=string(N), apple_us=t_apple*1e6, suite_us=t_suite*1e6, speedup=ratio))
            @printf("  %s N=%d: Apple=%.1fμs SuiteSparse=%.1fμs (%.2f×)\n", T, N, t_apple*1e6, t_suite*1e6, ratio)
        end
    end
    SPARSE_RESULTS["QR Factorize + Solve"] = results
end

# --- Cholesky Factorize + Solve ---
println("Benchmarking Cholesky Factorize + Solve (Apple)...")
let results = NamedTuple{(:type, :size, :apple_us, :suite_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}[]
    for T in spmv_types
        for N in solve_sizes
            A_spd = CHOL_MATRICES[(T,N)]
            b = CHOL_VECTORS[(T,N)]
            t_apple = bench_min(() -> begin
                f = AAFactorization(A_spd)
                factor!(f, SparseFactorizationCholesky)
                solve(f, b)
            end)
            t_suite = SUITE_CHOL[(T,N)]
            ratio = t_suite / t_apple
            push!(results, (type=string(T), size=string(N), apple_us=t_apple*1e6, suite_us=t_suite*1e6, speedup=ratio))
            @printf("  %s N=%d: Apple=%.1fμs SuiteSparse=%.1fμs (%.2f×)\n", T, N, t_apple*1e6, t_suite*1e6, ratio)
        end
    end
    SPARSE_RESULTS["Cholesky Factorize + Solve"] = results
end

# --- Print all results ---
println("\n" * "="^70)
println("SPARSE BENCHMARK RESULTS")
println("="^70)
for label in ("Sparse Matrix-Vector Multiply (density=0.01)", "QR Factorize + Solve", "Cholesky Factorize + Solve")
    print_sparse_table(label)
end
println()
