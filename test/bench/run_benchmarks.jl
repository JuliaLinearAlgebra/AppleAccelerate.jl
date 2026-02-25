# AppleAccelerate.jl Benchmark Suite
#
# Usage:
#   julia --project test/bench/run_benchmarks.jl           # run all
#   julia --project test/bench/run_benchmarks.jl fft        # run only FFT
#   julia --project test/bench/run_benchmarks.jl sparse      # run only sparse
#   julia --project test/bench/run_benchmarks.jl dense       # run only dense
#   julia --project test/bench/run_benchmarks.jl array       # run only array ops
#
# NOTE: bench_dense.jl and bench_sparse.jl must be run BEFORE loading
# AppleAccelerate — they load it internally after the OpenBLAS/SuiteSparse
# benchmarks. The other bench files require AppleAccelerate to already be
# loaded. The runner handles this ordering.
#
# All benchmarks run single-threaded (BLAS, FFTW, Accelerate) for fair,
# reproducible comparisons. Each bench file sets its own thread counts.

using Dates, Printf

const BENCH_DIR = @__DIR__

println("="^70)
println("AppleAccelerate.jl Benchmark Suite")
println("="^70)
println("Julia: ", VERSION)
println("Threads: ", Threads.nthreads())
println("CPU: ", Sys.cpu_info()[1].model)
println("Date: ", Dates.now())
println("="^70)

const requested = if length(ARGS) > 0
    Set(ARGS)
else
    Set(["fft", "sparse", "dense", "array"])
end

# Dense and sparse MUST run first (before AppleAccelerate is loaded) because
# they benchmark OpenBLAS/SuiteSparse and then load AppleAccelerate internally.
if "dense" in requested
    println("\n>>> Running Dense LA Benchmarks <<<\n")
    include(joinpath(BENCH_DIR, "bench_dense.jl"))
end

if "sparse" in requested
    println("\n>>> Running Sparse Benchmarks <<<\n")
    include(joinpath(BENCH_DIR, "bench_sparse.jl"))
end

# Ensure AppleAccelerate is loaded for remaining benchmarks
if !isdefined(Main, :AppleAccelerate)
    using AppleAccelerate
end

if "fft" in requested
    println("\n>>> Running FFT Benchmarks <<<\n")
    include(joinpath(BENCH_DIR, "bench_fft.jl"))
end

if "array" in requested
    println("\n>>> Running Array Operation Benchmarks <<<\n")
    include(joinpath(BENCH_DIR, "bench_array.jl"))
end

println("\n" * "="^70)
println("All requested benchmarks complete.")
println("="^70)
