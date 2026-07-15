# AppleAccelerate.jl Benchmark Suite
#
# Usage:
#   julia --project=test/bench test/bench/run_benchmarks.jl           # run all
#   julia --project=test/bench test/bench/run_benchmarks.jl fft        # run only FFT
#   julia --project=test/bench test/bench/run_benchmarks.jl sparse      # run only sparse
#   julia --project=test/bench test/bench/run_benchmarks.jl dense       # run only dense
#   julia --project=test/bench test/bench/run_benchmarks.jl array       # run only array ops
#
# Each bench file is launched in its OWN fresh Julia process. This is
# essential for a fair comparison: bench_dense.jl and bench_sparse.jl
# benchmark OpenBLAS/SuiteSparse *before* loading AppleAccelerate, then load
# it internally. If they shared one process, whichever ran second would find
# AppleAccelerate already loaded — its OpenBLAS/SuiteSparse "baseline" would
# then route BLAS through Accelerate (via libblastrampoline) and be
# contaminated. Separate processes guarantee each gets a clean pre-Accelerate
# baseline, and let the documented command reproduce the published numbers.
#
# All benchmarks run single-threaded (BLAS, FFTW, Accelerate) for fair,
# reproducible comparisons. Each bench file sets its own thread counts; set
# the BENCH_THREADS env var to sweep BLAS thread counts in the dense suite.

using Dates

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

# (key, banner, filename) in a stable, sensible display order.
const BENCHES = [
    ("dense",  "Dense LA Benchmarks",        "bench_dense.jl"),
    ("sparse", "Sparse Benchmarks",          "bench_sparse.jl"),
    ("fft",    "FFT Benchmarks",             "bench_fft.jl"),
    ("array",  "Array Operation Benchmarks", "bench_array.jl"),
]

# Reuse the current project (the bench environment) for each child process.
const PROJECT = something(Base.active_project(), joinpath(BENCH_DIR, "Project.toml"))

for (key, banner, file) in BENCHES
    key in requested || continue
    println("\n>>> Running $banner (isolated process) <<<\n")
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT) $(joinpath(BENCH_DIR, file))`
    run(pipeline(setenv(cmd, ENV); stdout, stderr))
end

println("\n" * "="^70)
println("All requested benchmarks complete.")
println("="^70)
