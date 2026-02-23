# FFT/DSP Benchmarks: Apple vDSP vs FFTW
#
# Compares planned FFT performance for complex and real transforms.
# All benchmarks run single-threaded for a fair, reproducible comparison.
#
# Benchmark count: ~36 problems
#   Complex 1D: 2 types × 7 sizes = 14
#   Real FFT:   2 types × 6 sizes = 12
#   Complex 2D: 2 types × 5 sizes = 10

using AppleAccelerate, FFTW, Printf, LinearAlgebra

# Force single-threaded for fair comparison
AppleAccelerate.set_num_threads(1)
FFTW.set_num_threads(1)
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

const FFT_RESULTS = Dict{String, Vector{NamedTuple{(:type, :size, :vdsp_us, :fftw_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}}}()

function print_fft_table(label)
    results = FFT_RESULTS[label]
    println("\n### $label\n")
    println("| Type | Size | vDSP (μs) | FFTW (μs) | Speedup |")
    println("|------|------|-----------|-----------|---------|")
    for r in results
        winner = r.speedup >= 1.0 ? "$(round(r.speedup, digits=2))×" : "$(round(1/r.speedup, digits=2))× slower"
        @printf("| %-10s | %10s | %9.1f | %9.1f | %s |\n",
                r.type, r.size, r.vdsp_us, r.fftw_us, winner)
    end
end

# --- Complex 1D FFT ---
println("Benchmarking Complex 1D FFT...")
let results = NamedTuple{(:type, :size, :vdsp_us, :fftw_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}[]
    for T in (ComplexF64, ComplexF32)
        F = real(T)
        for n in (256, 1024, 4096, 16384, 65536, 262144, 1048576)
            x = randn(T, n)
            setup = AppleAccelerate.plan_fft(n, F)
            x_copy = copy(x)
            plan_fw = FFTW.plan_fft(x_copy; flags=FFTW.MEASURE)

            t_vdsp = bench_min(() -> AppleAccelerate.fft(x, setup))
            t_fftw = bench_min(() -> plan_fw * x_copy)
            ratio = t_fftw / t_vdsp
            push!(results, (type=string(T), size=string(n), vdsp_us=t_vdsp*1e6, fftw_us=t_fftw*1e6, speedup=ratio))
            @printf("  %s n=%d: vDSP=%.1fμs FFTW=%.1fμs (%.2f×)\n", T, n, t_vdsp*1e6, t_fftw*1e6, ratio)
        end
    end
    FFT_RESULTS["Complex 1D FFT"] = results
end

# --- Real FFT (pre-planned) ---
println("Benchmarking Real FFT...")
let results = NamedTuple{(:type, :size, :vdsp_us, :fftw_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}[]
    for T in (Float64, Float32)
        for n in (256, 1024, 4096, 16384, 65536, 262144)
            x = randn(T, n)
            setup = AppleAccelerate.plan_rfft(x)
            x_copy = copy(x)
            plan_fw = FFTW.plan_rfft(x_copy; flags=FFTW.MEASURE)

            t_vdsp = bench_min(() -> AppleAccelerate.rfft(x, setup))
            t_fftw = bench_min(() -> plan_fw * x_copy)
            ratio = t_fftw / t_vdsp
            push!(results, (type=string(T), size=string(n), vdsp_us=t_vdsp*1e6, fftw_us=t_fftw*1e6, speedup=ratio))
            @printf("  %s n=%d: vDSP=%.1fμs FFTW=%.1fμs (%.2f×)\n", T, n, t_vdsp*1e6, t_fftw*1e6, ratio)
        end
    end
    FFT_RESULTS["Real FFT"] = results
end

# --- Complex 2D FFT ---
println("Benchmarking Complex 2D FFT...")
let results = NamedTuple{(:type, :size, :vdsp_us, :fftw_us, :speedup), Tuple{String, String, Float64, Float64, Float64}}[]
    for T in (ComplexF64, ComplexF32)
        F = real(T)
        for n in (16, 32, 64, 128, 256)
            x = randn(T, n, n)
            setup = AppleAccelerate.plan_fft(x)
            x_copy = copy(x)
            plan_fw = FFTW.plan_fft(x_copy; flags=FFTW.MEASURE)

            t_vdsp = bench_min(() -> AppleAccelerate.fft(x, setup))
            t_fftw = bench_min(() -> plan_fw * x_copy)
            ratio = t_fftw / t_vdsp
            push!(results, (type=string(T), size="$(n)×$(n)", vdsp_us=t_vdsp*1e6, fftw_us=t_fftw*1e6, speedup=ratio))
            @printf("  %s %d×%d: vDSP=%.1fμs FFTW=%.1fμs (%.2f×)\n", T, n, n, t_vdsp*1e6, t_fftw*1e6, ratio)
        end
    end
    FFT_RESULTS["Complex 2D FFT"] = results
end

# --- Print all results ---
println("\n" * "="^70)
println("FFT BENCHMARK RESULTS")
println("="^70)
for label in ("Complex 1D FFT", "Real FFT", "Complex 2D FFT")
    print_fft_table(label)
end
println()
