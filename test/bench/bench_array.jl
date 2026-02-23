# Array Operation Benchmarks: Apple vDSP vs Julia Base
#
# Compares AppleAccelerate's vDSP-backed element-wise operations against
# Julia Base (manual loops to bypass the broadcast override).
# All benchmarks run single-threaded for a fair, reproducible comparison.
#
# Benchmark count: ~84 problems
#   Unary math:     6 ops × 2 types × 3 sizes = 36
#   Reductions:     4 ops × 2 types × 3 sizes = 24
#   Binary ops:     2 ops × 2 types × 3 sizes = 12
#   Compound ops:   2 ops × 2 types × 3 sizes = 12

using AppleAccelerate, Printf, LinearAlgebra

# Force single-threaded for fair comparison
AppleAccelerate.set_num_threads(1)
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

# Julia reference implementations (bypass broadcast overrides)
julia_exp(X) = map(Base.exp, X)
julia_log(X) = map(Base.log, X)
julia_sin(X) = map(Base.sin, X)
julia_cos(X) = map(Base.cos, X)
julia_sqrt(X) = map(Base.sqrt, X)
julia_abs(X) = map(Base.abs, X)

julia_sum(X) = Base.sum(X)
julia_maximum(X) = Base.maximum(X)
julia_minimum(X) = Base.minimum(X)
julia_dot(X, Y) = LinearAlgebra.dot(X, Y)

function julia_vadd(X, Y)
    out = similar(X)
    @inbounds @simd for i in eachindex(X)
        out[i] = X[i] + Y[i]
    end
    out
end

function julia_vmul(X, Y)
    out = similar(X)
    @inbounds @simd for i in eachindex(X)
        out[i] = X[i] * Y[i]
    end
    out
end

function julia_vam(A, B, C)
    out = similar(A)
    @inbounds @simd for i in eachindex(A)
        out[i] = (A[i] + B[i]) * C[i]
    end
    out
end

function julia_vsma(A, b, C)
    out = similar(A)
    @inbounds @simd for i in eachindex(A)
        out[i] = A[i] * b + C[i]
    end
    out
end

const ARRAY_RESULTS = Dict{String, Vector{NamedTuple{(:op, :type, :size, :vdsp_us, :julia_us, :speedup), Tuple{String, String, String, Float64, Float64, Float64}}}}()

function print_array_table(label)
    results = ARRAY_RESULTS[label]
    println("\n### $label\n")
    println("| Op | Type | Size | vDSP (μs) | Julia (μs) | Speedup |")
    println("|----|------|------|-----------|------------|---------|")
    for r in results
        winner = r.speedup >= 1.0 ? "$(round(r.speedup, digits=2))×" : "$(round(1/r.speedup, digits=2))× slower"
        @printf("| %-6s | %-7s | %8s | %9.1f | %10.1f | %s |\n",
                r.op, r.type, r.size, r.vdsp_us, r.julia_us, winner)
    end
end

array_sizes = (10_000, 100_000, 1_000_000)

# --- Unary Math Functions ---
println("Benchmarking Unary Math Functions...")
let results = NamedTuple{(:op, :type, :size, :vdsp_us, :julia_us, :speedup), Tuple{String, String, String, Float64, Float64, Float64}}[]
    for T in (Float64, Float32)
        for N in array_sizes
            X = abs.(randn(T, N)) .+ T(0.1)  # positive for log/sqrt
            for (name, aa_f, jl_f) in [
                    ("exp",  AppleAccelerate.exp,  julia_exp),
                    ("log",  AppleAccelerate.log,  julia_log),
                    ("sin",  AppleAccelerate.sin,  julia_sin),
                    ("cos",  AppleAccelerate.cos,  julia_cos),
                    ("sqrt", AppleAccelerate.sqrt, julia_sqrt),
                    ("abs",  AppleAccelerate.abs,  julia_abs)]
                t_vdsp = bench_min(() -> aa_f(X))
                t_julia = bench_min(() -> jl_f(X))
                ratio = t_julia / t_vdsp
                push!(results, (op=name, type=string(T), size=string(N), vdsp_us=t_vdsp*1e6, julia_us=t_julia*1e6, speedup=ratio))
                @printf("  %s %s N=%d: vDSP=%.1fμs Julia=%.1fμs (%.2f×)\n", name, T, N, t_vdsp*1e6, t_julia*1e6, ratio)
            end
        end
    end
    ARRAY_RESULTS["Unary Math Functions"] = results
end

# --- Reductions ---
println("Benchmarking Reductions...")
let results = NamedTuple{(:op, :type, :size, :vdsp_us, :julia_us, :speedup), Tuple{String, String, String, Float64, Float64, Float64}}[]
    for T in (Float64, Float32)
        for N in array_sizes
            X = randn(T, N)
            Y = randn(T, N)
            for (name, aa_f, jl_f) in [
                    ("sum",     () -> AppleAccelerate.sum(X),     () -> julia_sum(X)),
                    ("maximum", () -> AppleAccelerate.maximum(X), () -> julia_maximum(X)),
                    ("minimum", () -> AppleAccelerate.minimum(X), () -> julia_minimum(X)),
                    ("dot",     () -> AppleAccelerate.dot(X, Y),  () -> julia_dot(X, Y))]
                t_vdsp = bench_min(aa_f)
                t_julia = bench_min(jl_f)
                ratio = t_julia / t_vdsp
                push!(results, (op=name, type=string(T), size=string(N), vdsp_us=t_vdsp*1e6, julia_us=t_julia*1e6, speedup=ratio))
                @printf("  %s %s N=%d: vDSP=%.1fμs Julia=%.1fμs (%.2f×)\n", name, T, N, t_vdsp*1e6, t_julia*1e6, ratio)
            end
        end
    end
    ARRAY_RESULTS["Reductions"] = results
end

# --- Binary Element-wise Ops ---
println("Benchmarking Binary Element-wise Ops...")
let results = NamedTuple{(:op, :type, :size, :vdsp_us, :julia_us, :speedup), Tuple{String, String, String, Float64, Float64, Float64}}[]
    for T in (Float64, Float32)
        for N in array_sizes
            X = randn(T, N)
            Y = randn(T, N)
            for (name, aa_f, jl_f) in [
                    ("vadd", () -> AppleAccelerate.vadd(X, Y), () -> julia_vadd(X, Y)),
                    ("vmul", () -> AppleAccelerate.vmul(X, Y), () -> julia_vmul(X, Y))]
                t_vdsp = bench_min(aa_f)
                t_julia = bench_min(jl_f)
                ratio = t_julia / t_vdsp
                push!(results, (op=name, type=string(T), size=string(N), vdsp_us=t_vdsp*1e6, julia_us=t_julia*1e6, speedup=ratio))
                @printf("  %s %s N=%d: vDSP=%.1fμs Julia=%.1fμs (%.2f×)\n", name, T, N, t_vdsp*1e6, t_julia*1e6, ratio)
            end
        end
    end
    ARRAY_RESULTS["Binary Element-wise Ops"] = results
end

# --- Compound Ops ---
println("Benchmarking Compound Ops...")
let results = NamedTuple{(:op, :type, :size, :vdsp_us, :julia_us, :speedup), Tuple{String, String, String, Float64, Float64, Float64}}[]
    for T in (Float64, Float32)
        for N in array_sizes
            A = randn(T, N)
            B = randn(T, N)
            C = randn(T, N)
            b = T(2.5)
            for (name, aa_f, jl_f) in [
                    ("vam",  () -> AppleAccelerate.vam(A, B, C),   () -> julia_vam(A, B, C)),
                    ("vsma", () -> AppleAccelerate.vsma(A, b, C),  () -> julia_vsma(A, b, C))]
                t_vdsp = bench_min(aa_f)
                t_julia = bench_min(jl_f)
                ratio = t_julia / t_vdsp
                push!(results, (op=name, type=string(T), size=string(N), vdsp_us=t_vdsp*1e6, julia_us=t_julia*1e6, speedup=ratio))
                @printf("  %s %s N=%d: vDSP=%.1fμs Julia=%.1fμs (%.2f×)\n", name, T, N, t_vdsp*1e6, t_julia*1e6, ratio)
            end
        end
    end
    ARRAY_RESULTS["Compound Ops"] = results
end

# --- Print all results ---
println("\n" * "="^70)
println("ARRAY OPERATION BENCHMARK RESULTS")
println("="^70)
for label in ("Unary Math Functions", "Reductions", "Binary Element-wise Ops", "Compound Ops")
    print_array_table(label)
end
println()
