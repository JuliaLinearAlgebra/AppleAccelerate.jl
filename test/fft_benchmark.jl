using AppleAccelerate, FFTW

function run_benchmark(label)
    println("=== $label ===")
    println("  Accelerate threads: ", AppleAccelerate.get_num_threads())
    println("  FFTW threads: ", FFTW.nthreads())
    println()

    for T in (ComplexF64, ComplexF32)
        F = real(T)
        println("--- $T ---")
        for n in (64, 256, 1024, 4096, 16384, 65536, 262144, 1048576)
            r = randn(T, n)

            # Apple vDSP
            setup = AppleAccelerate.plan_fft(n, F)
            AppleAccelerate.fft(r, setup)  # warmup
            t_apple = minimum([minimum([@elapsed(AppleAccelerate.fft(r, setup)) for _ in 1:50]) for _ in 1:3])

            # FFTW
            plan_fftw = FFTW.plan_fft(r; flags=FFTW.MEASURE)
            plan_fftw * r  # warmup
            t_fftw = minimum([minimum([@elapsed(plan_fftw * r) for _ in 1:50]) for _ in 1:3])

            ratio = t_fftw / t_apple
            faster = ratio > 1 ? "vDSP $(round(ratio, digits=2))x faster" : "FFTW $(round(1/ratio, digits=2))x faster"
            println("  n=$(lpad(n, 7)):  vDSP=$(lpad(round(t_apple*1e6, digits=1), 8))us  FFTW=$(lpad(round(t_fftw*1e6, digits=1), 8))us  $faster")
        end
        println()
    end
end

# Both single-threaded (fairest comparison — vDSP FFT is always single-threaded)
AppleAccelerate.set_num_threads(1)
FFTW.set_num_threads(1)
run_benchmark("Both single-threaded")

# FFTW multi-threaded
ncores = AppleAccelerate.get_num_threads()  # restore default
AppleAccelerate.set_num_threads(ncores)
FFTW.set_num_threads(Sys.CPU_THREADS)
run_benchmark("vDSP (single-threaded) vs FFTW multi-threaded ($(Sys.CPU_THREADS) threads)")
