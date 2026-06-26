# BLAS/LAPACK forwarding requires macOS >= 13.4. This file is the LAST one
# included by runtests.jl, so on older systems (or when the version can't be
# determined) we can cleanly `exit(0)` here without skipping any other test files
# — the non-forwarding suites (quadrature, bnns, nnlib) have already run.
let v = AppleAccelerate.get_macos_version()
if v === nothing || v < v"13.4"
    @info("AppleAccelerate.jl needs macOS >= 13.4 for BLAS forwarding. Not testing forwarding capabilities.")
    exit(0)
end

# Set up a debugging fallback function that prints out a stacktrace if the LinearAlgebra
# tests end up calling a function that we don't have forwarded.
accelerate_header = """
import LinearAlgebra
function debug_missing_function()
    println("Missing BLAS/LAPACK function!")
    display(stacktrace())
end
LinearAlgebra.BLAS.lbt_set_default_func(@cfunction(debug_missing_function, Cvoid, ()))
import AppleAccelerate
"""
eval(Meta.parseall(accelerate_header))

config = BLAS.get_config()
@info("Running with $(length(config.loaded_libs)) libraries loaded:")
display(config.loaded_libs)

@testset "Dense Linear Algebra" begin

@testset "Accelerate Forwarding Sanity Tests" begin
    ver = AppleAccelerate.get_macos_version()
    if ver >= v"13.4"
        @test LinearAlgebra.peakflops() > 0
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :lp64).libname, "Accelerate")
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :ilp64).libname, "Accelerate")

        # Accelerate has `_rook` symbols:
        @test endswith(BLAS.lbt_find_backing_library("dsytrf_rook_", :ilp64).libname, "Accelerate")
    end
end

@testset "CBLAS dot test" begin
    a = ComplexF64[
        1 + 1im,
        2 - 2im,
        3 + 3im
    ]
    @test BLAS.dotc(a, a) ≈ ComplexF64(28)
    @test BLAS.dotu(a, a) ≈ ComplexF64(12im)

    a = Float32[1, 2, 3]
    @test BLAS.dot(a, a) ≈ 14f0
end

@testset "BLAS threading tests" begin
    if AppleAccelerate.get_macos_version() >= v"15"
        @test AppleAccelerate.set_num_threads(1) == 1
        @test AppleAccelerate.get_num_threads() == 1
        @test AppleAccelerate.set_num_threads(4) > 1
        @test AppleAccelerate.get_num_threads() > 1
    else
        @test AppleAccelerate.get_num_threads() == 1
        @test AppleAccelerate.set_num_threads(1) == 1
    end
    # Regression: set_num_threads(n ≤ 0) used to fall through both branches,
    # leaving retval == -1 and tripping the BLASSetThreading @assert. It must
    # reject non-positive thread counts up front. Guarded like the cases above
    # so it only runs where the threading API is available (macOS ≥ 15).
    if AppleAccelerate.get_macos_version() >= v"15"
        @test_throws ArgumentError AppleAccelerate.set_num_threads(0)
        @test_throws ArgumentError AppleAccelerate.set_num_threads(-1)
    end
end

@testset "Forwarding spine coverage" begin
    # --- forward_accelerate argument validation -------------------------------
    # ILP64 without new_lapack is unsupported and must error up front.
    @test_throws ArgumentError AppleAccelerate.forward_accelerate(:ilp64; new_lapack = false)

    # --- BLAS.get_config reflects Accelerate forwarding -----------------------
    cfg = BLAS.get_config()
    acc_libs = filter(l -> endswith(l.libname, "Accelerate"), cfg.loaded_libs)
    # Accelerate must be forwarded in BOTH interfaces after load_accelerate().
    @test any(l -> l.interface === :lp64, acc_libs)
    @test any(l -> l.interface === :ilp64, acc_libs)
    @test endswith(BLAS.lbt_find_backing_library("dgemm_", :lp64).libname, "Accelerate")
    @test endswith(BLAS.lbt_find_backing_library("dgemm_", :ilp64).libname, "Accelerate")

    # --- forwarding is correct, not merely present ----------------------------
    # Reconcile Accelerate's results against independent references.
    let n = 32
        for T in (Float64, Float32, ComplexF64)
            A = rand(T, n, n) + n * I        # well-conditioned
            B = rand(T, n, n)
            x = rand(T, n)

            # gemm: compare against an explicit BLAS.gemm! on the same backend
            # and a hand-rolled triple loop (independent of the BLAS path).
            C = A * B
            Cref = fill(zero(T), n, n)
            BLAS.gemm!('N', 'N', one(T), A, B, zero(T), Cref)
            @test C ≈ Cref
            Chand = [sum(A[i, k] * B[k, j] for k in 1:n) for i in 1:n, j in 1:n]
            @test C ≈ Chand

            # lu / solve: residual must be tiny.
            F = lu(A)
            @test F.L * F.U ≈ A[F.p, :]
            xs = A \ x
            @test A * xs ≈ x
            rtol = T <: Union{Float32,ComplexF32} ? 1e-3 : 1e-9
            @test norm(A * xs - x) / norm(x) < rtol
        end
    end

    # --- get_macos_version returns a consistent VersionNumber (caching) -------
    v1 = AppleAccelerate.get_macos_version()
    v2 = AppleAccelerate.get_macos_version()
    @test v1 isa VersionNumber
    @test v1 == v2                     # cached: identical across calls
    @test AppleAccelerate._read_macos_version() isa VersionNumber

    # Exercise the precompile-time fallback branch in get_macos_version where the
    # cache Ref is still `nothing` and it must read the plist directly.
    saved_ver = AppleAccelerate._macos_version[]
    try
        AppleAccelerate._macos_version[] = nothing
        vfallback = AppleAccelerate.get_macos_version()   # hits the _read_macos_version() fallback
        @test vfallback isa VersionNumber
        @test vfallback == saved_ver
    finally
        AppleAccelerate._macos_version[] = saved_ver
    end

    # --- threading API on the synthetic "too old" path ------------------------
    # On macOS < 15 the threading toggles must warn and return 1 without touching
    # the (unavailable) Accelerate threading symbols. We can't downgrade the OS,
    # but we can spoof the cached version to drive those guarded branches.
    let saved = AppleAccelerate._macos_version[]
        try
            AppleAccelerate._macos_version[] = v"14.0.0"
            @test AppleAccelerate.get_num_threads() == 1
            @test AppleAccelerate.set_num_threads(1) == 1
            @test AppleAccelerate.set_num_threads(8) == 1
            # Argument validation still fires before the version guard.
            @test_throws ArgumentError AppleAccelerate.set_num_threads(0)
        finally
            AppleAccelerate._macos_version[] = saved
        end
    end
end

linalg_stdlib_test_path = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")

@testset verbose=false "LinearAlgebra.jl BLAS tests" begin
    blas_test_path = joinpath(linalg_stdlib_test_path, "blas.jl")
    if VERSION < v"1.12"
        # On Julia < 1.12, the stdlib `blas.jl` strided-interface tests assert
        # that `herk!`/`her2k!` leave the unreferenced (strictly-lower /
        # diagonal-imaginary) part of `C` untouched. Apple's Accelerate writes
        # into that part instead, so those two assertions fail even though the
        # hermitian result itself is correct. Skip just those two assertions
        # (matched by content, since their line numbers move between versions)
        # while still running the rest of the stdlib BLAS suite.
        src = read(blas_test_path, String)
        src = replace(src,
            "@test C == WrappedArray([23 50+38im; 35+27im 152])" =>
                "@test_skip C == WrappedArray([23 50+38im; 35+27im 152])",
            "@test C == WrappedArray([63 138+38im; 35+27im 352])" =>
                "@test_skip C == WrappedArray([63 138+38im; 35+27im 352])")
        include_string(@__MODULE__, src, blas_test_path)
    else
        include(blas_test_path)
    end
end

@testset verbose=false "LinearAlgebra.jl LAPACK tests" begin
    joinpath(linalg_stdlib_test_path, "lapack.jl") |> include
end

end # @testset "Dense Linear Algebra"

end # let v
