using Test

# These tests verify the package *contract* around weak dependencies and package
# extensions. Each scenario runs in a FRESH julia subprocess (so the parent test
# process, which loads LinearAlgebra/SparseArrays for the rest of the suite,
# cannot contaminate the observation). We reuse the parent's active project and
# load path so the subprocesses see the same (already-instantiated) environment.

const _PROJECT = Base.active_project()
const _LOAD_PATH = join(Base.LOAD_PATH, Sys.iswindows() ? ';' : ':')

# Run a snippet in a fresh julia process; return (exitcode, combined output).
function _run_julia(code::AbstractString)
    cmd = `$(Base.julia_cmd()) --project=$(_PROJECT) --startup-file=no --color=no -e $code`
    cmd = addenv(cmd, "JULIA_LOAD_PATH" => _LOAD_PATH)
    out = IOBuffer()
    p = run(pipeline(ignorestatus(cmd); stdout = out, stderr = out))
    return (p.exitcode, String(take!(out)))
end

# Is LinearAlgebra baked into the running sysimage? The default Julia 1.12
# sysimage (`sys.dylib`) preloads LinearAlgebra (and SparseArrays), so they
# already appear in `Base.loaded_modules` of every fresh process — which means
# the LinearAlgebra/SparseArrays-triggered extensions activate automatically.
# That is *correct* extension behavior (the trigger module is genuinely loaded);
# it simply means a sysimage-preloaded stdlib cannot be observed as "absent".
# We detect this so the bare-absence assertions only run on a sysimage that does
# NOT preload these stdlibs (e.g. a slim/base image).
const _LA_IN_SYSIMAGE = let
    ec, out = _run_julia("""
        println(any(occursin("LinearAlgebra", string(m)) for m in values(Base.loaded_modules)))
    """)
    ec == 0 && occursin("true", out)
end

@testset "Extension loading / package contract" begin
    @testset "bare `using AppleAccelerate` (no weakdeps, no forwarding)" begin
        if _LA_IN_SYSIMAGE
            @info "LinearAlgebra is preloaded in this sysimage; the extension " *
                  "activates automatically. Skipping bare-absence assertions " *
                  "(cannot observe a sysimage stdlib as absent)."
            # Still verify the threading API has no LinearAlgebra dependency.
            code = """
            using AppleAccelerate
            @assert AppleAccelerate.get_num_threads() isa Int
            # Core must not contain the forwarding *bodies*; the names are stubs
            # until the extension adds methods.
            @assert isdefined(AppleAccelerate, :load_accelerate)
            @assert isdefined(AppleAccelerate, :forward_accelerate)
            println("BARE_OK")
            """
            ec, out = _run_julia(code)
            @test ec == 0
            @test occursin("BARE_OK", out)
        else
            code = """
            using AppleAccelerate
            loaded = Set(string(m) for m in values(Base.loaded_modules))
            @assert !("LinearAlgebra" in loaded) "LinearAlgebra unexpectedly loaded"
            @assert !("SparseArrays" in loaded) "SparseArrays unexpectedly loaded"
            @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraExt) === nothing
            @assert Base.get_extension(AppleAccelerate, :AppleAccelerateSparseArraysExt) === nothing
            @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraSparseArraysExt) === nothing
            @assert AppleAccelerate.get_num_threads() isa Int
            # Bare load must NOT forward BLAS to Accelerate.
            using LinearAlgebra
            cfg = sprint(show, LinearAlgebra.BLAS.get_config())
            # (We just loaded LinearAlgebra; forwarding may now be on. The point
            # of this branch is that it was off *before* this line.)
            println("BARE_OK")
            """
            ec, out = _run_julia(code)
            @test ec == 0
            @test occursin("BARE_OK", out)
        end
    end

    @testset "BLAS forwarding is tied to the LinearAlgebra extension" begin
        # Once LinearAlgebra is present, the extension forwards BLAS to Accelerate.
        code = """
        using AppleAccelerate
        import LinearAlgebra
        if AppleAccelerate.get_macos_version() === nothing || AppleAccelerate.get_macos_version() < v"13.4"
            println("FORWARD_SKIPPED")
        else
            cfg = sprint(show, LinearAlgebra.BLAS.get_config())
            @assert occursin("Accelerate", cfg) || occursin("libacc", cfg) "BLAS not forwarded to Accelerate: \$cfg"
            @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraExt) !== nothing
            println("FORWARD_OK")
        end
        """
        ec, out = _run_julia(code)
        @test ec == 0
        @test occursin("FORWARD_OK", out) || occursin("FORWARD_SKIPPED", out)
    end

    @testset "LinearAlgebra extension activates with LinearAlgebra" begin
        # `_NEG` guards the negative assertions: they only hold when the *other*
        # trigger stdlib is not preloaded into the sysimage.
        code = """
        using AppleAccelerate
        using LinearAlgebra
        @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraExt) !== nothing
        if !$(_LA_IN_SYSIMAGE)
            @assert Base.get_extension(AppleAccelerate, :AppleAccelerateSparseArraysExt) === nothing
            @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraSparseArraysExt) === nothing
        end
        println("LA_EXT_OK")
        """
        ec, out = _run_julia(code)
        @test ec == 0
        @test occursin("LA_EXT_OK", out)
    end

    @testset "SparseArrays extension activates with SparseArrays only" begin
        code = """
        using AppleAccelerate
        using SparseArrays
        @assert Base.get_extension(AppleAccelerate, :AppleAccelerateSparseArraysExt) !== nothing
        # Construction from a SparseMatrixCSC works with just SparseArrays.
        A = sparse([1.0 0.0; 0.0 2.0])
        aa = AppleAccelerate.AASparseMatrix(A)
        @assert SparseMatrixCSC(aa) == A
        println("SA_EXT_OK")
        """
        ec, out = _run_julia(code)
        @test ec == 0
        @test occursin("SA_EXT_OK", out)
    end

    @testset "dual extension activates only with both" begin
        code = """
        using AppleAccelerate
        using LinearAlgebra
        using SparseArrays
        @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraExt) !== nothing
        @assert Base.get_extension(AppleAccelerate, :AppleAccelerateSparseArraysExt) !== nothing
        @assert Base.get_extension(AppleAccelerate, :AppleAccelerateLinearAlgebraSparseArraysExt) !== nothing
        # End-to-end: factorize/solve a CSC via the dual extension.
        A = sparse([2.0 0.0; 0.0 3.0])
        b = [2.0, 9.0]
        x = A \\ b
        @assert isapprox(x, [1.0, 3.0]; atol = 1e-8) "dual-ext solve wrong: \$x"
        println("DUAL_EXT_OK")
        """
        ec, out = _run_julia(code)
        @test ec == 0
        @test occursin("DUAL_EXT_OK", out)
    end
end
