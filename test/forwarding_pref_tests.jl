# Automatic BLAS/LAPACK forwarding on load can be switched off (issue #178). The switch is
# read in `__init__`, so every case runs in a fresh subprocess. Forwarding needs
# macOS >= 13.4; on older systems nothing is forwarded either way, so skip.
let v = AppleAccelerate.get_macos_version()
if v !== nothing && v >= v"13.4"

@testset "Automatic forwarding switch" begin
    # Prints "<forwarded after using> <forwarded after load_accelerate()>".
    probe = """
    using LinearAlgebra, AppleAccelerate
    forwarded() = any(l -> endswith(l.libname, "Accelerate"), BLAS.get_config().loaded_libs)
    before = forwarded()
    AppleAccelerate.load_accelerate()
    println(before, " ", forwarded())
    """
    function run_probe(code, env...)
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`
        return split(readchomp(addenv(cmd, env...)))
    end
    envname = AppleAccelerate._AUTO_FORWARD_ENV

    @testset "environment variable" begin
        # Default: forwards on load
        @test run_probe(probe, envname => nothing) == ["true", "true"]
        for off in ("0", "false", "OFF")
            # Disabled: BLAS untouched on load, load_accelerate() opts in
            @test run_probe(probe, envname => off) == ["false", "true"]
        end
        @test run_probe(probe, envname => "1") == ["true", "true"]
    end

    @testset "preference" begin
        set_pref(flag) = "using AppleAccelerate; AppleAccelerate.set_auto_forward!($flag)"
        try
            run_probe(set_pref(false), envname => nothing)
            @test run_probe(probe, envname => nothing) == ["false", "true"]
            # The environment variable takes precedence over the preference
            @test run_probe(probe, envname => "1") == ["true", "true"]
        finally
            run_probe(set_pref(true), envname => nothing)
        end
        @test run_probe(probe, envname => nothing) == ["true", "true"]
    end

    @testset "auto_forward_blas" begin
        @test withenv(() -> AppleAccelerate.auto_forward_blas(), envname => "no") == false
        @test withenv(() -> AppleAccelerate.auto_forward_blas(), envname => "yes") == true
        @test withenv(envname => "maybe") do
            @test_logs (:warn, r"unrecognized") AppleAccelerate.auto_forward_blas()
        end == true
    end
end

end
end
