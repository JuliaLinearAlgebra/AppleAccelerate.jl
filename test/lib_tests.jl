# Lock-in tests for the generated raw ABI layer (src/lib/LibAccelerate.jl).
#
# Clang.jl mis-resolves Apple's anonymous `_Complex` typedefs and emitted the
# double-precision complex element types as ComplexF32. The generator now
# post-processes these to ComplexF64; these assertions guard against regression
# (a wrong width would silently hand half-size buffers to double-complex symbols).
@testset "LibAccelerate complex typedefs" begin
    @test AppleAccelerate.LibAccelerate.__double_complex_t === ComplexF64
    @test AppleAccelerate.LibAccelerate.__SPARSE_double_complex === ComplexF64
    # Single-precision typedefs must stay ComplexF32.
    @test AppleAccelerate.LibAccelerate.__float_complex_t === ComplexF32
    @test AppleAccelerate.LibAccelerate.__SPARSE_float_complex === ComplexF32
end

# Drift guard for dead wrappers. Clang.jl turns every C declaration into a `function`
# wrapper, including inline-only / macro / non-exported "functions" (e.g. the high-level
# Sparse/Dense Solve inline API, header-only BNNS graph setters, the `CF_ENUM` macro).
# Those `@ccall libacc.<sym>(...)` blocks throw `could not load symbol` if ever called.
# The generator now strips any wrapper whose ccall symbol fails a `dlsym` probe; this test
# asserts the committed output is clean, so a future regen that re-introduces a dead wrapper
# fails CI rather than shipping a broken binding. We parse the generated file for every
# `@ccall libacc.<sym>(` and probe each symbol against the live framework.
@testset "LibAccelerate wrappers resolve to real symbols" begin
    import Libdl
    libpath = AppleAccelerate.LibAccelerate.libacc
    h = Libdl.dlopen(libpath)
    libsrc = joinpath(dirname(pathof(AppleAccelerate)), "lib", "LibAccelerate.jl")
    syms = Set{String}()
    for line in eachline(libsrc)
        m = match(r"@ccall\s+libacc\.([A-Za-z_][A-Za-z0-9_]*)\(", line)
        m === nothing || push!(syms, m.captures[1])
    end
    @test length(syms) > 800   # sanity: we actually scanned the wrapper layer
    dead = sort!([s for s in syms if Libdl.dlsym_e(h, s) == C_NULL])
    @test isempty(dead)
    isempty(dead) || @info "Generated wrappers ccall non-existent symbols" dead
end
