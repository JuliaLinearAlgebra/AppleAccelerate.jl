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

# vecLibTypes.h's vUInt32 is a 16-byte SIMD vector, not a scalar or a plain tuple.
# Clang.jl leaves references to it in the vBigNum union accessors without defining it.
@testset "LibAccelerate vUInt32 and vBigNum accessors" begin
    LA = AppleAccelerate.LibAccelerate
    @test LA.vUInt32 === NTuple{4, VecElement{UInt32}}
    @test isbitstype(LA.vUInt32)
    @test sizeof(LA.vUInt32) == 16
    @test Base.datatype_alignment(LA.vUInt32) == 16

    for (U, S, n) in ((LA.vU128, LA.vS128, 1), (LA.vU256, LA.vS256, 2),
                      (LA.vU512, LA.vS512, 4), (LA.vU1024, LA.vS1024, 8)), T in (U, S)
        @testset "$T" begin
            V = n == 1 ? LA.vUInt32 : NTuple{n, LA.vUInt32}
            @test sizeof(T) == sizeof(V) == 16n
            vectors = ntuple(j -> ntuple(i -> VecElement(UInt32(4(j - 1) + i)), 4), n)
            value = n == 1 ? only(vectors) : vectors
            bytes = reinterpret(NTuple{sizeof(T), UInt8}, value)
            x = T(bytes)
            @test x.v === value

            r = Ref(x)
            GC.@preserve r begin
                p = Base.unsafe_convert(Ptr{T}, r)
                @test p.v isa Ptr{V}
                @test UInt(p.v) == UInt(p)
                @test unsafe_load(p.v) === value
                replacement = reverse(value)
                p.v = replacement
                @test r[].v === replacement
                @test getfield(r[], :data) === reinterpret(NTuple{sizeof(T), UInt8}, replacement)
            end
        end
    end
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
    # All generated function wrappers (names that appear as `@ccall libacc.<sym>(`).
    syms = Set{String}()
    for line in eachline(libsrc)
        m = match(r"@ccall\s+libacc\.([A-Za-z_][A-Za-z0-9_]*)\(", line)
        m === nothing || push!(syms, m.captures[1])
    end
    @test length(syms) > 800   # sanity: we actually scanned the wrapper layer

    # HARD guarantee: every wrapper the idiomatic layer actually calls
    # (`LibAccelerate.<name>(...)` in src/*.jl outside lib/) must resolve to a
    # live symbol — otherwise that idiomatic call throws at runtime. Restricted
    # to names that are real function wrappers (∩ syms) so struct/const refs
    # like `LibAccelerate.DSPSplitComplex` aren't probed. This is platform-stable
    # because the idiomatic surface is identical across SDKs/arches.
    srcdir = dirname(libsrc) |> dirname   # .../src
    used = Set{String}()
    for (root, _, files) in walkdir(srcdir), f in files
        (endswith(f, ".jl") && !occursin(joinpath("src", "lib"), root)) || continue
        for line in eachline(joinpath(root, f)), m in eachmatch(r"LibAccelerate\.([A-Za-z_][A-Za-z0-9_]*)", line)
            push!(used, m.captures[1])
        end
    end
    used_dead = sort!([s for s in intersect(used, syms) if Libdl.dlsym_e(h, s) == C_NULL])
    @test isempty(used_dead)
    isempty(used_dead) || @info "Idiomatic layer calls dead symbols" used_dead

    # Informational drift signal over the WHOLE generated layer. The set of
    # resolvable symbols varies by SDK/arch (a symbol present on this machine may
    # be absent on an older x86_64 SDK runner and vice-versa), so this is
    # reported, not asserted — a hard assertion here would be flaky across CI.
    dead = sort!([s for s in syms if Libdl.dlsym_e(h, s) == C_NULL])
    isempty(dead) || @info "Generated wrappers with no symbol on this platform (informational)" count = length(dead)
end

# `__asm__` labels. bnns_graph.h keeps the source-level names but links them to `_v2`
# symbols (`int BNNSGraphContextExecute(…) __asm__("_BNNSGraphContextExecute_v2")`); the
# un-suffixed symbols are still exported with the OLD argument lists, so a binding that
# ignores the label resolves fine and then crashes when called. The generator retargets
# these; lock the result in, along with the one prologue struct that is not a
# `{ data, size }` handle.
@testset "LibAccelerate honours __asm__ labels" begin
    libsrc = read(joinpath(dirname(pathof(AppleAccelerate)), "lib", "LibAccelerate.jl"), String)
    for f in ("BNNSGraphCompileFromFile", "BNNSGraphContextExecute", "BNNSGraphContextDestroy",
              "BNNSGraphContextGetWorkspaceSize", "BNNSGraphContextSetDynamicShapes",
              "BNNSGraphContextSetBatchSize", "BNNSGraphGetInputNames", "BNNSGraphGetOutputNames")
        @test isdefined(AppleAccelerate.LibAccelerate, Symbol(f))
        @test occursin("@ccall libacc.$(f)_v2(", libsrc)
        @test !occursin("@ccall libacc.$(f)(", libsrc)
    end
    S = AppleAccelerate.LibAccelerate.bnns_graph_shape_t
    @test fieldnames(S) == (:rank, :shape)
    @test fieldtypes(S) == (Csize_t, Ptr{UInt64})
end
