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
