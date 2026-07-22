# vImage wrapper tests. Cross-validated against plain-Julia references where a
# reference is easy; every operation family is exercised. Images are deliberately
# NON-SQUARE so a width/height/rowBytes transpose bug fails loudly.

using AppleAccelerate, Test
const V = AppleAccelerate

# Non-square dimensions used throughout: Julia size (width, height) for planar.
const W, H = 7, 5

# Convenience constructors matching the module's layout convention.
planar8(w=W, h=H)  = rand(UInt8, w, h)
planarF(w=W, h=H)  = rand(Float32, w, h)
argb8(w=W, h=H)    = rand(UInt8, 4, w, h)
argbF(w=W, h=H)    = rand(Float32, 4, w, h)

@testset "vImage" begin

@testset "buffer / layout convention" begin
    A = reshape(Float32.(1:W*H), W, H)
    buf = V.vimage_buffer(A)
    @test buf.width == W          # first Julia dim == vImage width
    @test buf.height == H
    @test buf.rowBytes == W * sizeof(Float32)
    I = rand(UInt8, 4, W, H)
    bi = V.vimage_buffer(I)
    @test bi.width == W && bi.height == H
    @test bi.rowBytes == 4 * W * sizeof(UInt8)
    @test_throws V.vImageError throw(V.vImageError(-21772, "test"))
    @test occursin("NULL", V.vimage_error_string(-21772))
end

@testset "geometry: scale identity" begin
    for img in (planarF(), argbF())
        out = ndims(img) == 2 ? V.scale_PlanarF(img, W, H) : V.scale_ARGBFFFF(img, W, H)
        @test size(out) == size(img)
        @test out ≈ img rtol = 1e-4
    end
    # scale to a different (non-square) size runs and gives right shape
    big = V.scale_Planar8(planar8(), 11, 3)
    @test size(big) == (11, 3)
end

@testset "geometry: reflect (checks correct axis)" begin
    img = reshape(Float32.(1:W*H), W, H)
    @test V.horizontalReflect_PlanarF(img) ≈ reverse(img, dims = 1)
    @test V.verticalReflect_PlanarF(img)   ≈ reverse(img, dims = 2)
    a = argb8()
    @test V.horizontalReflect_ARGB8888(a) == reverse(a, dims = 2)  # width is dim 2 interleaved
end

@testset "geometry: rotate90 dims + rotate/shear/perspective/affine run" begin
    img = reshape(Float32.(1:W*H), W, H)
    r0 = V.rotate90_PlanarF(img, 0)   # no rotation
    @test size(r0) == (W, H)
    @test r0 ≈ img
    r4 = V.rotate90_PlanarF(img, 4)   # 4 & 3 == 0, wrapper masks the count
    @test r4 ≈ img
    r2 = V.rotate90_PlanarF(img, 2)   # 180° == reverse both axes
    @test size(r2) == (W, H)
    @test r2 ≈ reverse(reverse(img, dims = 1), dims = 2)
    r1 = V.rotate90_PlanarF(img, 1)   # 90° swaps dims
    @test size(r1) == (H, W)
    # affine identity reproduces input
    Id = V.vImage_AffineTransform(1, 0, 0, 1, 0, 0)
    @test V.affineWarp_PlanarF(img, Id; flags = V.kvImageBackgroundColorFill) ≈ img rtol = 1e-4
    IdD = V.vImage_AffineTransform_Double(1, 0, 0, 1, 0, 0)
    @test V.affineWarpD_PlanarF(img, IdD; flags = V.kvImageBackgroundColorFill) ≈ img rtol = 1e-4
    @test V.affineWarpCG_PlanarF(img, IdD; flags = V.kvImageBackgroundColorFill) ≈ img rtol = 1e-4
    # rotate / shear / perspective just need to run and give right shape
    @test size(V.rotate_Planar8(planar8(), 0.0f0)) == (W, H)
    @test size(V.horizontalShear_Planar8(planar8(), 0.0f0, 0.0f0)) == (W, H)
    @test size(V.verticalShearD_PlanarF(planarF(), 0.0, 0.0)) == (W, H)
    P = V.vImage_PerspectiveTransform(1, 0, 0, 1, 0, 0, 0, 0, 1)
    @test size(V.perspectiveWarp_Planar8(planar8(), P)) == (W, H)
end

@testset "morphology: constant image unchanged; vs manual max/min" begin
    c = fill(UInt8(50), W, H)
    ker = ones(UInt8, 3, 3)
    @test all(==(UInt8(50)), V.dilate_Planar8(c, ker))
    @test all(==(UInt8(50)), V.erode_Planar8(c, ker))
    @test all(==(UInt8(50)), V.max_Planar8(c, 3, 3))
    @test all(==(UInt8(50)), V.min_Planar8(c, 3, 3))
    # manual 3x3 max filter on interior vs vImage Max (edge-extended)
    img = rand(UInt8, W, H)
    mx = V.max_Planar8(img, 3, 3; flags = V.kvImageEdgeExtend)
    for x in 2:W-1, y in 2:H-1
        @test mx[x, y] == maximum(@view img[x-1:x+1, y-1:y+1])
    end
    mn = V.min_Planar8(img, 3, 3; flags = V.kvImageEdgeExtend)
    for x in 2:W-1, y in 2:H-1
        @test mn[x, y] == minimum(@view img[x-1:x+1, y-1:y+1])
    end
end

@testset "convolution: identity kernel, box/tent constant, sep" begin
    img = reshape(UInt8.(0:UInt8(W*H-1)), W, H)
    idk = Int16[0 0 0; 0 1 0; 0 0 0]
    @test V.convolve_Planar8(img, idk; flags = V.kvImageEdgeExtend) == img
    @test V.convolveWithBias_Planar8(img, idk; bias = 0, flags = V.kvImageEdgeExtend) == img
    c = fill(UInt8(40), W, H)
    @test all(==(UInt8(40)), V.boxConvolve_Planar8(c, 3, 3; flags = V.kvImageEdgeExtend))
    @test all(==(UInt8(40)), V.tentConvolve_Planar8(c, 3, 3; flags = V.kvImageEdgeExtend))
    # float convolve identity
    fimg = planarF()
    idf = Float32[0 0 0; 0 1 0; 0 0 0]
    @test V.convolve_PlanarF(fimg, idf; flags = V.kvImageEdgeExtend) ≈ fimg
    # separable identity (1-tap kernels)
    sep = V.sepConvolve_PlanarF(fimg, Float32[1], Float32[1]; flags = V.kvImageEdgeExtend)
    @test sep ≈ fimg rtol = 1e-4
end

@testset "histogram" begin
    hi = fill(UInt8(7), 4, 6)
    h = V.histogramCalculation_Planar8(hi)
    @test length(h) == 256
    @test h[8] == 24 && sum(h) == 24     # bin 7 (0-based) -> index 8
    ha = V.histogramCalculation_ARGB8888(fill(UInt8(3), 4, 4, 6))
    @test size(ha) == (256, 4)
    @test all(ha[4, :] .== 24)
    hf = V.histogramCalculation_PlanarF(fill(0.5f0, 4, 6), 10, 0f0, 1f0)
    @test sum(hf) == 24
    # equalization / contrast stretch run and preserve shape
    @test size(V.equalization_Planar8(planar8())) == (W, H)
    @test size(V.contrastStretch_Planar8(planar8())) == (W, H)
    @test size(V.endsInContrastStretch_Planar8(planar8(); percentLow = 1, percentHigh = 1)) == (W, H)
    @test size(V.equalization_PlanarF(planarF())) == (W, H)
end

@testset "alpha: premultiply/unpremultiply round trip; blends" begin
    a = argb8(); a[1, :, :] .= 0xff        # opaque alpha => exact round trip
    pm = V.premultiplyData_ARGB8888(a)
    @test V.unpremultiplyData_ARGB8888(pm) == a
    af = argbF(); af[1, :, :] .= 1f0
    @test V.unpremultiplyData_ARGBFFFF(V.premultiplyData_ARGBFFFF(af)) ≈ af rtol = 1e-5
    @test size(V.clipToAlpha_ARGB8888(a)) == (4, W, H)
    # planar premultiply
    src = planar8(); alpha = fill(0xff, W, H)
    @test V.premultiplyData_Planar8(src, alpha) == src
    # blends run
    top = argb8(); bot = argb8()
    @test size(V.premultipliedAlphaBlend_ARGB8888(top, bot)) == (4, W, H)
    @test size(V.alphaBlend_ARGB8888(top, bot)) == (4, W, H)
    rgba = rand(UInt8, 4, W, H)
    @test size(V.premultipliedAlphaBlendScreen_RGBA8888(rgba, rand(UInt8, 4, W, H))) == (4, W, H)
end

@testset "transform: matrix multiply, gamma, lookup, floodfill" begin
    # identity 4x4 colour matrix (divisor 256 -> *1)
    a = argb8()
    Imat = Int16[256 0 0 0; 0 256 0 0; 0 0 256 0; 0 0 0 256]
    out = V.matrixMultiply_ARGB8888(a, Imat; divisor = 256)
    @test out == a
    # RGB -> single plane weighted sum
    p = V.matrixMultiply_ARGB8888ToPlanar8!(similar(a[1, :, :]), a, Int16[0, 0, 0, 256]; divisor = 256)
    @test p == a[4, :, :]
    # identity piecewise gamma (linear a=1,b=0, gamma=1, boundary huge so linear branch used)
    fimg = rand(Float32, W, H)
    g = V.piecewiseGamma_PlanarF(fimg; exponentialCoeffs = (1f0, 0f0, 0f0), gamma = 1f0,
                                 linearCoeffs = (1f0, 0f0), boundary = 2f0)
    @test g ≈ fimg rtol = 1e-5
    # lookup table: identity table
    src = reshape(UInt8.(0:UInt8(W*H-1)), W, H)
    tbl = UInt16.(0:255)
    dst = Matrix{UInt16}(undef, W, H)
    V.lookupTable_Planar8toPlanar16!(dst, src, tbl)
    @test dst == UInt16.(src)
    # flood fill
    canvas = zeros(UInt8, W, H)
    V.floodFill_Planar8!(canvas, 0, 0, 9; connectivity = 4)
    @test all(==(UInt8(9)), canvas)
end

@testset "conversion: channels" begin
    a = argb8()
    # permute reverse then reverse == identity
    rev = UInt8[3, 2, 1, 0]
    p = V.permuteChannels_ARGB8888(a, rev)
    @test V.permuteChannels_ARGB8888(p, rev) == a
    @test p[1, :, :] == a[4, :, :]
    # extract channel matches manual slice
    @test V.extractChannel_ARGB8888(a, 2) == a[3, :, :]
    # buffer fill
    d = Array{UInt8,3}(undef, 4, W, H)
    V.bufferFill_ARGB8888!(d, UInt8[10, 20, 30, 40])
    @test all(d[:, x, y] == UInt8[10, 20, 30, 40] for x in 1:W, y in 1:H)
    # overwrite planar with scalar
    pl = planar8()
    V.overwriteChannelsWithScalar_Planar8!(pl, UInt8(77))
    @test all(==(UInt8(77)), pl)
    # table lookup planar identity
    src = reshape(UInt8.(0:UInt8(W*H-1)), W, H)
    @test V.tableLookUp_Planar8(src, UInt8.(0:255)) == src
    # deinterleave then interleave round trips
    A, R, G, B = (Matrix{UInt8}(undef, W, H) for _ in 1:4)
    V.convert_ARGB8888toPlanar8(a, A, R, G, B)
    @test R == a[2, :, :]
    dest = Array{UInt8,3}(undef, 4, W, H)
    V.convert_Planar8toARGB8888(dest, A, R, G, B)
    @test dest == a
end

@testset "conversion: format / bit depth" begin
    # Planar8 <-> PlanarF round trip
    p8 = reshape(UInt8.(0:UInt8(W*H-1)), W, H)
    pf = V.convert_Planar8toPlanarF(p8; maxFloat = 255f0, minFloat = 0f0)
    @test pf ≈ Float32.(p8)
    @test V.convert_PlanarFtoPlanar8(pf; maxFloat = 255f0, minFloat = 0f0) == p8
    # Planar8 <-> 16U round trip (shape checks)
    d16 = Matrix{UInt16}(undef, W, H)
    V.convert_Planar8To16U!(d16, p8)
    back = Matrix{UInt8}(undef, W, H)
    V.convert_16UToPlanar8!(back, d16)
    @test back == p8
    # ARGB8888 -> RGB888 drops the alpha (first) channel
    a = argb8()
    rgb = Array{UInt8,3}(undef, 3, W, H)
    V.convert_ARGB8888toRGB888!(rgb, a)
    @test rgb == a[2:4, :, :]
    # byte swap is its own inverse
    x = rand(UInt16, W, H); y = similar(x); z = similar(x)
    V.byteSwap_Planar16U!(y, x); V.byteSwap_Planar16U!(z, y)
    @test z == x
    # copyBuffer
    dst = similar(a)
    V.copyBuffer!(dst, a; pixelSize = 4)
    @test dst == a
    # flatten over opaque background: alpha=255 opaque => rgb unchanged
    ao = argb8(); ao[1, :, :] .= 0xff
    fl = Array{UInt8,3}(undef, 3, W, H)
    V.flatten_ARGB8888ToRGB888!(fl, ao, UInt8[0, 0, 0, 0]; isImagePremultiplied = true)
    @test fl == ao[2:4, :, :]
    # RGB101010 round trip: 8888 -> packed -> 8888
    packed = Matrix{UInt32}(undef, W, H)
    a2 = argb8()
    V.convert_ARGB8888ToARGB2101010!(packed, a2; rangeMin = 0, rangeMax = 1023)
    back2 = Array{UInt8,3}(undef, 4, W, H)
    V.convert_ARGB2101010ToARGB8888!(back2, packed; rangeMin = 0, rangeMax = 1023)
    @test back2[2:4, :, :] == a2[2:4, :, :]   # 2101010 keeps 10-bit RGB exactly; alpha is only 2-bit
end

@testset "error handling" begin
    # bad argument shapes throw a Julia error before the ccall
    @test_throws DimensionMismatch V.permuteChannels_ARGB8888(argb8(), UInt8[0, 1])
    @test_throws DimensionMismatch V.matrixMultiply_ARGB8888(argb8(), Int16[1 2; 3 4])
    @test_throws DimensionMismatch V.tableLookUp_Planar8(planar8(), UInt8[1, 2, 3])
end

@testset "YpCbCr <-> ARGB 4:4:4" begin
    pr = V.kvImageYpCbCrPixelRange_FullRange_8bit_Clamped
    toYpCbCr = V.convert_ARGBToYpCbCr_GenerateConversion(
        V.kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4, pr, V.kvImageARGB8888, V.kvImage444CrYpCb8)
    toARGB = V.convert_YpCbCrToARGB_GenerateConversion(
        V.kvImage_YpCbCrToARGBMatrix_ITU_R_601_4, pr, V.kvImage444CrYpCb8, V.kvImageARGB8888)
    @test toYpCbCr isa V.vImage_ARGBToYpCbCr
    @test toARGB isa V.vImage_YpCbCrToARGB

    # Round-trip through v308 (3-channel Cr,Y,Cb; no chroma subsampling => near-exact RGB)
    src = argb8(9, 6)
    ycc = Array{UInt8,3}(undef, 3, 9, 6)
    back = Array{UInt8,3}(undef, 4, 9, 6)
    @test V.convert_ARGB8888To444CrYpCb8!(ycc, src, toYpCbCr) === ycc
    V.convert_444CrYpCb8ToARGB8888!(back, ycc, toARGB; alpha = 0xff)
    @test size(ycc) == (3, 9, 6)
    @test maximum(abs.(Int.(src[2:4, :, :]) .- Int.(back[2:4, :, :]))) <= 3

    # 4-channel v408 (Cb,Y,Cr,A) round trip carries alpha through
    from408 = V.convert_ARGBToYpCbCr_GenerateConversion(
        V.kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4, pr, V.kvImageARGB8888, V.kvImage444CbYpCrA8)
    to408 = V.convert_YpCbCrToARGB_GenerateConversion(
        V.kvImage_YpCbCrToARGBMatrix_ITU_R_601_4, pr, V.kvImage444CbYpCrA8, V.kvImageARGB8888)
    y408 = Array{UInt8,3}(undef, 4, 9, 6)
    b408 = Array{UInt8,3}(undef, 4, 9, 6)
    V.convert_ARGB8888To444CbYpCrA8!(y408, src, from408)
    V.convert_444CbYpCrA8ToARGB8888!(b408, y408, to408)
    @test maximum(abs.(Int.(src[2:4, :, :]) .- Int.(b408[2:4, :, :]))) <= 3
    @test b408[1, :, :] == src[1, :, :]           # alpha preserved exactly

    # 16-bit y416 conversions run and produce the right shape
    from16 = V.convert_ARGBToYpCbCr_GenerateConversion(
        V.kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4, pr, V.kvImageARGB8888, V.kvImage444AYpCbCr16)
    y416 = Array{UInt16,3}(undef, 4, 9, 6)
    V.convert_ARGB8888To444AYpCbCr16!(y416, src, from16)
    @test size(y416) == (4, 9, 6)
    pr16 = V.vImage_YpCbCrPixelRange(0, 32768, 65535, 65535, 65535, 1, 65535, 0)  # full-range 16-bit
    to16 = V.convert_YpCbCrToARGB_GenerateConversion(
        V.kvImage_YpCbCrToARGBMatrix_ITU_R_601_4, pr16, V.kvImage444AYpCbCr16, V.kvImageARGB16U)
    argb16 = Array{UInt16,3}(undef, 4, 9, 6)
    V.convert_444AYpCbCr16ToARGB16U!(argb16, y416, to16)
    @test size(argb16) == (4, 9, 6)
end

@testset "multi-kernel convolution" begin
    src = argb8(9, 6)
    k = Int16[0 1 0; 1 4 1; 0 1 0]
    # identical kernel per channel == single-kernel convolve (same divisor, edging)
    single = V.convolve_ARGB8888(src, k; divisor = 8, flags = V.kvImageEdgeExtend)
    multi = V.convolveMultiKernel_ARGB8888(src, (k, k, k, k);
        divisors = (8, 8, 8, 8), biases = (0, 0, 0, 0), flags = V.kvImageEdgeExtend)
    @test size(multi) == size(src)
    # identical per-channel kernels reproduce the single-kernel result (up to a ±1 rounding
    # difference in the two code paths' final divide)
    @test maximum(abs.(Int.(multi) .- Int.(single))) <= 1

    # float multi-kernel runs and matches the float single-kernel convolve
    srcF = argbF(9, 6)
    kf = Float32[0 1 0; 1 4 1; 0 1 0] ./ 8
    singleF = V.convolve_ARGBFFFF(srcF, kf; flags = V.kvImageEdgeExtend)
    multiF = V.convolveMultiKernel_ARGBFFFF(srcF, (kf, kf, kf, kf); flags = V.kvImageEdgeExtend)
    @test multiF ≈ singleF rtol = 1e-4

    # float-kernel ARGB8888 runs and gives correct shape
    fk = V.convolveFloatKernel_ARGB8888(src, kf; flags = V.kvImageEdgeExtend)
    @test size(fk) == size(src)

    @test_throws DimensionMismatch V.convolveMultiKernel_ARGB8888(src, (k, k, k))
end

@testset "multi-plane matrix multiply" begin
    # RGB planes -> single luminance plane: out = 0.25 R + 0.5 G + 0.25 B
    r = rand(UInt8, 8, 5); g = rand(UInt8, 8, 5); b = rand(UInt8, 8, 5)
    y = Matrix{UInt8}(undef, 8, 5)
    M = reshape(Int16[64, 128, 64], 3, 1)          # (src_planes=3, dest_planes=1)
    V.matrixMultiply_Planar8([r, g, b], [y], M; divisor = 256)
    ref = round.(UInt8, (64 .* Int.(r) .+ 128 .* Int.(g) .+ 64 .* Int.(b) .+ 128) .÷ 256)
    @test maximum(abs.(Int.(y) .- Int.(ref))) <= 1

    # float form: exact matmul. out[j] = Σ_i in[i]·M[i,j]
    rf = rand(Float32, 8, 5); gf = rand(Float32, 8, 5)
    o1 = Matrix{Float32}(undef, 8, 5); o2 = Matrix{Float32}(undef, 8, 5)
    Mf = Float32[1.0 0.0; 0.5 0.5]                  # (2 src, 2 dest); in1->{o1}, in2->{o1/2,o2/2}
    V.matrixMultiply_PlanarF([rf, gf], [o1, o2], Mf)
    @test o1 ≈ (rf .+ 0.5f0 .* gf) rtol = 1e-5
    @test o2 ≈ (0.5f0 .* gf) rtol = 1e-5
    @test_throws DimensionMismatch V.matrixMultiply_PlanarF([rf, gf], [o1], Float32[1 0 0; 0 1 0])
end

@testset "additional alpha compositing" begin
    # premultiplied planar source-over: dest = top + (1-topA)*bottom.
    # With a zero bottom the result is just the (premultiplied) top.
    top = planar8(8, 5); topA = planar8(8, 5); bot = planar8(8, 5)
    out0 = V.premultipliedAlphaBlend_Planar8(top, topA, zeros(UInt8, 8, 5))
    @test out0 == top
    # With an opaque top (alpha 255) the bottom contributes nothing either.
    outO = V.premultipliedAlphaBlend_Planar8(top, fill(0xff, 8, 5), planar8(8, 5))
    @test outO == top
    @test size(out0) == (8, 5)

    # const-alpha interleaved blend runs and returns correct shape
    ct = argb8(8, 5); cb = argb8(8, 5)
    co = V.premultipliedConstAlphaBlend_ARGB8888(ct, 0x80, cb)
    @test size(co) == (4, 8, 5)

    # const-alpha planar blend
    cpo = V.premultipliedConstAlphaBlend_Planar8(top, 0x80, topA, bot)
    @test size(cpo) == (8, 5)

    # nonpremultiplied->premultiplied planar
    npo = V.alphaBlend_NonpremultipliedToPremultiplied_Planar8(top, topA, bot)
    @test size(npo) == (8, 5)

    # float planar 6-buffer blend
    fo = V.alphaBlend_PlanarF(planarF(8, 5), planarF(8, 5), planarF(8, 5), planarF(8, 5), planarF(8, 5))
    @test size(fo) == (8, 5)

    # with-permute premultiplied blend
    wp = V.premultipliedAlphaBlendWithPermute_ARGB8888(argb8(8, 5), argb8(8, 5); makeDestAlphaOpaque = true)
    @test size(wp) == (4, 8, 5)
    @test all(wp[1, :, :] .== 0xff)               # dest alpha forced opaque
end

@testset "histogram specification (float / interleaved)" begin
    src = planarF(8, 5)
    flat = fill(1, 256)                            # uniform desired histogram
    out = V.histogramSpecification_PlanarF(src, flat; entries = 256, minVal = 0f0, maxVal = 1f0)
    @test size(out) == (8, 5)

    srcA = argb8(8, 5)
    hists = [fill(Csize_t(1), 256) for _ in 1:4]
    outA = V.histogramSpecification_ARGB8888(srcA, hists)
    @test size(outA) == (4, 8, 5)

    srcAF = argbF(8, 5)
    histsF = [fill(Csize_t(1), 256) for _ in 1:4]
    outAF = V.histogramSpecification_ARGBFFFF(srcAF, histsF; entries = 256, minVal = 0f0, maxVal = 1f0)
    @test size(outAF) == (4, 8, 5)

    # ends-in contrast stretch, float forms
    es = V.endsInContrastStretch_PlanarF(planarF(8, 5); percentLow = 1, percentHigh = 1)
    @test size(es) == (8, 5)
    esA = V.endsInContrastStretch_ARGBFFFF(argbF(8, 5); percentLow = (1, 1, 1, 1), percentHigh = (1, 1, 1, 1))
    @test size(esA) == (4, 8, 5)
end

@testset "additional plane <-> interleaved conversions" begin
    # deinterleave 4-channel -> 3 planes dropping X. XRGB: src channels [X,R,G,B]
    src = argb8(8, 5)
    r = Matrix{UInt8}(undef, 8, 5); g = similar(r); b = similar(r)
    V.convert_XRGB8888ToPlanar8(src, r, g, b)
    @test r == src[2, :, :] && g == src[3, :, :] && b == src[4, :, :]
    # BGRX: dest order is (blue, green, red); src channels [B,G,R,X]
    rb = similar(r); gb = similar(r); bb = similar(r)
    V.convert_BGRX8888ToPlanar8(src, bb, gb, rb)
    @test bb == src[1, :, :] && gb == src[2, :, :] && rb == src[3, :, :]

    # 4 planes -> ARGBFFFF interleave (identity [0,1] mapping)
    a = planar8(8, 5); rr = planar8(8, 5); gg = planar8(8, 5); bb2 = planar8(8, 5)
    dst = Array{Float32,3}(undef, 4, 8, 5)
    V.convert_Planar8ToARGBFFFF!(dst, a, rr, gg, bb2; maxFloat = (255f0, 255f0, 255f0, 255f0))
    @test dst[1, :, :] ≈ Float32.(a) rtol = 1e-4
    @test dst[2, :, :] ≈ Float32.(rr) rtol = 1e-4

    # Planar16Q12 <-> ARGB8888 round trip through Int16 Q4.12 planes
    src8 = argb8(8, 5)
    aq = Matrix{Int16}(undef, 8, 5); rq = similar(aq); gq = similar(aq); bq = similar(aq)
    V.convert_ARGB8888toPlanar16Q12!(aq, rq, gq, bq, src8)
    back = Array{UInt8,3}(undef, 4, 8, 5)
    V.convert_Planar16Q12toARGB8888!(back, aq, rq, gq, bq)
    @test maximum(abs.(Int.(back) .- Int.(src8))) <= 1
end

end # testset vImage
