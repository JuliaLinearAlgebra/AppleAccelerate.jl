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

end # testset vImage
