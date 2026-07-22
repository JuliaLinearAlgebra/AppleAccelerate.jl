# vImage — Apple's image-processing subframework, hand-written idiomatic wrappers.
#
# vImage is a large, extremely regular C API: each *operation* (Scale, Rotate,
# Convolve, Dilate, Premultiply, PermuteChannels, format Convert, …) is provided
# for ~a dozen pixel formats (Planar8, PlanarF, ARGB8888, ARGBFFFF, …). There are
# ~480 operation functions. Clang.jl cannot parse the vImage headers (they drag in
# CoreGraphics/CoreVideo and an unparseable arm_neon.h), so — exactly like the C++
# libSparse symbols in sparse.jl — these bindings are written by hand with `ccall`
# straight against the Accelerate framework binary.
#
# The bindings are generated with `@eval` loops over (operation, pixel-format) so a
# handful of loops cover hundreds of C functions.
#
# ---------------------------------------------------------------------------
# MEMORY / LAYOUT CONVENTION  (verified empirically with non-square test images)
# ---------------------------------------------------------------------------
# vImage is ROW-major and describes a buffer by {data, height, width, rowBytes}.
# Julia is COLUMN-major. We map a vImage image onto a Julia array so that the
# in-memory byte order matches vImage's exactly (no transpose needed):
#
#   * PLANAR (single channel) image  -> `Matrix{T}` of size `(width, height)`.
#       vImage width  = size(A, 1)   (the fast / column dimension)
#       vImage height = size(A, 2)
#       rowBytes      = size(A, 1) * sizeof(T)
#
#   * INTERLEAVED (multi-channel) image -> `Array{T,3}` of size `(channels, width, height)`.
#       vImage width  = size(A, 2)
#       vImage height = size(A, 3)
#       rowBytes      = size(A, 1) * size(A, 2) * sizeof(T)
#     The channel dimension is first (fastest), matching interleaved pixel storage
#     e.g. ARGB8888 = `Array{UInt8,3}` with size(A,1) == 4.
#
# So: a `w×h` Julia matrix is a `w`-wide, `h`-tall vImage image. Because vImage's
# width is our first (fast) Julia dimension, NO transpose is required and the data
# can be shared copy-free.
#
# Apple docs: https://developer.apple.com/documentation/accelerate/vimage

# vImage symbols live in the Accelerate framework binary.
const vimage_lib = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

# --- Core types --------------------------------------------------------------

"""
    vImage_Buffer

Mirror of the C `vImage_Buffer` descriptor: `data::Ptr{Cvoid}`, `height`, `width`,
`rowBytes` (all `Csize_t`). Built from a Julia array by [`vimage_buffer`](@ref) and
passed to vImage functions by `Ref`. The struct only borrows the array's pointer, so
the backing array must be kept alive (`GC.@preserve`) for the duration of any call.
"""
struct vImage_Buffer
    data::Ptr{Cvoid}
    height::Csize_t
    width::Csize_t
    rowBytes::Csize_t
end

"""
    vImage_AffineTransform

Single-precision 3×2 affine transform (`a b; c d; tx ty`) used by
`affineWarp_*` geometry functions. See also [`vImage_AffineTransform_Double`](@ref).
"""
struct vImage_AffineTransform
    a::Cfloat; b::Cfloat; c::Cfloat; d::Cfloat; tx::Cfloat; ty::Cfloat
end

"""
    vImage_AffineTransform_Double

Double-precision variant of [`vImage_AffineTransform`](@ref), used by the
`affineWarpD_*` functions.
"""
struct vImage_AffineTransform_Double
    a::Cdouble; b::Cdouble; c::Cdouble; d::Cdouble; tx::Cdouble; ty::Cdouble
end

"""
    vImage_PerspectiveTransform

Single-precision 3×3 projective transform used by `perspectiveWarp_*` functions.
"""
struct vImage_PerspectiveTransform
    a::Cfloat; b::Cfloat; c::Cfloat; d::Cfloat; tx::Cfloat; ty::Cfloat
    vx::Cfloat; vy::Cfloat; v::Cfloat
end

# vImage_Error is a pointer-sized signed integer; 0 = success, <0 = error,
# >0 = temp-buffer size (when kvImageGetTempBufferSize is set).
const vImage_Error = Cssize_t
const vImage_Flags = UInt32

# --- Error codes -------------------------------------------------------------

const kvImageNoError                     = 0
const kvImageRoiLargerThanInputBuffer    = -21766
const kvImageInvalidKernelSize           = -21767
const kvImageInvalidEdgeStyle            = -21768
const kvImageInvalidOffset_X             = -21769
const kvImageInvalidOffset_Y             = -21770
const kvImageMemoryAllocationError       = -21771
const kvImageNullPointerArgument         = -21772
const kvImageInvalidParameter            = -21773
const kvImageBufferSizeMismatch          = -21774
const kvImageUnknownFlagsBit             = -21775
const kvImageInternalError               = -21776
const kvImageInvalidRowBytes             = -21777
const kvImageInvalidImageFormat          = -21778
const kvImageColorSyncIsAbsent           = -21779
const kvImageOutOfPlaceOperationRequired = -21780
const kvImageInvalidImageObject          = -21781
const kvImageInvalidCVImageFormat        = -21782
const kvImageUnsupportedConversion       = -21783
const kvImageCoreVideoIsAbsent           = -21784

const _VIMAGE_ERRSTR = Dict{Int,String}(
    kvImageNoError => "no error",
    kvImageRoiLargerThanInputBuffer => "ROI larger than input buffer (destination bigger than source?)",
    kvImageInvalidKernelSize => "invalid kernel size",
    kvImageInvalidEdgeStyle => "invalid edge style flag",
    kvImageInvalidOffset_X => "invalid X offset",
    kvImageInvalidOffset_Y => "invalid Y offset",
    kvImageMemoryAllocationError => "memory allocation failed",
    kvImageNullPointerArgument => "unexpected NULL pointer argument",
    kvImageInvalidParameter => "invalid parameter",
    kvImageBufferSizeMismatch => "buffer size mismatch (planar buffers differ in size?)",
    kvImageUnknownFlagsBit => "unknown/unsupported flag bit set",
    kvImageInternalError => "internal vImage error (please file a bug)",
    kvImageInvalidRowBytes => "invalid rowBytes",
    kvImageInvalidImageFormat => "invalid image format",
    kvImageColorSyncIsAbsent => "ColorSync is absent",
    kvImageOutOfPlaceOperationRequired => "out-of-place operation required",
    kvImageInvalidImageObject => "invalid image object",
    kvImageInvalidCVImageFormat => "invalid CV image format",
    kvImageUnsupportedConversion => "unsupported conversion",
    kvImageCoreVideoIsAbsent => "CoreVideo is absent",
)

"""
    vimage_error_string(code) -> String

Human-readable description of a vImage error code (a negative `vImage_Error`).
"""
vimage_error_string(code::Integer) = get(_VIMAGE_ERRSTR, Int(code), "unknown error code $code")

# --- Flags -------------------------------------------------------------------

const kvImageNoFlags                = vImage_Flags(0)
const kvImageLeaveAlphaUnchanged    = vImage_Flags(1)
const kvImageCopyInPlace            = vImage_Flags(2)
const kvImageBackgroundColorFill    = vImage_Flags(4)
const kvImageEdgeExtend             = vImage_Flags(8)
const kvImageDoNotTile              = vImage_Flags(16)
const kvImageHighQualityResampling  = vImage_Flags(32)
const kvImageTruncateKernel         = vImage_Flags(64)
const kvImageGetTempBufferSize      = vImage_Flags(128)
const kvImagePrintDiagnosticsToConsole = vImage_Flags(256)
const kvImageNoAllocate             = vImage_Flags(512)
const kvImageHDRContent             = vImage_Flags(1024)
const kvImageDoNotClamp             = vImage_Flags(2048)
const kvImageUseFP16Accumulator     = vImage_Flags(4096)

# --- Pixel type aliases (documentation / convenience) ------------------------

const Pixel_8       = UInt8
const Pixel_F       = Float32
const Pixel_16U     = UInt16
const Pixel_16S     = Int16
const Pixel_16Q12   = Int16
const Pixel_8888    = NTuple{4,UInt8}
const Pixel_FFFF    = NTuple{4,Float32}
const Pixel_ARGB_16U = NTuple{4,UInt16}
const Pixel_ARGB_16S = NTuple{4,Int16}

# --- Exceptions --------------------------------------------------------------

"""
    vImageError <: Exception

Thrown when a vImage function returns a negative `vImage_Error`. Carries the numeric
`code` and the C function name `fn`; see [`vimage_error_string`](@ref
AppleAccelerate.vimage_error_string) for the human-readable meaning.
"""
struct vImageError <: Exception
    code::Int
    fn::String
end
function Base.showerror(io::IO, e::vImageError)
    print(io, "vImageError: ", e.fn, " returned ", e.code, " (", vimage_error_string(e.code), ")")
end

@inline function _check(err::Integer, fn)
    err < 0 && throw(vImageError(Int(err), String(fn)))
    return err
end

# --- Buffer construction -----------------------------------------------------

"""
    vimage_buffer(A) -> vImage_Buffer

Build a [`vImage_Buffer`](@ref) descriptor that borrows the memory of Julia array `A`.

* `A::AbstractMatrix` is treated as a **planar** (single-channel) image of size
  `(width, height) == (size(A,1), size(A,2))`.
* `A::AbstractArray{T,3}` is treated as an **interleaved** image whose first
  dimension is the channel count, i.e. size `(channels, width, height)`.

The descriptor only stores a raw pointer; keep `A` alive with `GC.@preserve` while
the descriptor (or any `Ref` of it) is in use. See the layout note at the top of
`vimage.jl`.
"""
@inline function vimage_buffer(A::AbstractMatrix)
    T = eltype(A)
    return vImage_Buffer(Ptr{Cvoid}(pointer(A)), Csize_t(size(A, 2)),
                         Csize_t(size(A, 1)), Csize_t(size(A, 1) * sizeof(T)))
end
@inline function vimage_buffer(A::AbstractArray{T,3}) where {T}
    return vImage_Buffer(Ptr{Cvoid}(pointer(A)), Csize_t(size(A, 3)),
                         Csize_t(size(A, 2)), Csize_t(size(A, 1) * size(A, 2) * sizeof(T)))
end

# Per-pixel-format metadata: format-suffix => (Julia eltype, channel count).
const _VFMT = Dict{Symbol,Tuple{DataType,Int}}(
    :Planar8   => (UInt8, 1),  :PlanarF   => (Float32, 1), :Planar16U => (UInt16, 1),
    :Planar16S => (Int16, 1),  :Planar16F => (UInt16, 1),  :Planar16Q12 => (Int16, 1),
    :ARGB8888  => (UInt8, 4),  :ARGBFFFF  => (Float32, 4), :ARGB16U   => (UInt16, 4),
    :ARGB16S   => (Int16, 4),  :ARGB16F   => (UInt16, 4),  :ARGB16Q12 => (Int16, 4),
    :RGBA8888  => (UInt8, 4),  :RGBAFFFF  => (Float32, 4), :RGBA16U   => (UInt16, 4),
    :RGBA16F   => (UInt16, 4), :RGBA16Q12 => (Int16, 4),
    :BGRA8888  => (UInt8, 4),  :BGRAFFFF  => (Float32, 4), :BGRA16U   => (UInt16, 4),
    :BGRA16F   => (UInt16, 4),
    :CbCr8     => (UInt8, 2),  :CbCr16U   => (UInt16, 2),  :CbCr16S   => (Int16, 2),
    :CbCr16F   => (UInt16, 2),
)

# Julia function name from a C symbol: strip the leading "vImage" and lowercase
# the first character, e.g. :vImageScale_Planar8 -> :scale_Planar8.
_jlname(sym) = Symbol(lowercasefirst(string(sym)[7:end]))

# Allocate a fresh interleaved/planar array shaped like `src`.
_alloc_like(src::AbstractArray) = similar(src)

# =====================================================================================
# GEOMETRY
# =====================================================================================

# ---- Scale : (src, dest, tempBuffer, flags) -----------------------------------------
for sfx in (:Planar8, :Planar16S, :Planar16U, :PlanarF, :ARGB8888, :ARGB16U,
            :ARGB16S, :ARGBFFFF, :Planar16F, :ARGB16F, :CbCr16F, :CbCr16U, :CbCr8)
    sym  = Symbol("vImageScale_", sfx)
    bang = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray, src::AbstractArray; flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, vImage_Flags),
                    sb, db, C_NULL, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end

"""
    scale_Planar8(src, width, height; flags=kvImageNoFlags) -> Matrix
    scale_ARGB8888(src, width, height; flags=kvImageNoFlags) -> Array

Resample planar/interleaved image `src` to a new `width × height` (in vImage pixels),
allocating and returning the destination. Mutating `scale_*!(dest, src)` variants
scale into a caller-provided `dest` of the desired size. Available for formats
`Planar8/16S/16U/16F/F`, `ARGB8888/16U/16S/16F/FFFF`, `CbCr8/16U/16F`.
"""
scale_Planar8
for sfx in (:Planar8, :Planar16S, :Planar16U, :PlanarF, :ARGB8888, :ARGB16U,
            :ARGB16S, :ARGBFFFF, :Planar16F, :ARGB16F, :CbCr16F, :CbCr16U, :CbCr8)
    T, nch = _VFMT[sfx]
    alloc = _jlname(Symbol("vImageScale_", sfx))
    bang  = Symbol(alloc, "!")
    @eval begin
        function $alloc(src::AbstractArray, width::Integer, height::Integer; flags::Integer = kvImageNoFlags)
            dest = $(nch == 1 ? :(Matrix{$T}(undef, width, height)) :
                                 :(Array{$T,3}(undef, $nch, width, height)))
            return $bang(dest, src; flags = flags)
        end
    end
end

# ---- Horizontal / Vertical reflect : (src, dest, flags) -----------------------------
for op in (:HorizontalReflect, :VerticalReflect)
    for sfx in (:Planar8, :Planar16U, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S,
                :ARGBFFFF, :Planar16F, :CbCr16F, :ARGB16F)
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray; flags::Integer = kvImageNoFlags)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                        sb, db, vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractArray; flags::Integer = kvImageNoFlags) =
                $bang(_alloc_like(src), src; flags = flags)
        end
    end
end

"""
    horizontalReflect_Planar8(src; flags=kvImageNoFlags)
    verticalReflect_ARGB8888(src; flags=kvImageNoFlags)

Mirror an image across the vertical (`horizontalReflect_*`) or horizontal
(`verticalReflect_*`) axis. Allocating and `!` mutating variants; formats
`Planar8/16U/16F/F`, `ARGB8888/16U/16S/16F/FFFF`, `CbCr16F`.
"""
horizontalReflect_Planar8

# helper: does this format's back-color argument decay to a pointer (nch>1)?
# Emit the back-color prep + ccall for geometry ops that take a back color.

# ---- Rotate90 : (src, dest, rotationConstant, backColor, flags) ---------------------
for sfx in (:Planar8, :Planar16U, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S,
            :ARGBFFFF, :Planar16F, :CbCr16F, :ARGB16F)
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageRotate90_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    bctype = nch == 1 ? T : Ptr{T}
    if nch == 1
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray, rotationConstant::Integer;
                           backColor = zero($T), flags::Integer = kvImageNoFlags)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, UInt8, $T, vImage_Flags),
                        sb, db, UInt8(rotationConstant & 3), $T(backColor), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
        end
    else
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray, rotationConstant::Integer;
                           backColor = zeros($T, $nch), flags::Integer = kvImageNoFlags)
                bc = convert(Vector{$T}, collect(backColor))
                GC.@preserve src dest bc begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, UInt8, Ptr{$T}, vImage_Flags),
                        sb, db, UInt8(rotationConstant & 3), pointer(bc), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
        end
    end
    @eval begin
        function $alloc(src::AbstractArray, rotationConstant::Integer; kw...)
            r = Int(rotationConstant) & 3
            if ndims(src) == 2
                dest = isodd(r) ? similar(src, size(src, 2), size(src, 1)) : similar(src)
            else
                dest = isodd(r) ? similar(src, size(src, 1), size(src, 3), size(src, 2)) : similar(src)
            end
            return $bang(dest, src, rotationConstant; kw...)
        end
    end
end

"""
    rotate90_Planar8(src, rotationConstant; backColor, flags) -> array

Rotate by `rotationConstant * 90°` counter-clockwise (`rotationConstant ∈ 0:3`). The
allocating variant returns a correctly-shaped destination (dimensions swap for 90°/270°);
`rotate90_*!(dest, src, rotationConstant)` rotates into `dest`. `backColor` is a scalar
for planar formats or a length-4 vector for interleaved formats.
"""
rotate90_Planar8

# ---- Rotate (arbitrary angle) : (src, dest, tempBuffer, angle, backColor, flags) ----
for sfx in (:Planar8, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S, :ARGBFFFF,
            :Planar16F, :CbCr16F, :ARGB16F)
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageRotate_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    if nch == 1
        @eval function $bang(dest::AbstractArray, src::AbstractArray, angle::Real;
                             backColor = zero($T), flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Cfloat, $T, vImage_Flags),
                    sb, db, C_NULL, Cfloat(angle), $T(backColor), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    else
        @eval function $bang(dest::AbstractArray, src::AbstractArray, angle::Real;
                             backColor = zeros($T, $nch), flags::Integer = kvImageNoFlags)
            bc = convert(Vector{$T}, collect(backColor))
            GC.@preserve src dest bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Cfloat, Ptr{$T}, vImage_Flags),
                    sb, db, C_NULL, Cfloat(angle), pointer(bc), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
    @eval $alloc(src::AbstractArray, angle::Real; kw...) = $bang(_alloc_like(src), src, angle; kw...)
end

"""
    rotate_Planar8(src, angleInRadians; backColor, flags) -> array

Rotate `src` counter-clockwise by `angleInRadians` about the image centre, filling
exposed corners with `backColor`. Allocating (same-size canvas) and `!` mutating
variants. Formats `Planar8/16F/F`, `ARGB8888/16U/16S/16F/FFFF`, `CbCr16F`.
"""
rotate_Planar8

# ---- AffineWarp / AffineWarpD : (src, dest, temp, transform*, backColor, flags) -----
for (op, TT) in ((:AffineWarp, :vImage_AffineTransform), (:AffineWarpD, :vImage_AffineTransform_Double))
    for sfx in (:Planar8, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S, :ARGBFFFF,
                :Planar16F, :CbCr16F, :ARGB16F)
        T, nch = _VFMT[sfx]
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        if nch == 1
            @eval function $bang(dest::AbstractArray, src::AbstractArray, transform::$TT;
                                 backColor = zero($T), flags::Integer = kvImageBackgroundColorFill)
                tref = Ref(transform)
                GC.@preserve src dest tref begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{$TT}, $T, vImage_Flags),
                        sb, db, C_NULL, tref, $T(backColor), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
        else
            @eval function $bang(dest::AbstractArray, src::AbstractArray, transform::$TT;
                                 backColor = zeros($T, $nch), flags::Integer = kvImageBackgroundColorFill)
                tref = Ref(transform); bc = convert(Vector{$T}, collect(backColor))
                GC.@preserve src dest tref bc begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{$TT}, Ptr{$T}, vImage_Flags),
                        sb, db, C_NULL, tref, pointer(bc), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
        end
        @eval $alloc(src::AbstractArray, transform::$TT; kw...) =
            $bang(_alloc_like(src), src, transform; kw...)
    end
end

"""
    affineWarp_Planar8(src, transform::vImage_AffineTransform; backColor, flags) -> array
    affineWarpD_ARGB8888(src, transform::vImage_AffineTransform_Double; ...) -> array

Apply an affine warp. `affineWarp_*` uses a single-precision
[`vImage_AffineTransform`](@ref); `affineWarpD_*` uses the double-precision
[`vImage_AffineTransform_Double`](@ref). The identity transform reproduces the input.
Allocating (same-size) and `!` variants.
"""
affineWarp_Planar8

# ---- PerspectiveWarp : (src, dest, temp, transform*, interpolation, backColor, flags)
for sfx in (:Planar8, :Planar16U, :Planar16F, :ARGB8888, :ARGB16U, :ARGB16F)
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImagePerspectiveWarp_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    if nch == 1
        @eval function $bang(dest::AbstractArray, src::AbstractArray, transform::vImage_PerspectiveTransform;
                             interpolation::Integer = 1, backColor = zero($T),
                             flags::Integer = kvImageBackgroundColorFill)
            tref = Ref(transform)
            GC.@preserve src dest tref begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{vImage_PerspectiveTransform}, Int32, $T, vImage_Flags),
                    sb, db, C_NULL, tref, Int32(interpolation), $T(backColor), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    else
        @eval function $bang(dest::AbstractArray, src::AbstractArray, transform::vImage_PerspectiveTransform;
                             interpolation::Integer = 1, backColor = zeros($T, $nch),
                             flags::Integer = kvImageBackgroundColorFill)
            tref = Ref(transform); bc = convert(Vector{$T}, collect(backColor))
            GC.@preserve src dest tref bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{vImage_PerspectiveTransform}, Int32, Ptr{$T}, vImage_Flags),
                    sb, db, C_NULL, tref, Int32(interpolation), pointer(bc), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
    @eval $alloc(src::AbstractArray, transform::vImage_PerspectiveTransform; kw...) =
        $bang(_alloc_like(src), src, transform; kw...)
end

"""
    perspectiveWarp_Planar8(src, transform::vImage_PerspectiveTransform; interpolation=1, backColor, flags)

Apply a projective (perspective) warp. `interpolation` selects nearest (`0`) or linear
(`1`). Formats `Planar8/16U/16F`, `ARGB8888/16U/16F`.
"""
perspectiveWarp_Planar8

# ResamplingFilter (opaque, caller-owned) used by the shear functions.
function _new_resampling_filter(scale::Real, flags::Integer)
    f = ccall((:vImageNewResamplingFilter, vimage_lib), Ptr{Cvoid},
              (Cfloat, vImage_Flags), Cfloat(scale), vImage_Flags(flags))
    f == C_NULL && throw(vImageError(kvImageMemoryAllocationError, "vImageNewResamplingFilter"))
    return f
end
_destroy_resampling_filter(f::Ptr{Cvoid}) =
    ccall((:vImageDestroyResamplingFilter, vimage_lib), Cvoid, (Ptr{Cvoid},), f)

# ---- Horizontal / Vertical shear (float + double) -----------------------------------
# (src, dest, offX, offY, translate, shearSlope, ResamplingFilter, backColor, flags)
# The shear functions require a non-NULL ResamplingFilter; we create one at scale 1.0
# (pure shear, no resize), use it, and destroy it.
for (op, CT) in ((:HorizontalShear, Cfloat), (:VerticalShear, Cfloat),
                 (:HorizontalShearD, Cdouble), (:VerticalShearD, Cdouble))
    fmts = CT === Cfloat ?
        (:Planar8, :Planar16U, :Planar16S, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S,
         :ARGBFFFF, :Planar16F, :CbCr16F, :CbCr16U, :CbCr16S, :CbCr8, :ARGB16F) :
        (:Planar8, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S, :ARGBFFFF,
         :Planar16F, :CbCr16F, :CbCr16U, :CbCr16S, :ARGB16F)
    for sfx in fmts
        T, nch = _VFMT[sfx]
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        if nch == 1
            @eval function $bang(dest::AbstractArray, src::AbstractArray, translate::Real, shearSlope::Real;
                                 srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
                                 backColor = zero($T), flags::Integer = kvImageBackgroundColorFill)
                filt = _new_resampling_filter(1.0, flags)
                try
                    GC.@preserve src dest begin
                        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                        err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Csize_t, Csize_t, $CT, $CT, Ptr{Cvoid}, $T, vImage_Flags),
                            sb, db, Csize_t(srcOffsetX), Csize_t(srcOffsetY), $CT(translate), $CT(shearSlope),
                            filt, $T(backColor), vImage_Flags(flags))
                    end
                    _check(err, $(String(sym)))
                finally
                    _destroy_resampling_filter(filt)
                end
                return dest
            end
        else
            @eval function $bang(dest::AbstractArray, src::AbstractArray, translate::Real, shearSlope::Real;
                                 srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
                                 backColor = zeros($T, $nch), flags::Integer = kvImageBackgroundColorFill)
                bc = convert(Vector{$T}, collect(backColor))
                filt = _new_resampling_filter(1.0, flags)
                try
                    GC.@preserve src dest bc begin
                        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                        err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Csize_t, Csize_t, $CT, $CT, Ptr{Cvoid}, Ptr{$T}, vImage_Flags),
                            sb, db, Csize_t(srcOffsetX), Csize_t(srcOffsetY), $CT(translate), $CT(shearSlope),
                            filt, pointer(bc), vImage_Flags(flags))
                    end
                    _check(err, $(String(sym)))
                finally
                    _destroy_resampling_filter(filt)
                end
                return dest
            end
        end
        @eval $alloc(src::AbstractArray, translate::Real, shearSlope::Real; kw...) =
            $bang(_alloc_like(src), src, translate, shearSlope; kw...)
    end
end

"""
    horizontalShear_Planar8(src, xTranslate, shearSlope; backColor, flags)
    verticalShearD_ARGB8888(src, yTranslate, shearSlope; ...)

Shear an image using the default (Lanczos) resampling filter. `*Shear*` take a
single-precision translate/slope; the `*ShearD*` forms take double precision.
Allocating (same-size) and `!` mutating variants.
"""
horizontalShear_Planar8

# =====================================================================================
# MORPHOLOGY
# =====================================================================================

# ---- Dilate / Erode : (src, dest, offX, offY, kernel*, kh, kw, flags) ---------------
for op in (:Dilate, :Erode)
    for sfx in (:Planar8, :PlanarF, :ARGB8888, :ARGBFFFF)
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray, kernel::AbstractMatrix{UInt8};
                           srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
                           flags::Integer = kvImageNoFlags)
                k = kernel isa Matrix{UInt8} ? kernel : convert(Matrix{UInt8}, kernel)
                kh, kw = size(k)  # kernel stored (height, width) row-major; symmetric use is typical
                GC.@preserve src dest k begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Csize_t, Csize_t,
                         Ptr{UInt8}, Csize_t, Csize_t, vImage_Flags),
                        sb, db, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                        pointer(k), Csize_t(kh), Csize_t(kw), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractArray, kernel::AbstractMatrix{UInt8}; kw...) =
                $bang(_alloc_like(src), src, kernel; kw...)
        end
    end
end

"""
    dilate_Planar8(src, kernel::Matrix{UInt8}; flags) -> array
    erode_ARGB8888(src, kernel::Matrix{UInt8}; flags) -> array

Grayscale morphological dilation / erosion with an arbitrary `kernel` (structuring
element). `kernel` size is `(kernel_height, kernel_width)`. Allocating and `!` variants.
Formats `Planar8/F`, `ARGB8888/FFFF`.
"""
dilate_Planar8

# ---- Max / Min (rectangular) : (src, dest, temp, offX, offY, kh, kw, flags) ----------
for op in (:Max, :Min)
    for sfx in (:Planar8, :PlanarF, :ARGB8888, :ARGBFFFF)
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray, kernelHeight::Integer, kernelWidth::Integer;
                           srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
                           flags::Integer = kvImageNoFlags)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                         Csize_t, Csize_t, vImage_Flags),
                        sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                        Csize_t(kernelHeight), Csize_t(kernelWidth), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractArray, kernelHeight::Integer, kernelWidth::Integer; kw...) =
                $bang(_alloc_like(src), src, kernelHeight, kernelWidth; kw...)
        end
    end
end

"""
    max_Planar8(src, kernelHeight, kernelWidth; flags) -> array
    min_ARGB8888(src, kernelHeight, kernelWidth; flags) -> array

Rectangular maximum / minimum filter (dilation / erosion by a solid
`kernelHeight × kernelWidth` rectangle). Allocating and `!` variants; formats
`Planar8/F`, `ARGB8888/FFFF`.
"""
max_Planar8

# =====================================================================================
# CONVOLUTION
# =====================================================================================

# Formats using an integer (int16) kernel with a divisor:
const _CONV_INT = (:Planar8, :ARGB8888)
# Formats using a float kernel (no divisor):
const _CONV_FLT = (:PlanarF, :Planar16F, :ARGBFFFF, :ARGB16F)

# ---- Convolve -----------------------------------------------------------------------
for sfx in _CONV_INT
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageConvolve_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    bcarg = nch == 1 ? :($T(backColor)) : :(pointer(bc))
    bctyp = nch == 1 ? T : Ptr{T}
    prep  = nch == 1 ? :(bc = $T[]) : :(bc = convert(Vector{$T}, collect(backColor)))
    @eval function $bang(dest::AbstractArray, src::AbstractArray, kernel::AbstractMatrix{<:Integer};
                        divisor::Integer = 1, backColor = $(nch == 1 ? :(zero($T)) : :(zeros($T, $nch))),
                        srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageNoFlags)
        k = convert(Matrix{Int16}, kernel); kh, kw = size(k)
        $prep
        GC.@preserve src dest k bc begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                 Ptr{Int16}, UInt32, UInt32, Int32, $bctyp, vImage_Flags),
                sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                pointer(k), UInt32(kh), UInt32(kw), Int32(divisor), $bcarg, vImage_Flags(flags))
        end
        _check(err, $(String(sym)))
        return dest
    end
    @eval $alloc(src::AbstractArray, kernel::AbstractMatrix{<:Integer}; kw...) =
        $bang(_alloc_like(src), src, kernel; kw...)
end
for sfx in _CONV_FLT
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageConvolve_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    bcarg = nch == 1 ? :($T(backColor)) : :(pointer(bc))
    bctyp = nch == 1 ? T : Ptr{T}
    prep  = nch == 1 ? :(bc = $T[]) : :(bc = convert(Vector{$T}, collect(backColor)))
    @eval function $bang(dest::AbstractArray, src::AbstractArray, kernel::AbstractMatrix{<:Real};
                        backColor = $(nch == 1 ? :(zero($T)) : :(zeros($T, $nch))),
                        srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageNoFlags)
        k = convert(Matrix{Cfloat}, kernel); kh, kw = size(k)
        $prep
        GC.@preserve src dest k bc begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                 Ptr{Cfloat}, UInt32, UInt32, $bctyp, vImage_Flags),
                sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                pointer(k), UInt32(kh), UInt32(kw), $bcarg, vImage_Flags(flags))
        end
        _check(err, $(String(sym)))
        return dest
    end
    @eval $alloc(src::AbstractArray, kernel::AbstractMatrix{<:Real}; kw...) =
        $bang(_alloc_like(src), src, kernel; kw...)
end

"""
    convolve_Planar8(src, kernel::Matrix{<:Integer}; divisor=1, backColor, flags) -> array
    convolve_ARGBFFFF(src, kernel::Matrix{<:Real}; backColor, flags) -> array

General 2-D convolution. Integer-pixel formats (`Planar8`, `ARGB8888`) take an
integer `kernel` and an integer `divisor` (output = Σ(kernel .* window) / divisor);
floating-point formats (`PlanarF`, `Planar16F`, `ARGBFFFF`, `ARGB16F`) take a real
`kernel`. `kernel` is `(kernel_height, kernel_width)`. Allocating and `!` variants.
"""
convolve_Planar8

# ---- ConvolveWithBias ---------------------------------------------------------------
for sfx in _CONV_INT
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageConvolveWithBias_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    bcarg = nch == 1 ? :($T(backColor)) : :(pointer(bc))
    bctyp = nch == 1 ? T : Ptr{T}
    prep  = nch == 1 ? :(bc = $T[]) : :(bc = convert(Vector{$T}, collect(backColor)))
    @eval function $bang(dest::AbstractArray, src::AbstractArray, kernel::AbstractMatrix{<:Integer};
                        divisor::Integer = 1, bias::Integer = 0,
                        backColor = $(nch == 1 ? :(zero($T)) : :(zeros($T, $nch))),
                        srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageNoFlags)
        k = convert(Matrix{Int16}, kernel); kh, kw = size(k)
        $prep
        GC.@preserve src dest k bc begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                 Ptr{Int16}, UInt32, UInt32, Int32, Int32, $bctyp, vImage_Flags),
                sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                pointer(k), UInt32(kh), UInt32(kw), Int32(divisor), Int32(bias), $bcarg, vImage_Flags(flags))
        end
        _check(err, $(String(sym)))
        return dest
    end
    @eval $alloc(src::AbstractArray, kernel::AbstractMatrix{<:Integer}; kw...) =
        $bang(_alloc_like(src), src, kernel; kw...)
end
for sfx in _CONV_FLT
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageConvolveWithBias_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    bcarg = nch == 1 ? :($T(backColor)) : :(pointer(bc))
    bctyp = nch == 1 ? T : Ptr{T}
    prep  = nch == 1 ? :(bc = $T[]) : :(bc = convert(Vector{$T}, collect(backColor)))
    @eval function $bang(dest::AbstractArray, src::AbstractArray, kernel::AbstractMatrix{<:Real};
                        bias::Real = 0, backColor = $(nch == 1 ? :(zero($T)) : :(zeros($T, $nch))),
                        srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageNoFlags)
        k = convert(Matrix{Cfloat}, kernel); kh, kw = size(k)
        $prep
        GC.@preserve src dest k bc begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                 Ptr{Cfloat}, UInt32, UInt32, Cfloat, $bctyp, vImage_Flags),
                sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                pointer(k), UInt32(kh), UInt32(kw), Cfloat(bias), $bcarg, vImage_Flags(flags))
        end
        _check(err, $(String(sym)))
        return dest
    end
    @eval $alloc(src::AbstractArray, kernel::AbstractMatrix{<:Real}; kw...) =
        $bang(_alloc_like(src), src, kernel; kw...)
end

"""
    convolveWithBias_Planar8(src, kernel; divisor=1, bias=0, backColor, flags) -> array

Like [`convolve_Planar8`](@ref) but adds `bias` to the (divided) accumulator before
storing. Allocating and `!` variants for `Planar8`, `ARGB8888`, `PlanarF`, `Planar16F`,
`ARGBFFFF`, `ARGB16F`.
"""
convolveWithBias_Planar8

# ---- Box / Tent convolve : (src, dest, temp, offX, offY, kh, kw, backColor, flags) ---
for op in (:BoxConvolve, :TentConvolve)
    for sfx in (:Planar8, :ARGB8888)
        T, nch = _VFMT[sfx]
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        if nch == 1
            @eval function $bang(dest::AbstractArray, src::AbstractArray, kernelHeight::Integer, kernelWidth::Integer;
                                 backColor = zero($T), srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
                                 flags::Integer = kvImageEdgeExtend)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                         UInt32, UInt32, $T, vImage_Flags),
                        sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                        UInt32(kernelHeight), UInt32(kernelWidth), $T(backColor), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
        else
            @eval function $bang(dest::AbstractArray, src::AbstractArray, kernelHeight::Integer, kernelWidth::Integer;
                                 backColor = zeros($T, $nch), srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
                                 flags::Integer = kvImageEdgeExtend)
                bc = convert(Vector{$T}, collect(backColor))
                GC.@preserve src dest bc begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                         UInt32, UInt32, Ptr{$T}, vImage_Flags),
                        sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                        UInt32(kernelHeight), UInt32(kernelWidth), pointer(bc), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
        end
        @eval $alloc(src::AbstractArray, kernelHeight::Integer, kernelWidth::Integer; kw...) =
            $bang(_alloc_like(src), src, kernelHeight, kernelWidth; kw...)
    end
end

"""
    boxConvolve_Planar8(src, kernelHeight, kernelWidth; backColor, flags) -> array
    tentConvolve_ARGB8888(src, kernelHeight, kernelWidth; backColor, flags) -> array

Fast box (uniform average) and tent (triangular/linear) blur over a
`kernelHeight × kernelWidth` window. Convolving a constant image reproduces the
constant. Allocating and `!` variants for `Planar8`, `ARGB8888`.
"""
boxConvolve_Planar8

# =====================================================================================
# HISTOGRAM
# =====================================================================================

"""
    histogramCalculation_Planar8(src) -> Vector{Int}   # length 256

Compute the 256-bin histogram of an 8-bit planar image.
"""
function histogramCalculation_Planar8(src::AbstractMatrix{UInt8})
    hist = zeros(Csize_t, 256)
    GC.@preserve src hist begin
        sb = Ref(vimage_buffer(src))
        err = ccall((:vImageHistogramCalculation_Planar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{Csize_t}, vImage_Flags), sb, pointer(hist), kvImageNoFlags)
    end
    _check(err, "vImageHistogramCalculation_Planar8")
    return Int.(hist)
end

"""
    histogramCalculation_PlanarF(src, entries, minVal, maxVal) -> Vector{Int}

Compute an `entries`-bin histogram of a float planar image over `[minVal, maxVal]`.
"""
function histogramCalculation_PlanarF(src::AbstractMatrix{Float32}, entries::Integer,
                                      minVal::Real, maxVal::Real)
    hist = zeros(Csize_t, entries)
    GC.@preserve src hist begin
        sb = Ref(vimage_buffer(src))
        err = ccall((:vImageHistogramCalculation_PlanarF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{Csize_t}, Cuint, Cfloat, Cfloat, vImage_Flags),
            sb, pointer(hist), Cuint(entries), Cfloat(minVal), Cfloat(maxVal), kvImageNoFlags)
    end
    _check(err, "vImageHistogramCalculation_PlanarF")
    return Int.(hist)
end

"""
    histogramCalculation_ARGB8888(src) -> Matrix{Int}   # 256 × 4 (columns A,R,G,B)

Compute per-channel 256-bin histograms of an 8-bit interleaved 4-channel image.
"""
function histogramCalculation_ARGB8888(src::AbstractArray{UInt8,3})
    chans = [zeros(Csize_t, 256) for _ in 1:4]
    GC.@preserve src chans begin
        ptrs = [pointer(c) for c in chans]
        GC.@preserve ptrs begin
            sb = Ref(vimage_buffer(src))
            err = ccall((:vImageHistogramCalculation_ARGB8888, vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{Ptr{Csize_t}}, vImage_Flags),
                sb, pointer(ptrs), kvImageNoFlags)
        end
    end
    _check(err, "vImageHistogramCalculation_ARGB8888")
    return hcat((Int.(c) for c in chans)...)
end

"""
    histogramCalculation_ARGBFFFF(src, entries, minVal, maxVal) -> Matrix{Int}  # entries × 4

Per-channel histograms of a float interleaved 4-channel image over `[minVal, maxVal]`.
"""
function histogramCalculation_ARGBFFFF(src::AbstractArray{Float32,3}, entries::Integer,
                                       minVal::Real, maxVal::Real)
    chans = [zeros(Csize_t, entries) for _ in 1:4]
    GC.@preserve src chans begin
        ptrs = [pointer(c) for c in chans]
        GC.@preserve ptrs begin
            sb = Ref(vimage_buffer(src))
            err = ccall((:vImageHistogramCalculation_ARGBFFFF, vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{Ptr{Csize_t}}, Cuint, Cfloat, Cfloat, vImage_Flags),
                sb, pointer(ptrs), Cuint(entries), Cfloat(minVal), Cfloat(maxVal), kvImageNoFlags)
        end
    end
    _check(err, "vImageHistogramCalculation_ARGBFFFF")
    return hcat((Int.(c) for c in chans)...)
end

# Equalization / ContrastStretch : simple integer formats are (src, dest, flags).
for op in (:Equalization, :ContrastStretch)
    for sfx in (:Planar8, :ARGB8888)
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray; flags::Integer = kvImageNoFlags)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                        sb, db, vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractArray; flags::Integer = kvImageNoFlags) =
                $bang(_alloc_like(src), src; flags = flags)
        end
    end
    # Float formats: (src, dest, tempBuffer, entries, minVal, maxVal, flags)
    for sfx in (:PlanarF, :ARGBFFFF)
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray; entries::Integer = 4096,
                           minVal::Real = 0f0, maxVal::Real = 1f0, flags::Integer = kvImageNoFlags)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Cuint, Cfloat, Cfloat, vImage_Flags),
                        sb, db, C_NULL, Cuint(entries), Cfloat(minVal), Cfloat(maxVal), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractArray; kw...) = $bang(_alloc_like(src), src; kw...)
        end
    end
end

"""
    equalization_Planar8(src; flags) -> array
    contrastStretch_ARGBFFFF(src; entries=4096, minVal=0, maxVal=1, flags) -> array

Histogram equalization / linear contrast stretch. Integer formats (`Planar8`,
`ARGB8888`) take no extra arguments; float formats (`PlanarF`, `ARGBFFFF`) take
histogram `entries` and the value range `[minVal, maxVal]`. Allocating and `!` variants.
"""
equalization_Planar8

# EndsInContrastStretch
function endsInContrastStretch_Planar8!(dest::AbstractMatrix{UInt8}, src::AbstractMatrix{UInt8};
                                        percentLow::Integer = 0, percentHigh::Integer = 0,
                                        flags::Integer = kvImageNoFlags)
    GC.@preserve src dest begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageEndsInContrastStretch_Planar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Cuint, Cuint, vImage_Flags),
            sb, db, Cuint(percentLow), Cuint(percentHigh), vImage_Flags(flags))
    end
    _check(err, "vImageEndsInContrastStretch_Planar8")
    return dest
end
"""
    endsInContrastStretch_Planar8(src; percentLow=0, percentHigh=0, flags) -> array

Ends-in contrast stretch of an 8-bit planar image: clip `percentLow`% of the darkest
and `percentHigh`% of the brightest pixels, then rescale. Allocating and `!` variants.
"""
endsInContrastStretch_Planar8(src::AbstractMatrix{UInt8}; kw...) =
    endsInContrastStretch_Planar8!(similar(src), src; kw...)

function endsInContrastStretch_ARGB8888!(dest::AbstractArray{UInt8,3}, src::AbstractArray{UInt8,3};
                                         percentLow = (0, 0, 0, 0), percentHigh = (0, 0, 0, 0),
                                         flags::Integer = kvImageNoFlags)
    pl = convert(Vector{Cuint}, collect(percentLow)); ph = convert(Vector{Cuint}, collect(percentHigh))
    GC.@preserve src dest pl ph begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageEndsInContrastStretch_ARGB8888, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cuint}, Ptr{Cuint}, vImage_Flags),
            sb, db, pointer(pl), pointer(ph), vImage_Flags(flags))
    end
    _check(err, "vImageEndsInContrastStretch_ARGB8888")
    return dest
end
"""
    endsInContrastStretch_ARGB8888(src; percentLow, percentHigh, flags) -> array

Per-channel ends-in contrast stretch of an ARGB8888 image; `percentLow`/`percentHigh`
are length-4 collections. Allocating and `!` variants.
"""
endsInContrastStretch_ARGB8888(src::AbstractArray{UInt8,3}; kw...) =
    endsInContrastStretch_ARGB8888!(similar(src), src; kw...)

# HistogramSpecification (integer formats)
function histogramSpecification_Planar8!(dest::AbstractMatrix{UInt8}, src::AbstractMatrix{UInt8},
                                         desiredHistogram::AbstractVector; flags::Integer = kvImageNoFlags)
    dh = convert(Vector{Csize_t}, collect(desiredHistogram))
    length(dh) == 256 || throw(DimensionMismatch("desiredHistogram must have 256 entries"))
    GC.@preserve src dest dh begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageHistogramSpecification_Planar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Csize_t}, vImage_Flags),
            sb, db, pointer(dh), vImage_Flags(flags))
    end
    _check(err, "vImageHistogramSpecification_Planar8")
    return dest
end
"""
    histogramSpecification_Planar8(src, desiredHistogram; flags) -> array

Remap an 8-bit planar image so its histogram matches `desiredHistogram` (256 entries).
Allocating and `!` variants.
"""
histogramSpecification_Planar8(src::AbstractMatrix{UInt8}, dh::AbstractVector; kw...) =
    histogramSpecification_Planar8!(similar(src), src, dh; kw...)

# =====================================================================================
# ALPHA / COMPOSITING
# =====================================================================================

# Two-buffer premultiply/unpremultiply/clip-to-alpha for interleaved ARGB/RGBA formats
# with the shape (src, dest, flags).
for op in (:PremultiplyData, :UnpremultiplyData, :ClipToAlpha)
    fmts = op === :ClipToAlpha ? (:ARGB8888, :ARGBFFFF, :RGBA8888, :RGBAFFFF) :
        (:ARGB8888, :ARGBFFFF, :RGBA8888, :RGBAFFFF, :RGBA16F, :ARGB16U, :RGBA16U, :ARGB16Q12, :RGBA16Q12)
    for sfx in fmts
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, src::AbstractArray; flags::Integer = kvImageNoFlags)
                GC.@preserve src dest begin
                    sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                        sb, db, vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractArray; flags::Integer = kvImageNoFlags) =
                $bang(_alloc_like(src), src; flags = flags)
        end
    end
end

# Planar premultiply/unpremultiply/clip-to-alpha: (src, alpha, dest, flags)
for op in (:PremultiplyData, :UnpremultiplyData, :ClipToAlpha)
    for sfx in (:Planar8, :PlanarF)
        sym   = Symbol("vImage", op, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractMatrix, src::AbstractMatrix, alpha::AbstractMatrix;
                           flags::Integer = kvImageNoFlags)
                GC.@preserve src alpha dest begin
                    sb = Ref(vimage_buffer(src)); ab = Ref(vimage_buffer(alpha)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                        sb, ab, db, vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(src::AbstractMatrix, alpha::AbstractMatrix; flags::Integer = kvImageNoFlags) =
                $bang(similar(src), src, alpha; flags = flags)
        end
    end
end

"""
    premultiplyData_ARGB8888(src; flags) -> array
    unpremultiplyData_ARGBFFFF(src; flags) -> array
    premultiplyData_Planar8(src, alpha; flags) -> array

Premultiply / un-premultiply colour channels by alpha. Interleaved formats
(`ARGB8888/FFFF/16U/16Q12`, `RGBA8888/FFFF/16F/16U/16Q12`) carry alpha inside the
pixel; planar formats (`Planar8`, `PlanarF`) take a separate `alpha` plane.
`premultiply` then `unpremultiply` round-trips (within rounding). Allocating and `!`
variants. See also [`clipToAlpha_ARGB8888`](@ref).
"""
premultiplyData_ARGB8888

"""
    clipToAlpha_ARGB8888(src; flags) -> array

Clamp each colour channel to not exceed its alpha (produces valid premultiplied data).
Allocating and `!` variants; interleaved `ARGB8888/FFFF`, `RGBA8888/FFFF`, and planar
`Planar8/F` (which take a separate `alpha` plane).
"""
clipToAlpha_ARGB8888

# AlphaBlend / PremultipliedAlphaBlend / NonpremultipliedToPremultiplied for interleaved
# formats: (srcTop, srcBottom, dest, flags)
for (op, cname) in ((:PremultipliedAlphaBlend, :PremultipliedAlphaBlend),
                    (:NonpremultipliedToPremultiplied, :AlphaBlend_NonpremultipliedToPremultiplied),
                    (:PlainAlphaBlend, :AlphaBlend))
    fmts = op === :PremultipliedAlphaBlend ?
        (:ARGB8888, :BGRA8888, :ARGBFFFF, :BGRAFFFF) :
        (op === :PlainAlphaBlend ? (:ARGB8888, :ARGBFFFF) : (:ARGB8888, :ARGBFFFF))
    for sfx in fmts
        sym   = Symbol("vImage", cname, "_", sfx)
        bang  = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, srcTop::AbstractArray, srcBottom::AbstractArray;
                           flags::Integer = kvImageNoFlags)
                GC.@preserve srcTop srcBottom dest begin
                    tb = Ref(vimage_buffer(srcTop)); bb = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                        tb, bb, db, vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(srcTop::AbstractArray, srcBottom::AbstractArray; flags::Integer = kvImageNoFlags) =
                $bang(_alloc_like(srcTop), srcTop, srcBottom; flags = flags)
        end
    end
end

# Named premultiplied blend modes (RGBA8888): (srcTop, srcBottom, dest, flags)
for op in (:PremultipliedAlphaBlendMultiply, :PremultipliedAlphaBlendScreen,
           :PremultipliedAlphaBlendDarken, :PremultipliedAlphaBlendLighten)
    sym   = Symbol("vImage", op, "_RGBA8888")
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray{UInt8,3}, srcTop::AbstractArray{UInt8,3},
                       srcBottom::AbstractArray{UInt8,3}; flags::Integer = kvImageNoFlags)
            GC.@preserve srcTop srcBottom dest begin
                tb = Ref(vimage_buffer(srcTop)); bb = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    tb, bb, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(srcTop::AbstractArray{UInt8,3}, srcBottom::AbstractArray{UInt8,3}; flags::Integer = kvImageNoFlags) =
            $bang(_alloc_like(srcTop), srcTop, srcBottom; flags = flags)
    end
end

"""
    premultipliedAlphaBlend_ARGB8888(srcTop, srcBottom; flags) -> array

Composite premultiplied `srcTop` over `srcBottom` (source-over). Named Porter–Duff /
blend-mode variants are also available: `premultipliedAlphaBlendMultiply_RGBA8888`,
`…Screen…`, `…Darken…`, `…Lighten…`, plus `alphaBlend_ARGB8888` (non-premultiplied
over) and `alphaBlend_NonpremultipliedToPremultiplied_ARGB8888`. Allocating and `!`
variants.
"""
premultipliedAlphaBlend_ARGB8888

# Planar alpha blends: (srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha, dest, flags)
function alphaBlend_Planar8!(dest::AbstractMatrix{UInt8}, srcTop::AbstractMatrix{UInt8}, srcTopAlpha::AbstractMatrix{UInt8},
                             srcBottom::AbstractMatrix{UInt8}, srcBottomAlpha::AbstractMatrix{UInt8},
                             alpha::AbstractMatrix{UInt8}; flags::Integer = kvImageNoFlags)
    GC.@preserve srcTop srcTopAlpha srcBottom srcBottomAlpha alpha dest begin
        t = Ref(vimage_buffer(srcTop)); ta = Ref(vimage_buffer(srcTopAlpha))
        b = Ref(vimage_buffer(srcBottom)); ba = Ref(vimage_buffer(srcBottomAlpha))
        al = Ref(vimage_buffer(alpha)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageAlphaBlend_Planar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer},
             Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            t, ta, b, ba, al, db, vImage_Flags(flags))
    end
    _check(err, "vImageAlphaBlend_Planar8")
    return dest
end
"""
    alphaBlend_Planar8(srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha; flags) -> array

Non-premultiplied planar alpha blend of two images given their alpha planes and an
overall blend `alpha` plane. Allocating and `!` variants.
"""
alphaBlend_Planar8(srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha; kw...) =
    alphaBlend_Planar8!(similar(srcTop), srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha; kw...)

# =====================================================================================
# TRANSFORM
# =====================================================================================

# MatrixMultiply for interleaved ARGB (4x4 matrix). Integer & float.
function matrixMultiply_ARGB8888!(dest::AbstractArray{UInt8,3}, src::AbstractArray{UInt8,3},
                                  matrix::AbstractMatrix{<:Integer}; divisor::Integer = 256,
                                  preBias = nothing, postBias = nothing, flags::Integer = kvImageNoFlags)
    m = convert(Matrix{Int16}, matrix)  # row-major 4x4 as flat
    length(m) == 16 || throw(DimensionMismatch("matrix must be 4×4"))
    mflat = vec(permutedims(m))  # vImage expects row-major flattening
    prev = preBias === nothing ? Int16[] : convert(Vector{Int16}, collect(preBias))
    postv = postBias === nothing ? Int32[] : convert(Vector{Int32}, collect(postBias))
    GC.@preserve src dest mflat prev postv begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        prep = isempty(prev) ? Ptr{Int16}(C_NULL) : pointer(prev)
        postp = isempty(postv) ? Ptr{Int32}(C_NULL) : pointer(postv)
        err = ccall((:vImageMatrixMultiply_ARGB8888, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Int16}, Int32, Ptr{Int16}, Ptr{Int32}, vImage_Flags),
            sb, db, pointer(mflat), Int32(divisor), prep, postp, vImage_Flags(flags))
    end
    _check(err, "vImageMatrixMultiply_ARGB8888")
    return dest
end
"""
    matrixMultiply_ARGB8888(src, matrix::Matrix{<:Integer}; divisor=256, preBias, postBias, flags) -> array

Apply a 4×4 colour `matrix` to every ARGB8888 pixel: `out = (matrix*(in+preBias))/divisor + postBias`.
`matrix` is given in natural (row-per-output-channel) order. Allocating and `!` variants.
"""
matrixMultiply_ARGB8888(src::AbstractArray{UInt8,3}, matrix::AbstractMatrix{<:Integer}; kw...) =
    matrixMultiply_ARGB8888!(similar(src), src, matrix; kw...)

function matrixMultiply_ARGBFFFF!(dest::AbstractArray{Float32,3}, src::AbstractArray{Float32,3},
                                  matrix::AbstractMatrix{<:Real}; preBias = nothing, postBias = nothing,
                                  flags::Integer = kvImageNoFlags)
    m = convert(Matrix{Cfloat}, matrix)
    length(m) == 16 || throw(DimensionMismatch("matrix must be 4×4"))
    mflat = vec(permutedims(m))
    prev = preBias === nothing ? Cfloat[] : convert(Vector{Cfloat}, collect(preBias))
    postv = postBias === nothing ? Cfloat[] : convert(Vector{Cfloat}, collect(postBias))
    GC.@preserve src dest mflat prev postv begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        prep = isempty(prev) ? Ptr{Cfloat}(C_NULL) : pointer(prev)
        postp = isempty(postv) ? Ptr{Cfloat}(C_NULL) : pointer(postv)
        err = ccall((:vImageMatrixMultiply_ARGBFFFF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, vImage_Flags),
            sb, db, pointer(mflat), prep, postp, vImage_Flags(flags))
    end
    _check(err, "vImageMatrixMultiply_ARGBFFFF")
    return dest
end
"""
    matrixMultiply_ARGBFFFF(src, matrix::Matrix{<:Real}; preBias, postBias, flags) -> array

Floating-point 4×4 colour-matrix multiply of an ARGBFFFF image. Allocating and `!` variants.
"""
matrixMultiply_ARGBFFFF(src::AbstractArray{Float32,3}, matrix::AbstractMatrix{<:Real}; kw...) =
    matrixMultiply_ARGBFFFF!(similar(src), src, matrix; kw...)

# PiecewiseGamma (planar, single-channel; boundary type varies)
for (sfx, T, BT) in ((:Planar8, UInt8, UInt8), (:PlanarF, Float32, Cfloat),
                     (:Planar16Q12, Int16, Int16))
    sym   = Symbol("vImagePiecewiseGamma_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractMatrix{$T}, src::AbstractMatrix{$T};
                       exponentialCoeffs = (1f0, 0f0, 0f0), gamma::Real = 1f0,
                       linearCoeffs = (1f0, 0f0), boundary = zero($BT), flags::Integer = kvImageNoFlags)
            ec = convert(Vector{Cfloat}, collect(exponentialCoeffs)); lc = convert(Vector{Cfloat}, collect(linearCoeffs))
            GC.@preserve src dest ec lc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cfloat}, Cfloat, Ptr{Cfloat}, $BT, vImage_Flags),
                    sb, db, pointer(ec), Cfloat(gamma), pointer(lc), $BT(boundary), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(src::AbstractMatrix{$T}; kw...) = $bang(similar(src), src; kw...)
    end
end

"""
    piecewiseGamma_PlanarF(src; exponentialCoeffs, gamma, linearCoeffs, boundary, flags) -> array

Apply a piecewise gamma curve: for `x ≥ boundary`, `out = (a·x + b)^gamma + c` with
`exponentialCoeffs = (a, b, c)`; below the boundary a linear segment
`linearCoeffs = (d, e)` (`out = d·x + e`) is used. Allocating and `!` variants for
`Planar8`, `PlanarF`, `Planar16Q12`.
"""
piecewiseGamma_PlanarF

# LookupTable Planar8 -> various (table has 256 entries)
for (sfx, DT, TBL) in ((:Planar8toPlanar16, UInt16, UInt16), (:Planar8toPlanarF, Float32, Cfloat))
    sym   = Symbol("vImageLookupTable_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractMatrix{$DT}, src::AbstractMatrix{UInt8}, table::AbstractVector;
                       flags::Integer = kvImageNoFlags)
            t = convert(Vector{$TBL}, collect(table))
            length(t) == 256 || throw(DimensionMismatch("table must have 256 entries"))
            GC.@preserve src dest t begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{$TBL}, vImage_Flags),
                    sb, db, pointer(t), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    lookupTable_Planar8toPlanarF(dest, src, table) -> dest

Map each 8-bit source pixel through a 256-entry `table`. Mutating variants for
`Planar8toPlanar16` (UInt16 table/output) and `Planar8toPlanarF` (Float32).
"""
lookupTable_Planar8toPlanarF!

# FloodFill (in place on a single buffer)
for (sfx, T, nch) in ((:Planar8, UInt8, 1), (:Planar16U, UInt16, 1),
                      (:ARGB8888, UInt8, 4), (:ARGB16U, UInt16, 4))
    sym   = Symbol("vImageFloodFill_", sfx)
    fn    = _jlname(sym) # mutating in place; keep base name (operates on srcDest)
    if nch == 1
        @eval function $(Symbol(fn, "!"))(srcDest::AbstractMatrix{$T}, seedX::Integer, seedY::Integer,
                                          newValue; connectivity::Integer = 4, flags::Integer = kvImageNoFlags)
            GC.@preserve srcDest begin
                sb = Ref(vimage_buffer(srcDest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t, $T, Cint, vImage_Flags),
                    sb, C_NULL, Csize_t(seedX), Csize_t(seedY), $T(newValue), Cint(connectivity), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return srcDest
        end
    else
        @eval function $(Symbol(fn, "!"))(srcDest::AbstractArray{$T,3}, seedX::Integer, seedY::Integer,
                                          newValue; connectivity::Integer = 4, flags::Integer = kvImageNoFlags)
            nv = convert(Vector{$T}, collect(newValue))
            GC.@preserve srcDest nv begin
                sb = Ref(vimage_buffer(srcDest))
                # Pixel_8888/Pixel_ARGB_16U passed by value as a 4-element array -> pointer
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t, Ptr{$T}, Cint, vImage_Flags),
                    sb, C_NULL, Csize_t(seedX), Csize_t(seedY), pointer(nv), Cint(connectivity), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return srcDest
        end
    end
end
"""
    floodFill_Planar8!(srcDest, seedX, seedY, newValue; connectivity=4, flags)

Flood-fill (in place) the connected region of `srcDest` containing pixel
`(seedX, seedY)` (0-based, vImage coordinates) with `newValue`. `connectivity` is 4 or
8. Variants for `Planar8`, `Planar16U`, `ARGB8888`, `ARGB16U` (interleaved take a
length-4 `newValue`).
"""
floodFill_Planar8!

# =====================================================================================
# CONVERSION — channel manipulation
# =====================================================================================

# PermuteChannels : (src, dest, permuteMap, flags). ARGB (4) and RGB888 (3).
for (sfx, np) in ((:ARGB8888, 4), (:ARGB16U, 4), (:ARGBFFFF, 4), (:ARGB16F, 4), (:RGB888, 3))
    sym   = Symbol("vImagePermuteChannels_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray, src::AbstractArray, permuteMap; flags::Integer = kvImageNoFlags)
            pm = convert(Vector{UInt8}, collect(permuteMap))
            length(pm) == $np || throw(DimensionMismatch("permuteMap must have $($np) entries"))
            GC.@preserve src dest pm begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, vImage_Flags),
                    sb, db, pointer(pm), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(src::AbstractArray, permuteMap; flags::Integer = kvImageNoFlags) =
            $bang(_alloc_like(src), src, permuteMap; flags = flags)
    end
end
"""
    permuteChannels_ARGB8888(src, permuteMap; flags) -> array

Reorder the channels of an interleaved image. `permuteMap` (0-based, length 4 for ARGB
formats or 3 for `RGB888`) gives, for each destination channel, the source channel
index. Allocating and `!` variants for `ARGB8888/16U/16F/FFFF` and `RGB888`.
"""
permuteChannels_ARGB8888

# ExtractChannel : (src, dest, channelIndex, flags)
for sfx in (:ARGB8888, :ARGB16U, :ARGBFFFF)
    T, _ = _VFMT[sfx]
    sym   = Symbol("vImageExtractChannel_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractMatrix, src::AbstractArray{<:Any,3}, channelIndex::Integer;
                       flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Clong, vImage_Flags),
                    sb, db, Clong(channelIndex), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        function $alloc(src::AbstractArray{$T,3}, channelIndex::Integer; flags::Integer = kvImageNoFlags)
            dest = Matrix{$T}(undef, size(src, 2), size(src, 3))
            return $bang(dest, src, channelIndex; flags = flags)
        end
    end
end
"""
    extractChannel_ARGB8888(src, channelIndex; flags) -> Matrix

Extract one channel (0-based `channelIndex`) of an interleaved image into a planar
buffer. Allocating and `!` variants for `ARGB8888`, `ARGB16U`, `ARGBFFFF`.
"""
extractChannel_ARGB8888

# BufferFill : (dest, color, flags) — color decays to pointer for these formats.
for sfx in (:ARGB8888, :ARGB16U, :ARGB16S, :ARGBFFFF, :CbCr8, :CbCr16U, :CbCr16S)
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageBufferFill_", sfx)
    fn    = _jlname(sym)
    @eval begin
        function $(Symbol(fn, "!"))(dest::AbstractArray{$T,3}, color; flags::Integer = kvImageNoFlags)
            c = convert(Vector{$T}, collect(color))
            length(c) == $nch || throw(DimensionMismatch("color must have $($nch) entries"))
            GC.@preserve dest c begin
                db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{$T}, vImage_Flags),
                    db, pointer(c), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    bufferFill_ARGB8888!(dest, color; flags) -> dest

Fill every pixel of interleaved `dest` with the constant `color` (a length-`channels`
collection). Variants for `ARGB8888/16U/16S/FFFF` and `CbCr8/16U/16S`.
"""
bufferFill_ARGB8888!

# OverwriteChannelsWithScalar (planar) : (scalar, dest, flags)
for sfx in (:Planar8, :PlanarF, :Planar16S, :Planar16U, :Planar16F)
    T, _ = _VFMT[sfx]
    ST = sfx === :Planar16F ? UInt16 : T
    sym   = Symbol("vImageOverwriteChannelsWithScalar_", sfx)
    fn    = _jlname(sym)
    @eval begin
        function $(Symbol(fn, "!"))(dest::AbstractMatrix{$T}, scalar; flags::Integer = kvImageNoFlags)
            GC.@preserve dest begin
                db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    ($ST, Ptr{vImage_Buffer}, vImage_Flags),
                    $ST(scalar), db, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    overwriteChannelsWithScalar_Planar8!(dest, scalar; flags) -> dest

Overwrite an entire planar buffer with a constant `scalar`. Variants for
`Planar8/F/16S/16U/16F`.
"""
overwriteChannelsWithScalar_Planar8!

# TableLookUp_Planar8 : (src, dest, table[256], flags)
function tableLookUp_Planar8!(dest::AbstractMatrix{UInt8}, src::AbstractMatrix{UInt8},
                              table::AbstractVector{UInt8}; flags::Integer = kvImageNoFlags)
    length(table) == 256 || throw(DimensionMismatch("table must have 256 entries"))
    t = table isa Vector{UInt8} ? table : convert(Vector{UInt8}, table)
    GC.@preserve src dest t begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageTableLookUp_Planar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, vImage_Flags),
            sb, db, pointer(t), vImage_Flags(flags))
    end
    _check(err, "vImageTableLookUp_Planar8")
    return dest
end
"""
    tableLookUp_Planar8(src, table) -> Matrix

Map each 8-bit pixel through a 256-entry `table`. Allocating and `!` variants.
"""
tableLookUp_Planar8(src::AbstractMatrix{UInt8}, table::AbstractVector{UInt8}; kw...) =
    tableLookUp_Planar8!(similar(src), src, table; kw...)

# =====================================================================================
# CONVERSION — format / bit-depth (mutating; caller supplies correctly-typed dest)
# =====================================================================================

# Shape (src, dest, flags): dest element type/shape differs from src, so only a
# mutating `foo!(dest, src; flags)` is generated. Symbol list:
const _CONVERT_SS = (
    :vImageConvert_16Uto16F, :vImageConvert_16Fto16U, :vImageConvert_12UTo16U, :vImageConvert_16UTo12U,
    :vImageConvert_ARGBFFFFtoRGBFFF, :vImageConvert_RGBAFFFFtoRGBFFF, :vImageConvert_BGRAFFFFtoRGBFFF,
    :vImageConvert_ARGB1555toARGB8888, :vImageConvert_ARGB8888toARGB1555,
    :vImageConvert_RGBA5551toRGBA8888, :vImageConvert_RGBA8888toRGBA5551,
    :vImageConvert_RGB565toRGB888, :vImageConvert_ARGB8888toRGB565,
    :vImageConvert_RGBA8888toRGB565, :vImageConvert_BGRA8888toRGB565,
    :vImageConvert_RGBA5551toRGB565, :vImageConvert_ARGB1555toRGB565,
    :vImageConvert_Planar16FtoPlanarF, :vImageConvert_PlanarFtoPlanar16F,
    :vImageConvert_Planar8toPlanar16F, :vImageConvert_Planar16FtoPlanar8,
    :vImageConvert_16UToPlanar8, :vImageConvert_Planar8To16U,
    :vImageConvert_ARGB8888toRGB888, :vImageConvert_RGBA8888toRGB888, :vImageConvert_BGRA8888toRGB888,
    :vImageConvert_ARGB16UtoRGB16U, :vImageConvert_RGBA16UtoRGB16U, :vImageConvert_BGRA16UtoRGB16U,
    :vImageConvert_Planar1toPlanar8, :vImageConvert_Planar2toPlanar8, :vImageConvert_Planar4toPlanar8,
    :vImageConvert_8to16Q12, :vImageConvert_16Q12to8, :vImageConvert_16Q12to16F, :vImageConvert_16Fto16Q12,
    :vImageConvert_16Q12toF, :vImageConvert_Fto16Q12, :vImageConvert_16Q12to16U, :vImageConvert_16Uto16Q12,
    :vImageByteSwap_Planar16U,
)
for sym in _CONVERT_SS
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractArray, src::AbstractArray; flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    sb, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    convert_16UToPlanar8!(dest, src; flags) -> dest
    convert_ARGB8888toRGB888!(dest, src; flags) -> dest

Format / bit-depth conversions with the shape `(src, dest, flags)`. Because the
destination pixel format differs from the source, only the mutating `foo!(dest, src)`
form is provided — allocate a correctly-typed `dest` and pass it in. See the vImage
docs page for the full list of wrapped conversions.
"""
convert_16UToPlanar8!

# Shape (src, dest, maxFloat, minFloat, flags): float<->int planar scaling + clip.
for sym in (:vImageConvert_Planar8toPlanarF, :vImageConvert_PlanarFtoPlanar8, :vImageClip_PlanarF)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractMatrix, src::AbstractMatrix; maxFloat::Real = 1f0,
                       minFloat::Real = 0f0, flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Cfloat, Cfloat, vImage_Flags),
                    sb, db, Cfloat(maxFloat), Cfloat(minFloat), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    convert_Planar8toPlanarF!(dest, src; maxFloat=1, minFloat=0, flags) -> dest
    convert_PlanarFtoPlanar8!(dest, src; maxFloat=1, minFloat=0, flags) -> dest
    clip_PlanarF!(dest, src; maxFloat=1, minFloat=0, flags) -> dest

Planar 8-bit ↔ float conversion mapping `[0,255] ↔ [minFloat, maxFloat]`, and float
clamping (`clip_PlanarF!`). Also provided allocating: `convert_Planar8toPlanarF(src)`
returns a `Matrix{Float32}`; `convert_PlanarFtoPlanar8(src)` returns `Matrix{UInt8}`.
"""
convert_Planar8toPlanarF!
convert_Planar8toPlanarF(src::AbstractMatrix{UInt8}; kw...) =
    convert_Planar8toPlanarF!(Matrix{Float32}(undef, size(src)...), src; kw...)
convert_PlanarFtoPlanar8(src::AbstractMatrix{Float32}; kw...) =
    convert_PlanarFtoPlanar8!(Matrix{UInt8}(undef, size(src)...), src; kw...)

# Shape (src, dest, offset, scale, flags)
for sym in (:vImageConvert_16SToF, :vImageConvert_16UToF, :vImageConvert_FTo16S, :vImageConvert_FTo16U)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractMatrix, src::AbstractMatrix; offset::Real = 0f0,
                       scale::Real = 1f0, flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Cfloat, Cfloat, vImage_Flags),
                    sb, db, Cfloat(offset), Cfloat(scale), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    convert_16UToF!(dest, src; offset=0, scale=1, flags) -> dest

Integer ↔ float planar conversions with an affine map `out = (in - offset)/scale` (or
its inverse). Mutating variants for `16SToF`, `16UToF`, `FTo16S`, `FTo16U`.
"""
convert_16UToF!

# CopyBuffer : (src, dest, pixelSize, flags)
function copyBuffer!(dest::AbstractArray, src::AbstractArray; pixelSize::Integer, flags::Integer = kvImageNoFlags)
    GC.@preserve src dest begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageCopyBuffer, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Csize_t, vImage_Flags),
            sb, db, Csize_t(pixelSize), vImage_Flags(flags))
    end
    _check(err, "vImageCopyBuffer")
    return dest
end
"""
    copyBuffer!(dest, src; pixelSize, flags) -> dest

Copy `src` to `dest` accounting for `rowBytes` padding. `pixelSize` is the size of one
pixel in bytes.
"""
copyBuffer!

# Flatten : (src, dest, backgroundColor, isImagePremultiplied, flags) — RGB output.
for (sfx, DT) in ((:ARGB8888ToRGB888, UInt8), (:RGBA8888ToRGB888, UInt8), (:BGRA8888ToRGB888, UInt8),
                  (:ARGBFFFFToRGBFFF, Float32), (:RGBAFFFFToRGBFFF, Float32), (:BGRAFFFFToRGBFFF, Float32))
    T = DT
    sym  = Symbol("vImageFlatten_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractArray{$T,3}, src::AbstractArray{$T,3}, backgroundColor;
                       isImagePremultiplied::Bool = true, flags::Integer = kvImageNoFlags)
            bc = convert(Vector{$T}, collect(backgroundColor))
            GC.@preserve src dest bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{$T}, Bool, vImage_Flags),
                    sb, db, pointer(bc), isImagePremultiplied, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    flatten_ARGB8888ToRGB888!(dest, src, backgroundColor; isImagePremultiplied=true, flags) -> dest

Composite a 4-channel image over an opaque `backgroundColor` (length-4 collection),
producing a 3-channel RGB image. Mutating variants for `ARGB/RGBA/BGRA` × `8888/FFFF`.
"""
flatten_ARGB8888ToRGB888!

# Deinterleave / interleave (4-planar <-> ARGB) --------------------------------------
for (sfx, T) in ((:ARGB8888toPlanar8, UInt8), (:ARGBFFFFtoPlanarF, Float32))
    sym = Symbol("vImageConvert_", sfx)
    @eval begin
        function $(_jlname(sym))(src::AbstractArray{$T,3}, destA::AbstractMatrix{$T}, destR::AbstractMatrix{$T},
                                 destG::AbstractMatrix{$T}, destB::AbstractMatrix{$T}; flags::Integer = kvImageNoFlags)
            GC.@preserve src destA destR destG destB begin
                sb = Ref(vimage_buffer(src)); a = Ref(vimage_buffer(destA)); r = Ref(vimage_buffer(destR))
                g = Ref(vimage_buffer(destG)); b = Ref(vimage_buffer(destB))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    sb, a, r, g, b, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return (destA, destR, destG, destB)
        end
    end
end
for (sfx, T) in ((:Planar8toARGB8888, UInt8), (:PlanarFtoARGBFFFF, Float32))
    sym = Symbol("vImageConvert_", sfx)
    @eval begin
        function $(_jlname(sym))(dest::AbstractArray{$T,3}, srcA::AbstractMatrix{$T}, srcR::AbstractMatrix{$T},
                                 srcG::AbstractMatrix{$T}, srcB::AbstractMatrix{$T}; flags::Integer = kvImageNoFlags)
            GC.@preserve dest srcA srcR srcG srcB begin
                a = Ref(vimage_buffer(srcA)); r = Ref(vimage_buffer(srcR)); g = Ref(vimage_buffer(srcG))
                b = Ref(vimage_buffer(srcB)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    a, r, g, b, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    convert_ARGB8888toPlanar8(src, destA, destR, destG, destB) -> (destA,destR,destG,destB)
    convert_Planar8toARGB8888(dest, srcA, srcR, srcG, srcB) -> dest

Deinterleave an interleaved image into four planar channel buffers, or interleave four
planes into one buffer. Also available for the float variants
`convert_ARGBFFFFtoPlanarF` / `convert_PlanarFtoARGBFFFF`.
"""
convert_ARGB8888toPlanar8

# RGB101010 packed formats: (src, dest, rangeMin, rangeMax, permuteMap, flags).
# Packed 2101010 pixels are represented as a planar Matrix{UInt32}.
const _CONVERT_1010102 = (
    :vImageConvert_RGBA1010102ToARGB8888, :vImageConvert_ARGB8888ToRGBA1010102,
    :vImageConvert_RGBA1010102ToARGB16Q12, :vImageConvert_RGBA1010102ToARGB16U,
    :vImageConvert_ARGB16UToRGBA1010102, :vImageConvert_ARGB2101010ToARGB8888,
    :vImageConvert_ARGB8888ToXRGB2101010, :vImageConvert_ARGB8888ToARGB2101010,
    :vImageConvert_ARGB2101010ToARGB16Q12, :vImageConvert_ARGB2101010ToARGB16U,
    :vImageConvert_ARGB16UToXRGB2101010, :vImageConvert_ARGB16UToARGB2101010,
    :vImageConvert_ARGB2101010ToARGBFFFF, :vImageConvert_ARGBFFFFToXRGB2101010,
    :vImageConvert_ARGBFFFFToARGB2101010, :vImageConvert_ARGB2101010ToARGB16F,
)
for sym in _CONVERT_1010102
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractArray, src::AbstractArray; rangeMin::Integer = 0,
                       rangeMax::Integer = 1023, permuteMap = (0, 1, 2, 3), flags::Integer = kvImageNoFlags)
            pm = convert(Vector{UInt8}, collect(permuteMap))
            GC.@preserve src dest pm begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Int32, Int32, Ptr{UInt8}, vImage_Flags),
                    sb, db, Int32(rangeMin), Int32(rangeMax), pointer(pm), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    convert_ARGB8888ToARGB2101010!(dest, src; rangeMin=0, rangeMax=1023, permuteMap=(0,1,2,3), flags) -> dest

Conversions to/from packed 2101010 (10-bit-per-channel, 32-bit) pixel formats. Packed
buffers are represented as a planar `Matrix{UInt32}`. `rangeMin`/`rangeMax` set the
integer range mapping to `[0,1]` for the 10-bit channels; `permuteMap` reorders
channels. Mutating variants only. See the docs page for the full list.
"""
convert_ARGB8888ToARGB2101010!

# =====================================================================================
# ADDITIONAL FAMILIES (extended coverage)
# =====================================================================================

# ---- AffineWarpCG : uses vImage_CGAffineTransform == Double transform on LP64 --------
for sfx in (:Planar8, :PlanarF, :ARGB8888, :ARGB16U, :ARGB16S, :ARGBFFFF)
    T, nch = _VFMT[sfx]
    sym   = Symbol("vImageAffineWarpCG_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    if nch == 1
        @eval function $bang(dest::AbstractArray, src::AbstractArray, transform::vImage_AffineTransform_Double;
                             backColor = zero($T), flags::Integer = kvImageBackgroundColorFill)
            tref = Ref(transform)
            GC.@preserve src dest tref begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{vImage_AffineTransform_Double}, $T, vImage_Flags),
                    sb, db, C_NULL, tref, $T(backColor), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    else
        @eval function $bang(dest::AbstractArray, src::AbstractArray, transform::vImage_AffineTransform_Double;
                             backColor = zeros($T, $nch), flags::Integer = kvImageBackgroundColorFill)
            tref = Ref(transform); bc = convert(Vector{$T}, collect(backColor))
            GC.@preserve src dest tref bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{vImage_AffineTransform_Double}, Ptr{$T}, vImage_Flags),
                    sb, db, C_NULL, tref, pointer(bc), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
    @eval $alloc(src::AbstractArray, transform::vImage_AffineTransform_Double; kw...) =
        $bang(_alloc_like(src), src, transform; kw...)
end
"""
    affineWarpCG_ARGB8888(src, transform::vImage_AffineTransform_Double; backColor, flags)

Affine warp using a CoreGraphics-compatible (`vImage_CGAffineTransform`, double
precision) transform. Formats `Planar8/F`, `ARGB8888/16U/16S/FFFF`.
"""
affineWarpCG_ARGB8888

# ---- Separable convolution ----------------------------------------------------------
for (sfx, BT, nch) in ((:Planar8, UInt16, 1), (:Planar16U, UInt16, 1), (:PlanarF, Cfloat, 1),
                       (:Planar16F, UInt16, 1), (:ARGB8888, UInt8, 4))
    sym   = Symbol("vImageSepConvolve_", sfx)
    bang  = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    if nch == 1
        @eval function $bang(dest::AbstractArray, src::AbstractArray,
                             kernelX::AbstractVector{<:Real}, kernelY::AbstractVector{<:Real};
                             bias::Real = 0, backColor = zero($BT),
                             srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageEdgeExtend)
            kx = convert(Vector{Cfloat}, kernelX); ky = convert(Vector{Cfloat}, kernelY)
            GC.@preserve src dest kx ky begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                     Ptr{Cfloat}, UInt32, Ptr{Cfloat}, UInt32, Cfloat, $BT, vImage_Flags),
                    sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                    pointer(kx), UInt32(length(kx)), pointer(ky), UInt32(length(ky)),
                    Cfloat(bias), $BT(backColor), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    else
        @eval function $bang(dest::AbstractArray, src::AbstractArray,
                             kernelX::AbstractVector{<:Real}, kernelY::AbstractVector{<:Real};
                             bias::Real = 0, backColor = zeros(UInt8, 4),
                             srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageEdgeExtend)
            kx = convert(Vector{Cfloat}, kernelX); ky = convert(Vector{Cfloat}, kernelY)
            bc = convert(Vector{UInt8}, collect(backColor))
            GC.@preserve src dest kx ky bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                     Ptr{Cfloat}, UInt32, Ptr{Cfloat}, UInt32, Cfloat, Ptr{UInt8}, vImage_Flags),
                    sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                    pointer(kx), UInt32(length(kx)), pointer(ky), UInt32(length(ky)),
                    Cfloat(bias), pointer(bc), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
    @eval $alloc(src::AbstractArray, kernelX::AbstractVector{<:Real}, kernelY::AbstractVector{<:Real}; kw...) =
        $bang(_alloc_like(src), src, kernelX, kernelY; kw...)
end
"""
    sepConvolve_PlanarF(src, kernelX, kernelY; bias=0, backColor, flags)

Separable convolution: convolve rows by `kernelX` and columns by `kernelY`. Cheaper
than a full 2-D kernel when the kernel is separable. Formats `Planar8/16U/16F/F`,
`ARGB8888`. Allocating and `!` variants.
"""
sepConvolve_PlanarF

# ---- Flatten (4-channel interleaved output) -----------------------------------------
for sfx in (:ARGB8888, :RGBA8888, :ARGB16U, :RGBA16U, :ARGB16Q12, :RGBA16Q12, :ARGBFFFF, :RGBAFFFF)
    T, nch = _VFMT[sfx]
    sym  = Symbol("vImageFlatten_", sfx)
    bang = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray{$T,3}, src::AbstractArray{$T,3}, backgroundColor;
                       isImagePremultiplied::Bool = true, flags::Integer = kvImageNoFlags)
            bc = convert(Vector{$T}, collect(backgroundColor))
            GC.@preserve src dest bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{$T}, Bool, vImage_Flags),
                    sb, db, pointer(bc), isImagePremultiplied, vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(src::AbstractArray{$T,3}, backgroundColor; kw...) =
            $bang(_alloc_like(src), src, backgroundColor; kw...)
    end
end
"""
    flatten_ARGB8888(src, backgroundColor; isImagePremultiplied=true, flags) -> array

Composite a 4-channel image over an opaque `backgroundColor` (length-4), keeping the
4-channel layout. Formats `ARGB/RGBA` × `8888/16U/16Q12/FFFF`. Allocating and `!` variants.
"""
flatten_ARGB8888

# ---- OverwriteChannels / SelectChannels : (newSrc, origSrc, dest, copyMask, flags) ---
for op in (:OverwriteChannels, :SelectChannels)
    for sfx in (:ARGB8888, :ARGBFFFF)
        sym  = Symbol("vImage", op, "_", sfx)
        bang = Symbol(_jlname(sym), "!")
        alloc = _jlname(sym)
        @eval begin
            function $bang(dest::AbstractArray, newSrc::AbstractArray, origSrc::AbstractArray,
                           copyMask::Integer; flags::Integer = kvImageNoFlags)
                GC.@preserve newSrc origSrc dest begin
                    nb = Ref(vimage_buffer(newSrc)); ob = Ref(vimage_buffer(origSrc)); db = Ref(vimage_buffer(dest))
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, UInt8, vImage_Flags),
                        nb, ob, db, UInt8(copyMask), vImage_Flags(flags))
                end
                _check(err, $(String(sym)))
                return dest
            end
            $alloc(newSrc::AbstractArray, origSrc::AbstractArray, copyMask::Integer; flags::Integer = kvImageNoFlags) =
                $bang(_alloc_like(origSrc), newSrc, origSrc, copyMask; flags = flags)
        end
    end
end
"""
    overwriteChannels_ARGB8888(newSrc, origSrc, copyMask; flags) -> array
    selectChannels_ARGBFFFF(newSrc, origSrc, copyMask; flags) -> array

Copy the channels selected by the bitmask `copyMask` from `newSrc`, and the rest from
`origSrc`, into the destination. Formats `ARGB8888`, `ARGBFFFF`.
"""
overwriteChannels_ARGB8888

# ---- OverwriteChannelsWithPixel : (pixel, src, dest, copyMask, flags) ----------------
for sfx in (:ARGB8888, :ARGB16U, :ARGBFFFF)
    T, nch = _VFMT[sfx]
    sym  = Symbol("vImageOverwriteChannelsWithPixel_", sfx)
    bang = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray, src::AbstractArray, pixel, copyMask::Integer;
                       flags::Integer = kvImageNoFlags)
            px = convert(Vector{$T}, collect(pixel))
            GC.@preserve src dest px begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{$T}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, UInt8, vImage_Flags),
                    pointer(px), sb, db, UInt8(copyMask), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(src::AbstractArray, pixel, copyMask::Integer; flags::Integer = kvImageNoFlags) =
            $bang(_alloc_like(src), src, pixel, copyMask; flags = flags)
    end
end
"""
    overwriteChannelsWithPixel_ARGB8888(src, pixel, copyMask; flags) -> array

Replace the channels selected by `copyMask` with the constant value from `pixel`
(length-4), copying the rest from `src`. Formats `ARGB8888/16U/FFFF`.
"""
overwriteChannelsWithPixel_ARGB8888

# ---- OverwriteChannelsWithScalar (interleaved): (scalar, src, dest, copyMask, flags) -
for (sfx, ST) in ((:ARGB8888, UInt8), (:ARGBFFFF, Cfloat))
    T, _ = _VFMT[sfx]
    sym  = Symbol("vImageOverwriteChannelsWithScalar_", sfx)
    bang = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray{$T,3}, src::AbstractArray{$T,3}, scalar, copyMask::Integer;
                       flags::Integer = kvImageNoFlags)
            GC.@preserve src dest begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    ($ST, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, UInt8, vImage_Flags),
                    $ST(scalar), sb, db, UInt8(copyMask), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(src::AbstractArray{$T,3}, scalar, copyMask::Integer; flags::Integer = kvImageNoFlags) =
            $bang(_alloc_like(src), src, scalar, copyMask; flags = flags)
    end
end

# ---- PermuteChannelsWithMaskedInsert : (src,dest,permuteMap,copyMask,bgColor,flags) --
for sfx in (:ARGB8888, :ARGB16U, :ARGBFFFF)
    T, _ = _VFMT[sfx]
    sym  = Symbol("vImagePermuteChannelsWithMaskedInsert_", sfx)
    bang = Symbol(_jlname(sym), "!")
    alloc = _jlname(sym)
    @eval begin
        function $bang(dest::AbstractArray, src::AbstractArray, permuteMap, copyMask::Integer;
                       backgroundColor = zeros($T, 4), flags::Integer = kvImageNoFlags)
            pm = convert(Vector{UInt8}, collect(permuteMap)); bc = convert(Vector{$T}, collect(backgroundColor))
            GC.@preserve src dest pm bc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, UInt8, Ptr{$T}, vImage_Flags),
                    sb, db, pointer(pm), UInt8(copyMask), pointer(bc), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
        $alloc(src::AbstractArray, permuteMap, copyMask::Integer; kw...) =
            $bang(_alloc_like(src), src, permuteMap, copyMask; kw...)
    end
end
"""
    permuteChannelsWithMaskedInsert_ARGB8888(src, permuteMap, copyMask; backgroundColor, flags) -> array

Permute channels per `permuteMap`, but for channels selected by `copyMask` insert the
constant `backgroundColor` instead. Formats `ARGB8888/16U/FFFF`.
"""
permuteChannelsWithMaskedInsert_ARGB8888

# ---- TableLookUp_ARGB8888 : (src,dest,aTable,rTable,gTable,bTable,flags) --------------
function tableLookUp_ARGB8888!(dest::AbstractArray{UInt8,3}, src::AbstractArray{UInt8,3},
                               alphaTable, redTable, greenTable, blueTable; flags::Integer = kvImageNoFlags)
    a = convert(Vector{UInt8}, collect(alphaTable)); r = convert(Vector{UInt8}, collect(redTable))
    g = convert(Vector{UInt8}, collect(greenTable)); b = convert(Vector{UInt8}, collect(blueTable))
    all(t -> length(t) == 256, (a, r, g, b)) || throw(DimensionMismatch("each table must have 256 entries"))
    GC.@preserve src dest a r g b begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageTableLookUp_ARGB8888, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, vImage_Flags),
            sb, db, pointer(a), pointer(r), pointer(g), pointer(b), vImage_Flags(flags))
    end
    _check(err, "vImageTableLookUp_ARGB8888")
    return dest
end
"""
    tableLookUp_ARGB8888(src, alphaTable, redTable, greenTable, blueTable; flags) -> array

Independently remap each channel of an ARGB8888 image through its own 256-entry table.
Allocating and `!` variants.
"""
tableLookUp_ARGB8888(src::AbstractArray{UInt8,3}, a, r, g, b; kw...) =
    tableLookUp_ARGB8888!(similar(src), src, a, r, g, b; kw...)

# ---- Extra lookup tables ------------------------------------------------------------
function lookupTable_8to64U!(dest::AbstractMatrix{UInt64}, src::AbstractMatrix{UInt8}, table; flags::Integer = kvImageNoFlags)
    t = convert(Vector{UInt64}, collect(table)); length(t) == 256 || throw(DimensionMismatch("table must have 256 entries"))
    GC.@preserve src dest t begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageLookupTable_8to64U, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt64}, vImage_Flags), sb, db, pointer(t), vImage_Flags(flags))
    end
    _check(err, "vImageLookupTable_8to64U"); return dest
end
function lookupTable_PlanarFtoPlanar8!(dest::AbstractMatrix{UInt8}, src::AbstractMatrix{Float32}, table; flags::Integer = kvImageNoFlags)
    t = convert(Vector{UInt8}, collect(table)); length(t) == 4096 || throw(DimensionMismatch("table must have 4096 entries"))
    GC.@preserve src dest t begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageLookupTable_PlanarFtoPlanar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, vImage_Flags), sb, db, pointer(t), vImage_Flags(flags))
    end
    _check(err, "vImageLookupTable_PlanarFtoPlanar8"); return dest
end
function lookupTable_Planar16!(dest::AbstractMatrix{UInt16}, src::AbstractMatrix{UInt16}, table; flags::Integer = kvImageNoFlags)
    t = convert(Vector{UInt16}, collect(table)); length(t) == 0x10000 || throw(DimensionMismatch("table must have 65536 entries"))
    GC.@preserve src dest t begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageLookupTable_Planar16, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt16}, vImage_Flags), sb, db, pointer(t), vImage_Flags(flags))
    end
    _check(err, "vImageLookupTable_Planar16"); return dest
end
"""
    lookupTable_Planar16!(dest, src, table) -> dest

Additional lookup-table mappings: `lookupTable_8to64U!` (256-entry UInt64 table),
`lookupTable_PlanarFtoPlanar8!` (4096-entry UInt8 table over `[0,1]`), and
`lookupTable_Planar16!` (65536-entry UInt16 table).
"""
lookupTable_Planar16!

# ---- SymmetricPiecewiseGamma + cross-type PiecewiseGamma -----------------------------
for (sym, T, DT, BT) in (
        (:vImageSymmetricPiecewiseGamma_PlanarF, Float32, Float32, Cfloat),
        (:vImageSymmetricPiecewiseGamma_Planar16Q12, Int16, Int16, Int16),
        (:vImagePiecewiseGamma_Planar8toPlanarF, UInt8, Float32, UInt8),
        (:vImagePiecewiseGamma_Planar8toPlanar16Q12, UInt8, Int16, UInt8),
        (:vImagePiecewiseGamma_Planar16Q12toPlanar8, Int16, UInt8, Int16),
        (:vImagePiecewiseGamma_PlanarFtoPlanar8, Float32, UInt8, Cfloat))
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractMatrix{$DT}, src::AbstractMatrix{$T};
                       exponentialCoeffs = (1f0, 0f0, 0f0), gamma::Real = 1f0,
                       linearCoeffs = (1f0, 0f0), boundary = zero($BT), flags::Integer = kvImageNoFlags)
            ec = convert(Vector{Cfloat}, collect(exponentialCoeffs)); lc = convert(Vector{Cfloat}, collect(linearCoeffs))
            GC.@preserve src dest ec lc begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cfloat}, Cfloat, Ptr{Cfloat}, $BT, vImage_Flags),
                    sb, db, pointer(ec), Cfloat(gamma), pointer(lc), $BT(boundary), vImage_Flags(flags))
            end
            _check(err, $(String(sym)))
            return dest
        end
    end
end
"""
    symmetricPiecewiseGamma_PlanarF!(dest, src; exponentialCoeffs, gamma, linearCoeffs, boundary, flags)

Symmetric piecewise gamma (odd-symmetric about 0), plus cross-type piecewise gamma
mutating variants `piecewiseGamma_Planar8toPlanarF!`, `…Planar8toPlanar16Q12!`,
`…Planar16Q12toPlanar8!`, `…PlanarFtoPlanar8!`.
"""
symmetricPiecewiseGamma_PlanarF!

# ---- MatrixMultiply -> planar output ------------------------------------------------
function matrixMultiply_ARGB8888ToPlanar8!(dest::AbstractMatrix{UInt8}, src::AbstractArray{UInt8,3},
                                           matrix; divisor::Integer = 256, preBias = nothing,
                                           postBias::Integer = 0, flags::Integer = kvImageNoFlags)
    m = convert(Vector{Int16}, collect(matrix)); length(m) == 4 || throw(DimensionMismatch("matrix must have 4 entries"))
    prev = preBias === nothing ? Int16[] : convert(Vector{Int16}, collect(preBias))
    GC.@preserve src dest m prev begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        prep = isempty(prev) ? Ptr{Int16}(C_NULL) : pointer(prev)
        err = ccall((:vImageMatrixMultiply_ARGB8888ToPlanar8, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Int16}, Int32, Ptr{Int16}, Int32, vImage_Flags),
            sb, db, pointer(m), Int32(divisor), prep, Int32(postBias), vImage_Flags(flags))
    end
    _check(err, "vImageMatrixMultiply_ARGB8888ToPlanar8"); return dest
end
function matrixMultiply_ARGBFFFFToPlanarF!(dest::AbstractMatrix{Float32}, src::AbstractArray{Float32,3},
                                           matrix; preBias = nothing, postBias::Real = 0f0, flags::Integer = kvImageNoFlags)
    m = convert(Vector{Cfloat}, collect(matrix)); length(m) == 4 || throw(DimensionMismatch("matrix must have 4 entries"))
    prev = preBias === nothing ? Cfloat[] : convert(Vector{Cfloat}, collect(preBias))
    GC.@preserve src dest m prev begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        prep = isempty(prev) ? Ptr{Cfloat}(C_NULL) : pointer(prev)
        err = ccall((:vImageMatrixMultiply_ARGBFFFFToPlanarF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cfloat}, Ptr{Cfloat}, Cfloat, vImage_Flags),
            sb, db, pointer(m), prep, Cfloat(postBias), vImage_Flags(flags))
    end
    _check(err, "vImageMatrixMultiply_ARGBFFFFToPlanarF"); return dest
end
"""
    matrixMultiply_ARGB8888ToPlanar8!(dest, src, matrix; divisor=256, preBias, postBias, flags) -> dest

Reduce a 4-channel image to one planar channel via a length-4 dot-product `matrix`
(e.g. RGB→luminance). `!`-only; `matrixMultiply_ARGBFFFFToPlanarF!` is the float form.
"""
matrixMultiply_ARGB8888ToPlanar8!

# ---- Deinterleave / interleave: 3-plane RGB and 4-plane (16U) ------------------------
for (sfx, T) in ((:Planar8toRGB888, UInt8), (:PlanarFtoRGBFFF, Float32), (:Planar16UtoRGB16U, UInt16))
    sym = Symbol("vImageConvert_", sfx)
    @eval function $(_jlname(sym))(dest::AbstractArray{$T,3}, red::AbstractMatrix{$T}, green::AbstractMatrix{$T},
                                   blue::AbstractMatrix{$T}; flags::Integer = kvImageNoFlags)
        GC.@preserve dest red green blue begin
            r = Ref(vimage_buffer(red)); g = Ref(vimage_buffer(green)); b = Ref(vimage_buffer(blue)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                r, g, b, db, vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return dest
    end
end
for (sfx, T) in ((:RGB888toPlanar8, UInt8), (:RGBFFFtoPlanarF, Float32), (:RGB16UtoPlanar16U, UInt16))
    sym = Symbol("vImageConvert_", sfx)
    @eval function $(_jlname(sym))(src::AbstractArray{$T,3}, red::AbstractMatrix{$T}, green::AbstractMatrix{$T},
                                   blue::AbstractMatrix{$T}; flags::Integer = kvImageNoFlags)
        GC.@preserve src red green blue begin
            sb = Ref(vimage_buffer(src)); r = Ref(vimage_buffer(red)); g = Ref(vimage_buffer(green)); b = Ref(vimage_buffer(blue))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                sb, r, g, b, vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return (red, green, blue)
    end
end
function convert_Planar16UtoARGB16U(dest::AbstractArray{UInt16,3}, a::AbstractMatrix{UInt16}, r::AbstractMatrix{UInt16},
                                    g::AbstractMatrix{UInt16}, b::AbstractMatrix{UInt16}; flags::Integer = kvImageNoFlags)
    GC.@preserve dest a r g b begin
        ab = Ref(vimage_buffer(a)); rb = Ref(vimage_buffer(r)); gb = Ref(vimage_buffer(g)); bb = Ref(vimage_buffer(b)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageConvert_Planar16UtoARGB16U, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            ab, rb, gb, bb, db, vImage_Flags(flags))
    end
    _check(err, "vImageConvert_Planar16UtoARGB16U"); return dest
end
function convert_ARGB16UtoPlanar16U(src::AbstractArray{UInt16,3}, a::AbstractMatrix{UInt16}, r::AbstractMatrix{UInt16},
                                    g::AbstractMatrix{UInt16}, b::AbstractMatrix{UInt16}; flags::Integer = kvImageNoFlags)
    GC.@preserve src a r g b begin
        sb = Ref(vimage_buffer(src)); ab = Ref(vimage_buffer(a)); rb = Ref(vimage_buffer(r)); gb = Ref(vimage_buffer(g)); bb = Ref(vimage_buffer(b))
        err = ccall((:vImageConvert_ARGB16UtoPlanar16U, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            sb, ab, rb, gb, bb, vImage_Flags(flags))
    end
    _check(err, "vImageConvert_ARGB16UtoPlanar16U"); return (a, r, g, b)
end
"""
    convert_Planar8toRGB888(dest, red, green, blue) -> dest
    convert_RGB888toPlanar8(src, red, green, blue) -> (red, green, blue)

Interleave three planar channels into an RGB image, or split an RGB image into planes.
Variants for `Planar8/RGB888`, `PlanarF/RGBFFF`, `Planar16U/RGB16U`, plus 4-plane
`convert_Planar16UtoARGB16U` / `convert_ARGB16UtoPlanar16U`.
"""
convert_Planar8toRGB888

# ---- RGB (+alpha) -> 4-channel converts ---------------------------------------------
for sfx in (:RGB565toARGB8888, :RGB565toRGBA8888, :RGB565toBGRA8888)
    sym = Symbol("vImageConvert_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval function $bang(dest::AbstractArray{UInt8,3}, src::AbstractMatrix{UInt16}; alpha::Integer = 255,
                         flags::Integer = kvImageNoFlags)
        GC.@preserve src dest begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (UInt8, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                UInt8(alpha), sb, db, vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return dest
    end
end
for (sfx, T, AT) in ((:RGB888toARGB8888, UInt8, UInt8), (:RGB888toRGBA8888, UInt8, UInt8), (:RGB888toBGRA8888, UInt8, UInt8),
                     (:RGBFFFtoARGBFFFF, Float32, Cfloat), (:RGBFFFtoRGBAFFFF, Float32, Cfloat), (:RGBFFFtoBGRAFFFF, Float32, Cfloat),
                     (:RGB16UtoARGB16U, UInt16, UInt16), (:RGB16UtoRGBA16U, UInt16, UInt16), (:RGB16UtoBGRA16U, UInt16, UInt16))
    sym = Symbol("vImageConvert_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval function $bang(dest::AbstractArray{$T,3}, rgbSrc::AbstractArray{$T,3};
                         alpha = one($AT), alphaPlane::Union{Nothing,AbstractMatrix{$T}} = nothing,
                         premultiply::Bool = false, flags::Integer = kvImageNoFlags)
        if alphaPlane === nothing
            GC.@preserve rgbSrc dest begin
                sb = Ref(vimage_buffer(rgbSrc)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, $AT, Ptr{vImage_Buffer}, Bool, vImage_Flags),
                    sb, C_NULL, $AT(alpha), db, premultiply, vImage_Flags(flags))
            end
        else
            GC.@preserve rgbSrc dest alphaPlane begin
                sb = Ref(vimage_buffer(rgbSrc)); ab = Ref(vimage_buffer(alphaPlane)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, $AT, Ptr{vImage_Buffer}, Bool, vImage_Flags),
                    sb, ab, $AT(alpha), db, premultiply, vImage_Flags(flags))
            end
        end
        _check(err, $(String(sym))); return dest
    end
end
"""
    convert_RGB888toARGB8888!(dest, rgbSrc; alpha=1, alphaPlane=nothing, premultiply=false, flags)
    convert_RGB565toARGB8888!(dest, src; alpha=255, flags)

Add an alpha channel to a 3-channel image, producing a 4-channel result. A constant
`alpha` (or an `alphaPlane`) is inserted; `premultiply` optionally premultiplies. `!`-only.
"""
convert_RGB888toARGB8888!

# ---- Permute + masked convert (bit-depth) : (src,dest,permuteMap,copyMask,bg,flags) --
for (sfx, T) in ((:ARGB16UToARGB8888, UInt8), (:ARGB8888ToARGB16U, UInt16), (:RGB16UToARGB8888, UInt8))
    sym = Symbol("vImageConvert_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval function $bang(dest::AbstractArray, src::AbstractArray; permuteMap = (0, 1, 2, 3), copyMask::Integer = 0,
                         backgroundColor = zeros($T, 4), flags::Integer = kvImageNoFlags)
        pm = convert(Vector{UInt8}, collect(permuteMap)); bc = convert(Vector{$T}, collect(backgroundColor))
        GC.@preserve src dest pm bc begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, UInt8, Ptr{$T}, vImage_Flags),
                sb, db, pointer(pm), UInt8(copyMask), pointer(bc), vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return dest
    end
end
"""
    convert_ARGB16UToARGB8888!(dest, src; permuteMap=(0,1,2,3), copyMask=0, backgroundColor, flags)

Bit-depth conversion between 8- and 16-bit interleaved 4-channel formats with optional
channel permutation and masked background insertion. `!`-only variants
`convert_ARGB16UToARGB8888!`, `convert_ARGB8888ToARGB16U!`, `convert_RGB16UToARGB8888!`.
"""
convert_ARGB16UToARGB8888!

# ---- 4-channel float<->int with per-channel maxFloat/minFloat -----------------------
function convert_ARGB8888toPlanarF(src::AbstractArray{UInt8,3}, a::AbstractMatrix{Float32}, r::AbstractMatrix{Float32},
                                   g::AbstractMatrix{Float32}, b::AbstractMatrix{Float32};
                                   maxFloat = (1f0,1f0,1f0,1f0), minFloat = (0f0,0f0,0f0,0f0), flags::Integer = kvImageNoFlags)
    mx = convert(Vector{Cfloat}, collect(maxFloat)); mn = convert(Vector{Cfloat}, collect(minFloat))
    GC.@preserve src a r g b mx mn begin
        sb = Ref(vimage_buffer(src)); ab = Ref(vimage_buffer(a)); rb = Ref(vimage_buffer(r)); gb = Ref(vimage_buffer(g)); bb = Ref(vimage_buffer(b))
        err = ccall((:vImageConvert_ARGB8888toPlanarF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cfloat}, Ptr{Cfloat}, vImage_Flags),
            sb, ab, rb, gb, bb, pointer(mx), pointer(mn), vImage_Flags(flags))
    end
    _check(err, "vImageConvert_ARGB8888toPlanarF"); return (a, r, g, b)
end
"""
    convert_ARGB8888toPlanarF(src, a, r, g, b; maxFloat, minFloat, flags) -> (a,r,g,b)

Deinterleave an ARGB8888 image into four float planes, mapping `[0,255]` to
`[minFloat, maxFloat]` per channel (length-4 collections).
"""
convert_ARGB8888toPlanarF

# ---- RGB101010 packed with alpha ----------------------------------------------------
for (sfx, AT) in ((:XRGB2101010ToARGB8888, UInt8), (:XRGB2101010ToARGB16Q12, Int16),
                  (:XRGB2101010ToARGB16U, UInt16), (:XRGB2101010ToARGBFFFF, Cfloat),
                  (:XRGB2101010ToARGB16F, Cfloat))
    sym = Symbol("vImageConvert_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval function $bang(dest::AbstractArray, src::AbstractMatrix{UInt32}; alpha = zero($AT),
                         rangeMin::Integer = 0, rangeMax::Integer = 1023, permuteMap = (0,1,2,3),
                         flags::Integer = kvImageNoFlags)
        pm = convert(Vector{UInt8}, collect(permuteMap))
        GC.@preserve src dest pm begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, $AT, Ptr{vImage_Buffer}, Int32, Int32, Ptr{UInt8}, vImage_Flags),
                sb, $AT(alpha), db, Int32(rangeMin), Int32(rangeMax), pointer(pm), vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return dest
    end
end
"""
    convert_XRGB2101010ToARGB8888!(dest, src; alpha=0, rangeMin=0, rangeMax=1023, permuteMap, flags)

Expand a packed XRGB2101010 image (planar `Matrix{UInt32}`) into a 4-channel image,
supplying the alpha channel from the scalar `alpha`. `!`-only.
"""
convert_XRGB2101010ToARGB8888!

# =====================================================================================
# Y'CbCr ⇄ ARGB CONVERSION  (4:4:4 chroma layouts, 8- and 16-bit)
# =====================================================================================
#
# Y'CbCr conversion is a "setup + execute" family: a *matrix* + *pixel range* + input
# and output *format codes* are baked once into an opaque info struct by a
# `GenerateConversion` call, and that info struct is then handed to the per-pixel
# convert routines. See the struct byte-layouts (verified against vImage_Types.h and
# empirically):
#   * vImage_YpCbCrToARGBMatrix : 5 × Float32           (20 bytes)
#   * vImage_ARGBToYpCbCrMatrix : 8 × Float32           (32 bytes)
#   * vImage_YpCbCrPixelRange   : 8 × Int32             (32 bytes)
#   * vImage_YpCbCrToARGB / vImage_ARGBToYpCbCr : opaque[128], 16-byte aligned.
# The opaque info blobs are represented as NTuple{8,UInt128} (128 bytes, 16-aligned).

"""
    vImage_YpCbCrToARGBMatrix(Yp, Cr_R, Cr_G, Cb_G, Cb_B)

Sparse 3×3 Y'CbCr→RGB conversion matrix (5 non-trivial Float32 coefficients). Feed to
[`convert_YpCbCrToARGB_GenerateConversion`](@ref AppleAccelerate.convert_YpCbCrToARGB_GenerateConversion).
The standard ITU matrices are available as
[`kvImage_YpCbCrToARGBMatrix_ITU_R_601_4`](@ref AppleAccelerate.kvImage_YpCbCrToARGBMatrix_ITU_R_601_4)
and `…_709_2`.
"""
struct vImage_YpCbCrToARGBMatrix
    Yp::Cfloat; Cr_R::Cfloat; Cr_G::Cfloat; Cb_G::Cfloat; Cb_B::Cfloat
end

"""
    vImage_ARGBToYpCbCrMatrix(R_Yp, G_Yp, B_Yp, R_Cb, G_Cb, B_Cb_R_Cr, G_Cr, B_Cr)

3×3 RGB→Y'CbCr conversion matrix (8 Float32 coefficients). Inverse of
[`vImage_YpCbCrToARGBMatrix`](@ref). Standard matrices:
[`kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4`](@ref AppleAccelerate.kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4)
and `…_709_2`.
"""
struct vImage_ARGBToYpCbCrMatrix
    R_Yp::Cfloat; G_Yp::Cfloat; B_Yp::Cfloat; R_Cb::Cfloat
    G_Cb::Cfloat; B_Cb_R_Cr::Cfloat; G_Cr::Cfloat; B_Cr::Cfloat
end

"""
    vImage_YpCbCrPixelRange(Yp_bias, CbCr_bias, YpRangeMax, CbCrRangeMax, YpMax, YpMin, CbCrMax, CbCrMin)

Range and clamping information for a Y'CbCr pixel format (8 × `Int32`). Handy presets:
`kvImageYpCbCrPixelRange_VideoRange_8bit_Clamped`,
`kvImageYpCbCrPixelRange_FullRange_8bit_Clamped`.
"""
struct vImage_YpCbCrPixelRange
    Yp_bias::Int32; CbCr_bias::Int32; YpRangeMax::Int32; CbCrRangeMax::Int32
    YpMax::Int32; YpMin::Int32; CbCrMax::Int32; CbCrMin::Int32
end

"""
    vImage_YpCbCrToARGB

Opaque 128-byte (16-byte-aligned) info blob produced by
[`convert_YpCbCrToARGB_GenerateConversion`](@ref AppleAccelerate.convert_YpCbCrToARGB_GenerateConversion)
and consumed by the `convert_*ToARGB*` Y'CbCr routines.
"""
struct vImage_YpCbCrToARGB
    opaque::NTuple{8,UInt128}
end

"""
    vImage_ARGBToYpCbCr

Opaque 128-byte info blob produced by
[`convert_ARGBToYpCbCr_GenerateConversion`](@ref AppleAccelerate.convert_ARGBToYpCbCr_GenerateConversion)
and consumed by the `convert_ARGB*To444*` Y'CbCr routines.
"""
struct vImage_ARGBToYpCbCr
    opaque::NTuple{8,UInt128}
end

# --- vImageARGBType / vImageYpCbCrType format codes --------------------------
const kvImageARGB8888   = Cint(0)
const kvImageARGB16U    = Cint(1)
const kvImageARGB16Q12  = Cint(2)

const kvImage422CbYpCrYp8                  = Cint(0)
const kvImage422YpCbYpCr8                  = Cint(1)
const kvImage422CbYpCrYp8_AA8              = Cint(2)
const kvImage420Yp8_Cb8_Cr8                = Cint(3)
const kvImage420Yp8_CbCr8                  = Cint(4)
const kvImage444AYpCbCr8                   = Cint(5)
const kvImage444CrYpCb8                    = Cint(6)
const kvImage444CbYpCrA8                   = Cint(7)
const kvImage444CrYpCb10                   = Cint(8)
const kvImage422CrYpCbYpCbYpCbYpCrYpCrYp10 = Cint(9)
const kvImage422CbYpCrYp16                 = Cint(13)
const kvImage444AYpCbCr16                  = Cint(14)

# --- Standard conversion matrix constants (loaded from the framework) --------
@inline function _load_const(::Type{T}, name::Symbol) where {T}
    pp = cglobal((name, vimage_lib), Ptr{T})   # symbol is a `const T *` variable
    return unsafe_load(unsafe_load(pp))
end

"""
    kvImage_YpCbCrToARGBMatrix_ITU_R_601_4 :: vImage_YpCbCrToARGBMatrix

ITU-R BT.601-4 Y'CbCr→RGB conversion matrix (standard-definition video). `…_709_2` is
the BT.709-2 (high-definition) variant.
"""
const kvImage_YpCbCrToARGBMatrix_ITU_R_601_4 = _load_const(vImage_YpCbCrToARGBMatrix, :kvImage_YpCbCrToARGBMatrix_ITU_R_601_4)
const kvImage_YpCbCrToARGBMatrix_ITU_R_709_2 = _load_const(vImage_YpCbCrToARGBMatrix, :kvImage_YpCbCrToARGBMatrix_ITU_R_709_2)

"""
    kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4 :: vImage_ARGBToYpCbCrMatrix

ITU-R BT.601-4 RGB→Y'CbCr conversion matrix. `…_709_2` is the BT.709-2 variant.
"""
const kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4 = _load_const(vImage_ARGBToYpCbCrMatrix, :kvImage_ARGBToYpCbCrMatrix_ITU_R_601_4)
const kvImage_ARGBToYpCbCrMatrix_ITU_R_709_2 = _load_const(vImage_ARGBToYpCbCrMatrix, :kvImage_ARGBToYpCbCrMatrix_ITU_R_709_2)

# Convenience pixel-range presets (8-bit).
const kvImageYpCbCrPixelRange_VideoRange_8bit_Clamped = vImage_YpCbCrPixelRange(16, 128, 235, 240, 235, 16, 240, 16)
const kvImageYpCbCrPixelRange_FullRange_8bit_Clamped  = vImage_YpCbCrPixelRange(0, 128, 255, 255, 255, 1, 255, 0)

"""
    convert_YpCbCrToARGB_GenerateConversion(matrix, pixelRange, ypCbCrType, argbType; flags) -> vImage_YpCbCrToARGB

Bake a `matrix::vImage_YpCbCrToARGBMatrix` + `pixelRange::vImage_YpCbCrPixelRange` and
the input `ypCbCrType` / output `argbType` format codes into a reusable
[`vImage_YpCbCrToARGB`](@ref) info blob for the `convert_*ToARGB*` routines. Reuse the
returned info across many conversions rather than regenerating it.
"""
function convert_YpCbCrToARGB_GenerateConversion(matrix::vImage_YpCbCrToARGBMatrix,
        pixelRange::vImage_YpCbCrPixelRange, ypCbCrType::Integer, argbType::Integer;
        flags::Integer = kvImageNoFlags)
    out = Ref{vImage_YpCbCrToARGB}(); rm = Ref(matrix); rp = Ref(pixelRange)
    GC.@preserve out rm rp begin
        err = ccall((:vImageConvert_YpCbCrToARGB_GenerateConversion, vimage_lib), vImage_Error,
            (Ptr{vImage_YpCbCrToARGBMatrix}, Ptr{vImage_YpCbCrPixelRange}, Ptr{vImage_YpCbCrToARGB}, Cint, Cint, vImage_Flags),
            rm, rp, out, Cint(ypCbCrType), Cint(argbType), vImage_Flags(flags))
    end
    _check(err, "vImageConvert_YpCbCrToARGB_GenerateConversion")
    return out[]
end

"""
    convert_ARGBToYpCbCr_GenerateConversion(matrix, pixelRange, argbType, ypCbCrType; flags) -> vImage_ARGBToYpCbCr

Bake a `matrix::vImage_ARGBToYpCbCrMatrix` + `pixelRange` and the input `argbType` /
output `ypCbCrType` format codes into a reusable [`vImage_ARGBToYpCbCr`](@ref) info blob
for the `convert_ARGB*To444*` routines.
"""
function convert_ARGBToYpCbCr_GenerateConversion(matrix::vImage_ARGBToYpCbCrMatrix,
        pixelRange::vImage_YpCbCrPixelRange, argbType::Integer, ypCbCrType::Integer;
        flags::Integer = kvImageNoFlags)
    out = Ref{vImage_ARGBToYpCbCr}(); rm = Ref(matrix); rp = Ref(pixelRange)
    GC.@preserve out rm rp begin
        err = ccall((:vImageConvert_ARGBToYpCbCr_GenerateConversion, vimage_lib), vImage_Error,
            (Ptr{vImage_ARGBToYpCbCrMatrix}, Ptr{vImage_YpCbCrPixelRange}, Ptr{vImage_ARGBToYpCbCr}, Cint, Cint, vImage_Flags),
            rm, rp, out, Cint(argbType), Cint(ypCbCrType), vImage_Flags(flags))
    end
    _check(err, "vImageConvert_ARGBToYpCbCr_GenerateConversion")
    return out[]
end

# ---- Y'CbCr -> ARGB conversions (info::vImage_YpCbCrToARGB) --------------------------
# (sym, hasAlphaArg, alphaJuliaType)
for (sym, hasA, AT) in ((:vImageConvert_444AYpCbCr8ToARGB8888,  false, UInt8),
                        (:vImageConvert_444CbYpCrA8ToARGB8888,  false, UInt8),
                        (:vImageConvert_444CrYpCb8ToARGB8888,   true,  UInt8),
                        (:vImageConvert_444AYpCbCr16ToARGB8888, false, UInt8),
                        (:vImageConvert_444AYpCbCr16ToARGB16U,  false, UInt16))
    bang = Symbol(_jlname(sym), "!")
    if hasA
        @eval function $bang(dest::AbstractArray, src::AbstractArray, info::vImage_YpCbCrToARGB;
                             permuteMap = (0, 1, 2, 3), alpha = zero($AT), flags::Integer = kvImageNoFlags)
            pm = convert(Vector{UInt8}, collect(permuteMap)); ri = Ref(info)
            GC.@preserve src dest pm ri begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_YpCbCrToARGB}, Ptr{UInt8}, $AT, vImage_Flags),
                    sb, db, ri, pointer(pm), $AT(alpha), vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
    else
        @eval function $bang(dest::AbstractArray, src::AbstractArray, info::vImage_YpCbCrToARGB;
                             permuteMap = (0, 1, 2, 3), flags::Integer = kvImageNoFlags)
            pm = convert(Vector{UInt8}, collect(permuteMap)); ri = Ref(info)
            GC.@preserve src dest pm ri begin
                sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_YpCbCrToARGB}, Ptr{UInt8}, vImage_Flags),
                    sb, db, ri, pointer(pm), vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
    end
end

# ---- ARGB -> Y'CbCr conversions (info::vImage_ARGBToYpCbCr) --------------------------
for sym in (:vImageConvert_ARGB8888To444AYpCbCr8, :vImageConvert_ARGB8888To444CbYpCrA8,
            :vImageConvert_ARGB8888To444CrYpCb8, :vImageConvert_ARGB8888To444AYpCbCr16,
            :vImageConvert_ARGB16UTo444AYpCbCr16)
    bang = Symbol(_jlname(sym), "!")
    @eval function $bang(dest::AbstractArray, src::AbstractArray, info::vImage_ARGBToYpCbCr;
                         permuteMap = (0, 1, 2, 3), flags::Integer = kvImageNoFlags)
        pm = convert(Vector{UInt8}, collect(permuteMap)); ri = Ref(info)
        GC.@preserve src dest pm ri begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_ARGBToYpCbCr}, Ptr{UInt8}, vImage_Flags),
                sb, db, ri, pointer(pm), vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return dest
    end
end

"""
    convert_444CrYpCb8ToARGB8888!(dest, src, info; permuteMap=(0,1,2,3), alpha=0, flags)
    convert_ARGB8888To444CrYpCb8!(dest, src, info; permuteMap=(0,1,2,3), flags)

4:4:4 Y'CbCr ⇄ ARGB conversions. Build `info` once with
[`convert_YpCbCrToARGB_GenerateConversion`](@ref AppleAccelerate.convert_YpCbCrToARGB_GenerateConversion)
(for `*ToARGB*`) or
[`convert_ARGBToYpCbCr_GenerateConversion`](@ref AppleAccelerate.convert_ARGBToYpCbCr_GenerateConversion)
(for `ARGB*To444*`), then reuse it. Images are interleaved `Array{T,3}` sized
`(channels, width, height)`: 3-channel `v308` (`444CrYpCb8`), 4-channel `v408`/`y408`
(`444CbYpCrA8` / `444AYpCbCr8`) and 4-channel 16-bit `y416` (`444AYpCbCr16`). `permuteMap`
reorders the ARGB channels; `alpha` supplies the constant alpha when the Y'CbCr side has
none. `!`-only. Wrapped: `convert_444AYpCbCr8ToARGB8888!`, `convert_444CbYpCrA8ToARGB8888!`,
`convert_444CrYpCb8ToARGB8888!`, `convert_444AYpCbCr16ToARGB8888!`,
`convert_444AYpCbCr16ToARGB16U!`, and the five `convert_ARGB*To444*!` inverses.
"""
convert_444CrYpCb8ToARGB8888!

# =====================================================================================
# MULTI-KERNEL / FLOAT-KERNEL CONVOLUTION
# =====================================================================================

"""
    convolveMultiKernel_ARGB8888!(dest, src, kernels; divisors, biases, backgroundColor,
                                  srcOffsetX=0, srcOffsetY=0, flags=kvImageEdgeExtend) -> dest

Convolve a 4-channel ARGB8888 image applying a **separate** integer `kernel` to each
channel. `kernels` is a length-4 collection of equal-size `(kernelHeight, kernelWidth)`
integer matrices (both dims odd); `divisors` and `biases` are length-4. One edging-mode
flag is required (`kvImageEdgeExtend`, `kvImageBackgroundColorFill`, `kvImageCopyInPlace`
or `kvImageTruncateKernel`). Allocating and `!` variants. With four identical kernels this
matches [`convolve_Planar8`](@ref).
"""
function convolveMultiKernel_ARGB8888!(dest::AbstractArray{UInt8,3}, src::AbstractArray{UInt8,3},
        kernels; divisors = (1, 1, 1, 1), biases = (0, 0, 0, 0),
        backgroundColor = (0, 0, 0, 0), srcOffsetX::Integer = 0, srcOffsetY::Integer = 0,
        flags::Integer = kvImageEdgeExtend)
    ks = [convert(Matrix{Int16}, k) for k in kernels]
    length(ks) == 4 || throw(DimensionMismatch("need 4 kernels"))
    kh, kw = size(ks[1])
    all(size(k) == (kh, kw) for k in ks) || throw(DimensionMismatch("all kernels must match in size"))
    div = convert(Vector{Int32}, collect(divisors)); bia = convert(Vector{Int32}, collect(biases))
    bc = convert(Vector{UInt8}, collect(backgroundColor))
    GC.@preserve src dest ks div bia bc begin
        kptrs = Ptr{Int16}[pointer(k) for k in ks]
        GC.@preserve kptrs begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall((:vImageConvolveMultiKernel_ARGB8888, vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                 Ptr{Ptr{Int16}}, UInt32, UInt32, Ptr{Int32}, Ptr{Int32}, Ptr{UInt8}, vImage_Flags),
                sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                pointer(kptrs), UInt32(kh), UInt32(kw), pointer(div), pointer(bia), pointer(bc), vImage_Flags(flags))
        end
    end
    _check(err, "vImageConvolveMultiKernel_ARGB8888"); return dest
end
convolveMultiKernel_ARGB8888(src::AbstractArray{UInt8,3}, kernels; kw...) =
    convolveMultiKernel_ARGB8888!(similar(src), src, kernels; kw...)

"""
    convolveMultiKernel_ARGBFFFF!(dest, src, kernels; biases, backgroundColor,
                                  srcOffsetX=0, srcOffsetY=0, flags=kvImageEdgeExtend) -> dest

Float per-channel convolution of an ARGBFFFF image: `kernels` is a length-4 collection of
equal-size float matrices, `biases` length-4 (no divisor). See
[`convolveMultiKernel_ARGB8888!`](@ref). Allocating and `!` variants.
"""
function convolveMultiKernel_ARGBFFFF!(dest::AbstractArray{Float32,3}, src::AbstractArray{Float32,3},
        kernels; biases = (0f0, 0f0, 0f0, 0f0), backgroundColor = (0f0, 0f0, 0f0, 0f0),
        srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageEdgeExtend)
    ks = [convert(Matrix{Cfloat}, k) for k in kernels]
    length(ks) == 4 || throw(DimensionMismatch("need 4 kernels"))
    kh, kw = size(ks[1])
    all(size(k) == (kh, kw) for k in ks) || throw(DimensionMismatch("all kernels must match in size"))
    bia = convert(Vector{Cfloat}, collect(biases)); bc = convert(Vector{Cfloat}, collect(backgroundColor))
    GC.@preserve src dest ks bia bc begin
        kptrs = Ptr{Cfloat}[pointer(k) for k in ks]
        GC.@preserve kptrs begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall((:vImageConvolveMultiKernel_ARGBFFFF, vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
                 Ptr{Ptr{Cfloat}}, UInt32, UInt32, Ptr{Cfloat}, Ptr{Cfloat}, vImage_Flags),
                sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
                pointer(kptrs), UInt32(kh), UInt32(kw), pointer(bia), pointer(bc), vImage_Flags(flags))
        end
    end
    _check(err, "vImageConvolveMultiKernel_ARGBFFFF"); return dest
end
convolveMultiKernel_ARGBFFFF(src::AbstractArray{Float32,3}, kernels; kw...) =
    convolveMultiKernel_ARGBFFFF!(similar(src), src, kernels; kw...)

"""
    convolveFloatKernel_ARGB8888!(dest, src, kernel; bias=0, backColor=(0,0,0,0),
                                  srcOffsetX=0, srcOffsetY=0, flags) -> dest

Convolve an ARGB8888 image with a single **floating-point** `kernel` (higher precision
than the integer-kernel [`convolve_Planar8`](@ref)). `kernel` is
`(kernelHeight, kernelWidth)`. Allocating and `!` variants.
"""
function convolveFloatKernel_ARGB8888!(dest::AbstractArray{UInt8,3}, src::AbstractArray{UInt8,3},
        kernel::AbstractMatrix{<:Real}; bias::Real = 0, backColor = (0, 0, 0, 0),
        srcOffsetX::Integer = 0, srcOffsetY::Integer = 0, flags::Integer = kvImageNoFlags)
    k = convert(Matrix{Cfloat}, kernel); kh, kw = size(k)
    bc = convert(Vector{UInt8}, collect(backColor))
    GC.@preserve src dest k bc begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageConvolveFloatKernel_ARGB8888, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Csize_t, Csize_t,
             Ptr{Cfloat}, UInt32, UInt32, Cfloat, Ptr{UInt8}, vImage_Flags),
            sb, db, C_NULL, Csize_t(srcOffsetX), Csize_t(srcOffsetY),
            pointer(k), UInt32(kh), UInt32(kw), Cfloat(bias), pointer(bc), vImage_Flags(flags))
    end
    _check(err, "vImageConvolveFloatKernel_ARGB8888"); return dest
end
convolveFloatKernel_ARGB8888(src::AbstractArray{UInt8,3}, kernel::AbstractMatrix{<:Real}; kw...) =
    convolveFloatKernel_ARGB8888!(similar(src), src, kernel; kw...)

# =====================================================================================
# MULTI-PLANE MATRIX MULTIPLY
# =====================================================================================
# srcs[]/dests[] are C arrays of pointers-to-vImage_Buffer; matrix is (src_planes ×
# dest_planes) with out[j] = Σ_i in[i]·matrix[i,j], flattened row-major.

for (sym, T, KT, integer) in ((:vImageMatrixMultiply_Planar8,  UInt8, Int16,   true),
                              (:vImageMatrixMultiply_Planar16S, Int16, Int16,   true),
                              (:vImageMatrixMultiply_PlanarF,   Float32, Cfloat, false))
    fn = _jlname(sym)
    if integer
        @eval function $fn(srcs::AbstractVector{<:AbstractMatrix{$T}}, dests::AbstractVector{<:AbstractMatrix{$T}},
                           matrix::AbstractMatrix{<:Integer}; divisor::Integer = 1,
                           preBias = nothing, postBias = nothing, flags::Integer = kvImageNoFlags)
            sp = length(srcs); dp = length(dests)
            size(matrix) == (sp, dp) || throw(DimensionMismatch("matrix must be (src_planes, dest_planes)"))
            m = convert(Vector{$KT}, vec(permutedims(convert(Matrix{$KT}, matrix))))
            prev = preBias === nothing ? $KT[] : convert(Vector{$KT}, collect(preBias))
            postv = postBias === nothing ? Int32[] : convert(Vector{Int32}, collect(postBias))
            srefs = [Ref(vimage_buffer(s)) for s in srcs]; drefs = [Ref(vimage_buffer(d)) for d in dests]
            GC.@preserve srcs dests m prev postv srefs drefs begin
                sptrs = Ptr{vImage_Buffer}[Base.unsafe_convert(Ptr{vImage_Buffer}, r) for r in srefs]
                dptrs = Ptr{vImage_Buffer}[Base.unsafe_convert(Ptr{vImage_Buffer}, r) for r in drefs]
                prep = isempty(prev) ? Ptr{$KT}(C_NULL) : pointer(prev)
                postp = isempty(postv) ? Ptr{Int32}(C_NULL) : pointer(postv)
                GC.@preserve sptrs dptrs begin
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{Ptr{vImage_Buffer}}, Ptr{Ptr{vImage_Buffer}}, UInt32, UInt32,
                         Ptr{$KT}, Int32, Ptr{$KT}, Ptr{Int32}, vImage_Flags),
                        pointer(sptrs), pointer(dptrs), UInt32(sp), UInt32(dp),
                        pointer(m), Int32(divisor), prep, postp, vImage_Flags(flags))
                end
            end
            _check(err, $(String(sym))); return dests
        end
    else
        @eval function $fn(srcs::AbstractVector{<:AbstractMatrix{$T}}, dests::AbstractVector{<:AbstractMatrix{$T}},
                           matrix::AbstractMatrix{<:Real}; preBias = nothing, postBias = nothing,
                           flags::Integer = kvImageNoFlags)
            sp = length(srcs); dp = length(dests)
            size(matrix) == (sp, dp) || throw(DimensionMismatch("matrix must be (src_planes, dest_planes)"))
            m = convert(Vector{Cfloat}, vec(permutedims(convert(Matrix{Cfloat}, matrix))))
            prev = preBias === nothing ? Cfloat[] : convert(Vector{Cfloat}, collect(preBias))
            postv = postBias === nothing ? Cfloat[] : convert(Vector{Cfloat}, collect(postBias))
            srefs = [Ref(vimage_buffer(s)) for s in srcs]; drefs = [Ref(vimage_buffer(d)) for d in dests]
            GC.@preserve srcs dests m prev postv srefs drefs begin
                sptrs = Ptr{vImage_Buffer}[Base.unsafe_convert(Ptr{vImage_Buffer}, r) for r in srefs]
                dptrs = Ptr{vImage_Buffer}[Base.unsafe_convert(Ptr{vImage_Buffer}, r) for r in drefs]
                prep = isempty(prev) ? Ptr{Cfloat}(C_NULL) : pointer(prev)
                postp = isempty(postv) ? Ptr{Cfloat}(C_NULL) : pointer(postv)
                GC.@preserve sptrs dptrs begin
                    err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                        (Ptr{Ptr{vImage_Buffer}}, Ptr{Ptr{vImage_Buffer}}, UInt32, UInt32,
                         Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, vImage_Flags),
                        pointer(sptrs), pointer(dptrs), UInt32(sp), UInt32(dp),
                        pointer(m), prep, postp, vImage_Flags(flags))
                end
            end
            _check(err, $(String(sym))); return dests
        end
    end
end

"""
    matrixMultiply_Planar8(srcs, dests, matrix; divisor=1, preBias, postBias, flags) -> dests
    matrixMultiply_PlanarF(srcs, dests, matrix; preBias, postBias, flags) -> dests

Multiply the `M` source planes by an `M×N` `matrix` to produce the `N` destination
planes: `out[j] = Σ_i (in[i]+preBias[i])·matrix[i,j]` (then `+postBias[j]`, and `/divisor`
for the integer forms `Planar8`/`Planar16S`). `srcs`/`dests` are vectors of equally-sized
planar matrices. Writes into `dests`.
"""
matrixMultiply_Planar8

# =====================================================================================
# ADDITIONAL ALPHA / COMPOSITING
# =====================================================================================

# Premultiplied planar blend: (srcTop, srcTopAlpha, srcBottom, dest, flags)
for sfx in (:Planar8, :PlanarF)
    sym  = Symbol("vImagePremultipliedAlphaBlend_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractMatrix, srcTop::AbstractMatrix, srcTopAlpha::AbstractMatrix,
                       srcBottom::AbstractMatrix; flags::Integer = kvImageNoFlags)
            GC.@preserve srcTop srcTopAlpha srcBottom dest begin
                t = Ref(vimage_buffer(srcTop)); ta = Ref(vimage_buffer(srcTopAlpha))
                b = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    t, ta, b, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
        $(_jlname(sym))(srcTop::AbstractMatrix, srcTopAlpha::AbstractMatrix, srcBottom::AbstractMatrix; kw...) =
            $bang(similar(srcTop), srcTop, srcTopAlpha, srcBottom; kw...)
    end
end

# Nonpremultiplied-to-premultiplied planar blend: same 4-buffer shape
for sfx in (:Planar8, :PlanarF)
    sym  = Symbol("vImageAlphaBlend_NonpremultipliedToPremultiplied_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractMatrix, srcTop::AbstractMatrix, srcTopAlpha::AbstractMatrix,
                       srcBottom::AbstractMatrix; flags::Integer = kvImageNoFlags)
            GC.@preserve srcTop srcTopAlpha srcBottom dest begin
                t = Ref(vimage_buffer(srcTop)); ta = Ref(vimage_buffer(srcTopAlpha))
                b = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    t, ta, b, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
        $(_jlname(sym))(srcTop::AbstractMatrix, srcTopAlpha::AbstractMatrix, srcBottom::AbstractMatrix; kw...) =
            $bang(similar(srcTop), srcTop, srcTopAlpha, srcBottom; kw...)
    end
end

# Const-alpha premultiplied blend. Planar: (srcTop, constAlpha, srcTopAlpha, srcBottom, dest);
# interleaved: (srcTop, constAlpha, srcBottom, dest).
for (sfx, T, AT) in ((:Planar8, UInt8, UInt8), (:PlanarF, Float32, Cfloat))
    sym  = Symbol("vImagePremultipliedConstAlphaBlend_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractMatrix{$T}, srcTop::AbstractMatrix{$T}, constAlpha,
                       srcTopAlpha::AbstractMatrix{$T}, srcBottom::AbstractMatrix{$T}; flags::Integer = kvImageNoFlags)
            GC.@preserve srcTop srcTopAlpha srcBottom dest begin
                t = Ref(vimage_buffer(srcTop)); ta = Ref(vimage_buffer(srcTopAlpha))
                b = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, $AT, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    t, $AT(constAlpha), ta, b, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
        $(_jlname(sym))(srcTop::AbstractMatrix{$T}, constAlpha, srcTopAlpha::AbstractMatrix{$T}, srcBottom::AbstractMatrix{$T}; kw...) =
            $bang(similar(srcTop), srcTop, constAlpha, srcTopAlpha, srcBottom; kw...)
    end
end
for (sfx, T, AT) in ((:ARGB8888, UInt8, UInt8), (:ARGBFFFF, Float32, Cfloat))
    sym  = Symbol("vImagePremultipliedConstAlphaBlend_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractArray{$T,3}, srcTop::AbstractArray{$T,3}, constAlpha,
                       srcBottom::AbstractArray{$T,3}; flags::Integer = kvImageNoFlags)
            GC.@preserve srcTop srcBottom dest begin
                t = Ref(vimage_buffer(srcTop)); b = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, $AT, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                    t, $AT(constAlpha), b, db, vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
        $(_jlname(sym))(srcTop::AbstractArray{$T,3}, constAlpha, srcBottom::AbstractArray{$T,3}; kw...) =
            $bang(similar(srcTop), srcTop, constAlpha, srcBottom; kw...)
    end
end

"""
    premultipliedAlphaBlend_Planar8(srcTop, srcTopAlpha, srcBottom; flags) -> plane
    premultipliedConstAlphaBlend_ARGB8888(srcTop, constAlpha, srcBottom; flags) -> array

Additional alpha-compositing operators. `premultipliedAlphaBlend_Planar8/PlanarF` blend a
premultiplied planar top (with its alpha plane) over a premultiplied bottom;
`alphaBlend_NonpremultipliedToPremultiplied_Planar8/PlanarF` do the same from a
non-premultiplied top; `premultipliedConstAlphaBlend_*` scale the top's alpha by a scalar
`constAlpha` (planar forms take a separate `srcTopAlpha` plane, interleaved
`ARGB8888/FFFF` forms use the in-pixel alpha). Allocating and `!` variants.
"""
premultipliedAlphaBlend_Planar8

# AlphaBlend_PlanarF: (srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha, dest, flags)
function alphaBlend_PlanarF!(dest::AbstractMatrix{Float32}, srcTop::AbstractMatrix{Float32}, srcTopAlpha::AbstractMatrix{Float32},
                            srcBottom::AbstractMatrix{Float32}, srcBottomAlpha::AbstractMatrix{Float32},
                            alpha::AbstractMatrix{Float32}; flags::Integer = kvImageNoFlags)
    GC.@preserve srcTop srcTopAlpha srcBottom srcBottomAlpha alpha dest begin
        t = Ref(vimage_buffer(srcTop)); ta = Ref(vimage_buffer(srcTopAlpha))
        b = Ref(vimage_buffer(srcBottom)); ba = Ref(vimage_buffer(srcBottomAlpha))
        al = Ref(vimage_buffer(alpha)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageAlphaBlend_PlanarF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer},
             Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            t, ta, b, ba, al, db, vImage_Flags(flags))
    end
    _check(err, "vImageAlphaBlend_PlanarF"); return dest
end
"""
    alphaBlend_PlanarF(srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha; flags) -> plane

Float non-premultiplied planar alpha blend (float analogue of
[`alphaBlend_Planar8`](@ref)). Allocating and `!` variants.
"""
alphaBlend_PlanarF(srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha; kw...) =
    alphaBlend_PlanarF!(similar(srcTop), srcTop, srcTopAlpha, srcBottom, srcBottomAlpha, alpha; kw...)

# PremultipliedAlphaBlendWithPermute: (srcTop, srcBottom, dest, permuteMap[4], makeDestAlphaOpaque, flags)
for sfx in (:ARGB8888, :RGBA8888)
    sym  = Symbol("vImagePremultipliedAlphaBlendWithPermute_", sfx)
    bang = Symbol(_jlname(sym), "!")
    @eval begin
        function $bang(dest::AbstractArray{UInt8,3}, srcTop::AbstractArray{UInt8,3}, srcBottom::AbstractArray{UInt8,3};
                       permuteMap = (0, 1, 2, 3), makeDestAlphaOpaque::Bool = false, flags::Integer = kvImageNoFlags)
            pm = convert(Vector{UInt8}, collect(permuteMap))
            GC.@preserve srcTop srcBottom dest pm begin
                t = Ref(vimage_buffer(srcTop)); b = Ref(vimage_buffer(srcBottom)); db = Ref(vimage_buffer(dest))
                err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                    (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{UInt8}, Bool, vImage_Flags),
                    t, b, db, pointer(pm), makeDestAlphaOpaque, vImage_Flags(flags))
            end
            _check(err, $(String(sym))); return dest
        end
        $(_jlname(sym))(srcTop::AbstractArray{UInt8,3}, srcBottom::AbstractArray{UInt8,3}; kw...) =
            $bang(similar(srcTop), srcTop, srcBottom; kw...)
    end
end
"""
    premultipliedAlphaBlendWithPermute_ARGB8888(srcTop, srcBottom; permuteMap=(0,1,2,3),
                                                makeDestAlphaOpaque=false, flags) -> array

Source-over premultiplied blend that first permutes the top image's channels by
`permuteMap`; set `makeDestAlphaOpaque` to force the result alpha to opaque. Variants for
`ARGB8888` and `RGBA8888`. Allocating and `!` variants.
"""
premultipliedAlphaBlendWithPermute_ARGB8888

# =====================================================================================
# HISTOGRAM SPECIFICATION (float + interleaved)
# =====================================================================================

function histogramSpecification_PlanarF!(dest::AbstractMatrix{Float32}, src::AbstractMatrix{Float32},
        desiredHistogram::AbstractVector; entries::Integer = length(desiredHistogram),
        minVal::Real = 0f0, maxVal::Real = 1f0, flags::Integer = kvImageNoFlags)
    dh = convert(Vector{Csize_t}, collect(desiredHistogram))
    GC.@preserve src dest dh begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageHistogramSpecification_PlanarF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{Csize_t}, Cuint, Cfloat, Cfloat, vImage_Flags),
            sb, db, C_NULL, pointer(dh), Cuint(entries), Cfloat(minVal), Cfloat(maxVal), vImage_Flags(flags))
    end
    _check(err, "vImageHistogramSpecification_PlanarF"); return dest
end
histogramSpecification_PlanarF(src::AbstractMatrix{Float32}, dh::AbstractVector; kw...) =
    histogramSpecification_PlanarF!(similar(src), src, dh; kw...)

function histogramSpecification_ARGB8888!(dest::AbstractArray{UInt8,3}, src::AbstractArray{UInt8,3},
        desiredHistogram; flags::Integer = kvImageNoFlags)
    hs = [convert(Vector{Csize_t}, collect(h)) for h in desiredHistogram]
    length(hs) == 4 || throw(DimensionMismatch("need 4 desired histograms"))
    GC.@preserve src dest hs begin
        hptrs = Ptr{Csize_t}[pointer(h) for h in hs]
        GC.@preserve hptrs begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall((:vImageHistogramSpecification_ARGB8888, vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Ptr{Csize_t}}, vImage_Flags),
                sb, db, pointer(hptrs), vImage_Flags(flags))
        end
    end
    _check(err, "vImageHistogramSpecification_ARGB8888"); return dest
end
histogramSpecification_ARGB8888(src::AbstractArray{UInt8,3}, dh; kw...) =
    histogramSpecification_ARGB8888!(similar(src), src, dh; kw...)

function histogramSpecification_ARGBFFFF!(dest::AbstractArray{Float32,3}, src::AbstractArray{Float32,3},
        desiredHistogram; entries::Integer = length(first(desiredHistogram)),
        minVal::Real = 0f0, maxVal::Real = 1f0, flags::Integer = kvImageNoFlags)
    hs = [convert(Vector{Csize_t}, collect(h)) for h in desiredHistogram]
    length(hs) == 4 || throw(DimensionMismatch("need 4 desired histograms"))
    GC.@preserve src dest hs begin
        hptrs = Ptr{Csize_t}[pointer(h) for h in hs]
        GC.@preserve hptrs begin
            sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
            err = ccall((:vImageHistogramSpecification_ARGBFFFF, vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{Ptr{Csize_t}}, Cuint, Cfloat, Cfloat, vImage_Flags),
                sb, db, C_NULL, pointer(hptrs), Cuint(entries), Cfloat(minVal), Cfloat(maxVal), vImage_Flags(flags))
        end
    end
    _check(err, "vImageHistogramSpecification_ARGBFFFF"); return dest
end
histogramSpecification_ARGBFFFF(src::AbstractArray{Float32,3}, dh; kw...) =
    histogramSpecification_ARGBFFFF!(similar(src), src, dh; kw...)

"""
    histogramSpecification_PlanarF(src, desiredHistogram; entries, minVal=0, maxVal=1, flags) -> plane
    histogramSpecification_ARGB8888(src, desiredHistogram; flags) -> array

Remap an image so its histogram matches `desiredHistogram`. `PlanarF`/`ARGBFFFF` take the
histogram bin count `entries` and value range `[minVal, maxVal]` (the `ARGB*` forms take a
length-4 collection of per-channel histograms). Allocating and `!` variants; see also
[`histogramSpecification_Planar8`](@ref).
"""
histogramSpecification_PlanarF

function endsInContrastStretch_PlanarF!(dest::AbstractMatrix{Float32}, src::AbstractMatrix{Float32};
        percentLow::Integer = 0, percentHigh::Integer = 0, entries::Integer = 4096,
        minVal::Real = 0f0, maxVal::Real = 1f0, flags::Integer = kvImageNoFlags)
    GC.@preserve src dest begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageEndsInContrastStretch_PlanarF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Cuint, Cuint, Cuint, Cfloat, Cfloat, vImage_Flags),
            sb, db, C_NULL, Cuint(percentLow), Cuint(percentHigh), Cuint(entries), Cfloat(minVal), Cfloat(maxVal), vImage_Flags(flags))
    end
    _check(err, "vImageEndsInContrastStretch_PlanarF"); return dest
end
endsInContrastStretch_PlanarF(src::AbstractMatrix{Float32}; kw...) =
    endsInContrastStretch_PlanarF!(similar(src), src; kw...)

function endsInContrastStretch_ARGBFFFF!(dest::AbstractArray{Float32,3}, src::AbstractArray{Float32,3};
        percentLow = (0, 0, 0, 0), percentHigh = (0, 0, 0, 0), entries::Integer = 4096,
        minVal::Real = 0f0, maxVal::Real = 1f0, flags::Integer = kvImageNoFlags)
    pl = convert(Vector{Cuint}, collect(percentLow)); ph = convert(Vector{Cuint}, collect(percentHigh))
    GC.@preserve src dest pl ph begin
        sb = Ref(vimage_buffer(src)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageEndsInContrastStretch_ARGBFFFF, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cvoid}, Ptr{Cuint}, Ptr{Cuint}, Cuint, Cfloat, Cfloat, vImage_Flags),
            sb, db, C_NULL, pointer(pl), pointer(ph), Cuint(entries), Cfloat(minVal), Cfloat(maxVal), vImage_Flags(flags))
    end
    _check(err, "vImageEndsInContrastStretch_ARGBFFFF"); return dest
end
endsInContrastStretch_ARGBFFFF(src::AbstractArray{Float32,3}; kw...) =
    endsInContrastStretch_ARGBFFFF!(similar(src), src; kw...)

"""
    endsInContrastStretch_PlanarF(src; percentLow=0, percentHigh=0, entries=4096, minVal, maxVal, flags) -> plane

Float ends-in contrast stretch (clip `percentLow`%/`percentHigh`% of the darkest/brightest
pixels, then rescale). `endsInContrastStretch_ARGBFFFF` takes length-4 percentages. See the
integer [`endsInContrastStretch_Planar8`](@ref). Allocating and `!` variants.
"""
endsInContrastStretch_PlanarF

# =====================================================================================
# ADDITIONAL PLANE ⇄ INTERLEAVED CONVERSIONS
# =====================================================================================

# 4-channel interleaved -> 3 planes (drop the X channel): (src, r, g, b, flags) with
# XRGB order (dest planes r,g,b) and BGRX order (dest planes b,g,r).
for (sym, T) in ((:vImageConvert_XRGB8888ToPlanar8, UInt8), (:vImageConvert_BGRX8888ToPlanar8, UInt8),
                 (:vImageConvert_XRGBFFFFToPlanarF, Float32), (:vImageConvert_BGRXFFFFToPlanarF, Float32))
    fn = _jlname(sym)
    @eval function $fn(src::AbstractArray{$T,3}, c1::AbstractMatrix{$T}, c2::AbstractMatrix{$T},
                       c3::AbstractMatrix{$T}; flags::Integer = kvImageNoFlags)
        GC.@preserve src c1 c2 c3 begin
            sb = Ref(vimage_buffer(src)); p1 = Ref(vimage_buffer(c1)); p2 = Ref(vimage_buffer(c2)); p3 = Ref(vimage_buffer(c3))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
                sb, p1, p2, p3, vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return (c1, c2, c3)
    end
end
"""
    convert_XRGB8888ToPlanar8(src, red, green, blue) -> (red, green, blue)
    convert_BGRX8888ToPlanar8(src, blue, green, red) -> (blue, green, red)

Deinterleave a 4-channel image into three planes, discarding the padding (X) channel.
`XRGB*` writes planes in `r,g,b` order, `BGRX*` in `b,g,r` order (`Planar8`/`PlanarF`).
"""
convert_XRGB8888ToPlanar8

# 4 planes (A,R,G,B) -> 4-channel float, mapping [0,255]/[0,1] to [minFloat,maxFloat] per channel.
for (sym, DT, T) in ((:vImageConvert_Planar8ToARGBFFFF, Float32, UInt8),
                     (:vImageConvert_PlanarFToARGB8888, UInt8, Float32))
    bang = Symbol(_jlname(sym), "!")
    @eval function $bang(dest::AbstractArray{$DT,3}, alpha::AbstractMatrix{$T}, red::AbstractMatrix{$T},
                         green::AbstractMatrix{$T}, blue::AbstractMatrix{$T};
                         maxFloat = (1f0, 1f0, 1f0, 1f0), minFloat = (0f0, 0f0, 0f0, 0f0), flags::Integer = kvImageNoFlags)
        mx = convert(Vector{Cfloat}, collect(maxFloat)); mn = convert(Vector{Cfloat}, collect(minFloat))
        GC.@preserve dest alpha red green blue mx mn begin
            ab = Ref(vimage_buffer(alpha)); rb = Ref(vimage_buffer(red)); gb = Ref(vimage_buffer(green)); bb = Ref(vimage_buffer(blue)); db = Ref(vimage_buffer(dest))
            err = ccall(($(QuoteNode(sym)), vimage_lib), vImage_Error,
                (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{Cfloat}, Ptr{Cfloat}, vImage_Flags),
                ab, rb, gb, bb, db, pointer(mx), pointer(mn), vImage_Flags(flags))
        end
        _check(err, $(String(sym))); return dest
    end
end
"""
    convert_Planar8ToARGBFFFF!(dest, alpha, red, green, blue; maxFloat, minFloat, flags) -> dest
    convert_PlanarFToARGB8888!(dest, alpha, red, green, blue; maxFloat, minFloat, flags) -> dest

Interleave four planes `(alpha, red, green, blue)` into a 4-channel image, rescaling each
channel between `minFloat` and `maxFloat` (length-4 collections). `Planar8ToARGBFFFF`
produces float output from 8-bit planes; `PlanarFToARGB8888` produces 8-bit output from
float planes. `!`-only.
"""
convert_Planar8ToARGBFFFF!

# Planar16Q12 (Q4.12 fixed-point) <-> 8-bit interleaved, via separate planes.
function convert_Planar16Q12toARGB8888!(dest::AbstractArray{UInt8,3}, alpha::AbstractMatrix{Int16},
        red::AbstractMatrix{Int16}, green::AbstractMatrix{Int16}, blue::AbstractMatrix{Int16}; flags::Integer = kvImageNoFlags)
    GC.@preserve dest alpha red green blue begin
        ab = Ref(vimage_buffer(alpha)); rb = Ref(vimage_buffer(red)); gb = Ref(vimage_buffer(green)); bb = Ref(vimage_buffer(blue)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageConvert_Planar16Q12toARGB8888, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            ab, rb, gb, bb, db, vImage_Flags(flags))
    end
    _check(err, "vImageConvert_Planar16Q12toARGB8888"); return dest
end
function convert_Planar16Q12toRGB888!(dest::AbstractArray{UInt8,3}, red::AbstractMatrix{Int16},
        green::AbstractMatrix{Int16}, blue::AbstractMatrix{Int16}; flags::Integer = kvImageNoFlags)
    GC.@preserve dest red green blue begin
        rb = Ref(vimage_buffer(red)); gb = Ref(vimage_buffer(green)); bb = Ref(vimage_buffer(blue)); db = Ref(vimage_buffer(dest))
        err = ccall((:vImageConvert_Planar16Q12toRGB888, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            rb, gb, bb, db, vImage_Flags(flags))
    end
    _check(err, "vImageConvert_Planar16Q12toRGB888"); return dest
end
function convert_ARGB8888toPlanar16Q12!(alpha::AbstractMatrix{Int16}, red::AbstractMatrix{Int16},
        green::AbstractMatrix{Int16}, blue::AbstractMatrix{Int16}, src::AbstractArray{UInt8,3}; flags::Integer = kvImageNoFlags)
    GC.@preserve src alpha red green blue begin
        sb = Ref(vimage_buffer(src)); ab = Ref(vimage_buffer(alpha)); rb = Ref(vimage_buffer(red)); gb = Ref(vimage_buffer(green)); bb = Ref(vimage_buffer(blue))
        err = ccall((:vImageConvert_ARGB8888toPlanar16Q12, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            sb, ab, rb, gb, bb, vImage_Flags(flags))
    end
    _check(err, "vImageConvert_ARGB8888toPlanar16Q12"); return (alpha, red, green, blue)
end
function convert_RGB888toPlanar16Q12!(red::AbstractMatrix{Int16}, green::AbstractMatrix{Int16},
        blue::AbstractMatrix{Int16}, src::AbstractArray{UInt8,3}; flags::Integer = kvImageNoFlags)
    GC.@preserve src red green blue begin
        sb = Ref(vimage_buffer(src)); rb = Ref(vimage_buffer(red)); gb = Ref(vimage_buffer(green)); bb = Ref(vimage_buffer(blue))
        err = ccall((:vImageConvert_RGB888toPlanar16Q12, vimage_lib), vImage_Error,
            (Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, Ptr{vImage_Buffer}, vImage_Flags),
            sb, rb, gb, bb, vImage_Flags(flags))
    end
    _check(err, "vImageConvert_RGB888toPlanar16Q12"); return (red, green, blue)
end
"""
    convert_Planar16Q12toARGB8888!(dest, alpha, red, green, blue) -> dest
    convert_ARGB8888toPlanar16Q12!(alpha, red, green, blue, src) -> (alpha,red,green,blue)

Convert between 8-bit interleaved images and `Planar16Q12` (signed Q4.12 fixed-point, where
`4096 == 1.0`) held as separate `Int16` planes. `!`-only; the `RGB888` variants omit the
alpha plane.
"""
convert_Planar16Q12toARGB8888!
