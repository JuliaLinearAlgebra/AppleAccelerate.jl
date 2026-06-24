## complexarray.jl — Complex-valued vDSP operations using split-complex format ##

# ============================================================
# Split-complex structs (used here and by dsp.jl for FFT)
# ============================================================
struct DSPDoubleSplitComplex
    realp::Ptr{Float64}
    imagp::Ptr{Float64}
end

DSPDoubleSplitComplex(realp::Vector{Float64}, imagp::Vector{Float64}) = DSPDoubleSplitComplex(pointer(realp), pointer(imagp))

struct DSPSplitComplex
    realp::Ptr{Float32}
    imagp::Ptr{Float32}
end

DSPSplitComplex(realp::Vector{Float32}, imagp::Vector{Float32}) = DSPSplitComplex(pointer(realp), pointer(imagp))

# ------------------------------------------------------------
# Interleaved ⇄ split-complex bridge
#
# A Julia `Vector{Complex{T}}` is stored interleaved as [re, im, re, im, …].
# A DSPSplitComplex whose `realp` points at the base and whose `imagp` points
# one real element later describes exactly that layout *if* the vDSP routine is
# given a stride of 2 (`_CSTRIDE`): element i is then at realp[2i] / imagp[2i],
# i.e. base[2i] / base[2i+1]. This lets the split-complex vDSP routines read and
# write Julia complex buffers directly — no deinterleave/reinterleave copies, and
# the mutating (`!`) variants become genuinely in-place and allocation-free.
#
# The returned struct only borrows the pointer; callers MUST keep the backing
# array alive with `GC.@preserve` for the duration of the ccall.
# ------------------------------------------------------------
const _CSTRIDE = 2

@inline _split_view(::Type{DSPSplitComplex}, p::Ptr) =
    DSPSplitComplex(Ptr{Float32}(p), Ptr{Float32}(p) + sizeof(Float32))
@inline _split_view(::Type{DSPDoubleSplitComplex}, p::Ptr) =
    DSPDoubleSplitComplex(Ptr{Float64}(p), Ptr{Float64}(p) + sizeof(Float64))

# ============================================================
# Complex → Complex unary operations (operate on interleaved storage)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))

    # vneg: negate complex vector
    @eval begin
        function vneg!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}})
            length(X) == length(result) || throw(DimensionMismatch("vneg!: X and result must have equal lengths"))
            n = length(X)
            GC.@preserve X result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvneg", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function vneg(X::Vector{Complex{$T}})
            result = similar(X)
            vneg!(result, X)
        end
    end

    # vconj: conjugate complex vector
    @eval begin
        function vconj!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}})
            length(X) == length(result) || throw(DimensionMismatch("vconj!: X and result must have equal lengths"))
            n = length(X)
            GC.@preserve X result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvconj", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function vconj(X::Vector{Complex{$T}})
            result = similar(X)
            vconj!(result, X)
        end
    end

    # vcopy: copy complex vector (via vDSP_zvmov)
    @eval begin
        function vcopy(X::Vector{Complex{$T}})
            n = length(X)
            result = Vector{Complex{$T}}(undef, n)
            GC.@preserve X result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvmov", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
    end
end

@doc """
    vconj(X::Vector{Complex{T}}) -> Vector{Complex{T}}
    vconj!(result, X)

Complex conjugate: `result[i] = conj(X[i])`.
Wraps [`vDSP_zvconj`](https://developer.apple.com/documentation/accelerate/vdsp_zvconj).
""" vconj

@doc """
    vcopy(X::Vector{Complex{T}}) -> Vector{Complex{T}}

Copy complex vector via split-complex move.
Wraps [`vDSP_zvmov`](https://developer.apple.com/documentation/accelerate/vdsp_zvmov).
""" vcopy

# ============================================================
# Complex → Complex binary operations (operate on interleaved storage)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))

    # vmul: element-wise complex multiply (conjugate flag = +1 for normal multiply)
    @eval begin
        function vmul!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}}, Y::Vector{Complex{$T}})
            length(X) == length(Y) == length(result) || throw(DimensionMismatch("vmul!: X, Y, and result must have equal lengths"))
            n = length(X)
            GC.@preserve X Y result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ysplit = _split_view($DSPSplit, pointer(Y))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64, Cint),
                      xsplit, _CSTRIDE, ysplit, _CSTRIDE, osplit, _CSTRIDE, n, Cint(1))
            end
            return result
        end
        function vmul(X::Vector{Complex{$T}}, Y::Vector{Complex{$T}})
            result = similar(X)
            vmul!(result, X, Y)
        end
    end

    # vdiv: element-wise complex divide
    # Note: vDSP_zvdiv computes B/A, so we swap: pass Y as A and X as B to get X/Y
    @eval begin
        function vdiv!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}}, Y::Vector{Complex{$T}})
            length(X) == length(Y) == length(result) || throw(DimensionMismatch("vdiv!: X, Y, and result must have equal lengths"))
            n = length(X)
            GC.@preserve X Y result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ysplit = _split_view($DSPSplit, pointer(Y))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvdiv", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      ysplit, _CSTRIDE, xsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function vdiv(X::Vector{Complex{$T}}, Y::Vector{Complex{$T}})
            result = similar(X)
            vdiv!(result, X, Y)
        end
    end

    # vsmul: complex vector * complex scalar
    # vDSP_zvzsml: A * B → C where B is a single split-complex element
    @eval begin
        function vsmul!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}}, c::Complex{$T})
            length(X) == length(result) || throw(DimensionMismatch("vsmul!: X and result must have equal lengths"))
            n = length(X)
            c2 = $T[real(c), imag(c)]
            GC.@preserve X result c2 begin
                xsplit = _split_view($DSPSplit, pointer(X))
                csplit = _split_view($DSPSplit, pointer(c2))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvzsml", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, _CSTRIDE, csplit, osplit, _CSTRIDE, n)
            end
            return result
        end
        function vsmul(X::Vector{Complex{$T}}, c::Complex{$T})
            result = similar(X)
            vsmul!(result, X, c)
        end
    end
end

# ============================================================
# Complex → Real operations (interleaved complex in, real out)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))

    # vabs: complex absolute value (modulus)
    @eval begin
        function vabs!(result::Vector{$T}, X::Vector{Complex{$T}})
            length(X) == length(result) || throw(DimensionMismatch("vabs!: X and result must have equal lengths"))
            n = length(X)
            GC.@preserve X result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ccall(($(string("vDSP_zvabs", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, _CSTRIDE, result, 1, n)
            end
            return result
        end
        function vabs(X::Vector{Complex{$T}})
            result = Vector{$T}(undef, length(X))
            vabs!(result, X)
        end
    end

    # vphase: complex phase (angle)
    @eval begin
        function vphase!(result::Vector{$T}, X::Vector{Complex{$T}})
            length(X) == length(result) || throw(DimensionMismatch("vphase!: X and result must have equal lengths"))
            n = length(X)
            GC.@preserve X result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ccall(($(string("vDSP_zvphas", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, _CSTRIDE, result, 1, n)
            end
            return result
        end
        function vphase(X::Vector{Complex{$T}})
            result = Vector{$T}(undef, length(X))
            vphase!(result, X)
        end
    end

    # vmags: squared magnitude (abs2)
    @eval begin
        function vmags!(result::Vector{$T}, X::Vector{Complex{$T}})
            length(X) == length(result) || throw(DimensionMismatch("vmags!: X and result must have equal lengths"))
            n = length(X)
            GC.@preserve X result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ccall(($(string("vDSP_zvmags", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, _CSTRIDE, result, 1, n)
            end
            return result
        end
        function vmags(X::Vector{Complex{$T}})
            result = Vector{$T}(undef, length(X))
            vmags!(result, X)
        end
    end

    # vmagsa: squared magnitude + accumulate: abs2(X) + B
    @eval begin
        function vmagsa!(result::Vector{$T}, X::Vector{Complex{$T}}, B::Vector{$T})
            length(X) == length(B) == length(result) || throw(DimensionMismatch("vmagsa!: X, B, and result must have equal lengths"))
            n = length(X)
            GC.@preserve X B result begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ccall(($(string("vDSP_zvmgsa", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, _CSTRIDE, B, 1, result, 1, n)
            end
            return result
        end
        function vmagsa(X::Vector{Complex{$T}}, B::Vector{$T})
            result = Vector{$T}(undef, length(X))
            vmagsa!(result, X, B)
        end
    end
end

@doc """
    vphase(X::Vector{Complex{T}}) -> Vector{T}
    vphase!(result, X)

Complex phase (angle): `result[i] = atan(imag(X[i]), real(X[i]))`.
Wraps [`vDSP_zvphas`](https://developer.apple.com/documentation/accelerate/vdsp_zvphas).
""" vphase

@doc """
    vmags(X::Vector{Complex{T}}) -> Vector{T}
    vmags!(result, X)

Squared magnitude: `result[i] = abs2(X[i]) = real(X[i])^2 + imag(X[i])^2`.
Wraps [`vDSP_zvmags`](https://developer.apple.com/documentation/accelerate/vdsp_zvmags).
""" vmags

@doc """
    vmagsa(X::Vector{Complex{T}}, B::Vector{T}) -> Vector{T}
    vmagsa!(result, X, B)

Squared magnitude and accumulate: `result[i] = abs2(X[i]) + B[i]`.
Wraps [`vDSP_zvmgsa`](https://developer.apple.com/documentation/accelerate/vdsp_zvmgsa).
""" vmagsa

# ============================================================
# Complex dot product (returns scalar)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))
    @eval begin
        function dot(X::Vector{Complex{$T}}, Y::Vector{Complex{$T}})
            length(X) == length(Y) || throw(DimensionMismatch("dot: X and Y must have equal lengths"))
            n = length(X)
            out = $T[0, 0]                 # interleaved [re, im] scalar result
            GC.@preserve X Y out begin
                xsplit = _split_view($DSPSplit, pointer(X))
                ysplit = _split_view($DSPSplit, pointer(Y))
                osplit = _split_view($DSPSplit, pointer(out))
                ccall(($(string("vDSP_zdotpr", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, UInt64),
                      xsplit, _CSTRIDE, ysplit, _CSTRIDE, osplit, n)
            end
            return complex(out[1], out[2])
        end
    end
end

# ============================================================
# Coordinate conversion (interleaved format)
# ============================================================
for (T, suff) in ((Float32, ""), (Float64, "D"))

    # polar: Cartesian (re, im interleaved) → (magnitude, angle) interleaved
    # vDSP_polar takes interleaved input [re, im, re, im, ...] and writes [mag, angle, mag, angle, ...]
    @eval begin
        function polar(X::Vector{Complex{$T}})
            n = length(X)
            interleaved_out = Vector{$T}(undef, 2 * n)
            GC.@preserve X interleaved_out begin
                ccall(($(string("vDSP_polar", suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      pointer(X), 2, interleaved_out, 2, n)
            end
            magnitudes = interleaved_out[1:2:end]
            angles = interleaved_out[2:2:end]
            return (magnitudes, angles)
        end
    end

    # rect: (magnitude, angle) → Cartesian complex
    # vDSP_rect takes interleaved [mag, angle, mag, angle, ...] and writes [re, im, re, im, ...]
    @eval begin
        function rect(magnitudes::Vector{$T}, angles::Vector{$T})
            length(magnitudes) == length(angles) || throw(DimensionMismatch("rect: magnitudes and angles must have equal lengths"))
            n = length(magnitudes)
            interleaved_in = Vector{$T}(undef, 2 * n)
            interleaved_in[1:2:end] .= magnitudes
            interleaved_in[2:2:end] .= angles
            result = Vector{Complex{$T}}(undef, n)
            GC.@preserve interleaved_in result begin
                ccall(($(string("vDSP_rect", suff)), libacc), Cvoid,
                      (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      interleaved_in, 2, pointer(result), 2, n)
            end
            return result
        end
    end
end

@doc """
    polar(X::Vector{Complex{T}}) -> (magnitudes::Vector{T}, angles::Vector{T})

Convert complex Cartesian coordinates to polar form.
Wraps [`vDSP_polar`](https://developer.apple.com/documentation/accelerate/vdsp_polar).
""" polar

@doc """
    rect(magnitudes::Vector{T}, angles::Vector{T}) -> Vector{Complex{T}}

Convert polar coordinates to complex Cartesian form.
Wraps [`vDSP_rect`](https://developer.apple.com/documentation/accelerate/vdsp_rect).
""" rect

# ============================================================
# Batch 5: Additional Complex Vector Operations
# ============================================================

for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                             (Float64, "D", :DSPDoubleSplitComplex))

    # --- Complex-complex binary ops ---

    # zvadd: C = A + B
    @eval begin
        function zvadd!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zvadd!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvadd", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, bsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zvadd(A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            result = similar(A)
            zvadd!(result, A, B)
        end
    end

    # zvsub: C = A - B  (NOTE: vDSP_zvsub computes __A - __B, pass A as first, B as second)
    @eval begin
        function zvsub!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zvsub!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvsub", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, bsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zvsub(A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            result = similar(A)
            zvsub!(result, A, B)
        end
    end

    # zvcmul: C = conj(A) * B
    @eval begin
        function zvcmul!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zvcmul!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvcmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, bsplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zvcmul(A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            result = similar(A)
            zvcmul!(result, A, B)
        end
    end

    # --- Complex-real binary ops ---

    # zrvmul: C = A * B (complex * real)
    @eval begin
        function zrvmul!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{$T})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zrvmul!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zrvmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, B, 1, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zrvmul(A::Vector{Complex{$T}}, B::Vector{$T})
            result = similar(A)
            zrvmul!(result, A, B)
        end
    end

    # zrvdiv: C = A / B (complex / real) — NOTE: vDSP_zrvdiv computes B/A, swap
    @eval begin
        function zrvdiv!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{$T})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zrvdiv!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zrvdiv", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, B, 1, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zrvdiv(A::Vector{Complex{$T}}, B::Vector{$T})
            result = similar(A)
            zrvdiv!(result, A, B)
        end
    end

    # zrvadd: C = A + B (add real to complex, adds to real part)
    @eval begin
        function zrvadd!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{$T})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zrvadd!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zrvadd", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, B, 1, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zrvadd(A::Vector{Complex{$T}}, B::Vector{$T})
            result = similar(A)
            zrvadd!(result, A, B)
        end
    end

    # zrvsub: C = A - B (complex - real) — NOTE: vDSP_zrvsub computes B - A, swap
    @eval begin
        function zrvsub!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{$T})
            length(A) == length(B) == length(result) || throw(DimensionMismatch("zrvsub!: A, B, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B result begin
                asplit = _split_view($DSPSplit, pointer(A))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zrvsub", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, B, 1, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zrvsub(A::Vector{Complex{$T}}, B::Vector{$T})
            result = similar(A)
            zrvsub!(result, A, B)
        end
    end

    # --- Complex compound ops ---

    # zvcma: D = conj(A)*B + C
    @eval begin
        function zvcma!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{Complex{$T}}, C::Vector{Complex{$T}})
            length(A) == length(B) == length(C) == length(result) || throw(DimensionMismatch("zvcma!: A, B, C, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B C result begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                csplit = _split_view($DSPSplit, pointer(C))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvcma", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, bsplit, _CSTRIDE, csplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zvcma(A::Vector{Complex{$T}}, B::Vector{Complex{$T}}, C::Vector{Complex{$T}})
            result = similar(A)
            zvcma!(result, A, B, C)
        end
    end

    # zvma: D = A*B + C
    @eval begin
        function zvma!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, B::Vector{Complex{$T}}, C::Vector{Complex{$T}})
            length(A) == length(B) == length(C) == length(result) || throw(DimensionMismatch("zvma!: A, B, C, and result must have equal lengths"))
            n = length(A)
            GC.@preserve A B C result begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                csplit = _split_view($DSPSplit, pointer(C))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvma", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, bsplit, _CSTRIDE, csplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zvma(A::Vector{Complex{$T}}, B::Vector{Complex{$T}}, C::Vector{Complex{$T}})
            result = similar(A)
            zvma!(result, A, B, C)
        end
    end

    # zvsma: D = A*b[scalar] + C
    @eval begin
        function zvsma!(result::Vector{Complex{$T}}, A::Vector{Complex{$T}}, b::Complex{$T}, C::Vector{Complex{$T}})
            length(A) == length(C) == length(result) || throw(DimensionMismatch("zvsma!: A, C, and result must have equal lengths"))
            n = length(A)
            b2 = $T[real(b), imag(b)]
            GC.@preserve A C result b2 begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(b2))
                csplit = _split_view($DSPSplit, pointer(C))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvsma", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, _CSTRIDE, bsplit, csplit, _CSTRIDE, osplit, _CSTRIDE, n)
            end
            return result
        end
        function zvsma(A::Vector{Complex{$T}}, b::Complex{$T}, C::Vector{Complex{$T}})
            result = similar(A)
            zvsma!(result, A, b, C)
        end
    end

    # --- Complex dot products ---

    # zidotpr: conjugate dot product — sum(conj(A) .* B)
    @eval begin
        function zidotpr(A::Vector{Complex{$T}}, B::Vector{Complex{$T}})
            length(A) == length(B) || throw(DimensionMismatch("zidotpr: A and B must have equal lengths"))
            n = length(A)
            out = $T[0, 0]
            GC.@preserve A B out begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                osplit = _split_view($DSPSplit, pointer(out))
                ccall(($(string("vDSP_zidotpr", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, UInt64),
                      asplit, _CSTRIDE, bsplit, _CSTRIDE, osplit, n)
            end
            return complex(out[1], out[2])
        end
    end

    # zrdotpr: complex-real dot product — sum(A .* B) where B is real
    @eval begin
        function zrdotpr(A::Vector{Complex{$T}}, B::Vector{$T})
            length(A) == length(B) || throw(DimensionMismatch("zrdotpr: A and B must have equal lengths"))
            n = length(A)
            out = $T[0, 0]
            GC.@preserve A B out begin
                asplit = _split_view($DSPSplit, pointer(A))
                osplit = _split_view($DSPSplit, pointer(out))
                ccall(($(string("vDSP_zrdotpr", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, UInt64),
                      asplit, _CSTRIDE, B, 1, osplit, n)
            end
            return complex(out[1], out[2])
        end
    end

    # --- Complex fill ---

    # zvfill: fill complex vector with complex scalar
    @eval begin
        function zvfill!(result::Vector{Complex{$T}}, c::Complex{$T})
            n = length(result)
            c2 = $T[real(c), imag(c)]
            GC.@preserve result c2 begin
                csplit = _split_view($DSPSplit, pointer(c2))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zvfill", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Ref{$DSPSplit}, Int64, UInt64),
                      csplit, osplit, _CSTRIDE, n)
            end
            return result
        end
    end

    # --- Complex convolution ---

    # zconv: complex convolution
    @eval begin
        function zconv!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}}, K::Vector{Complex{$T}})
            xn = length(X)
            kn = length(K)
            rn = length(result)
            # Pad X with (kn-1) leading and kn trailing zeros, matching the real conv path.
            xpad = zeros(Complex{$T}, xn + 2kn - 1)
            copyto!(xpad, kn, X, 1, xn)
            GC.@preserve xpad K result begin
                xsplit = _split_view($DSPSplit, pointer(xpad))
                ksplit = _split_view($DSPSplit, pointer(K))
                osplit = _split_view($DSPSplit, pointer(result))
                ccall(($(string("vDSP_zconv", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64, UInt64),
                      xsplit, _CSTRIDE, ksplit, _CSTRIDE, osplit, _CSTRIDE, rn, kn)
            end
            return result
        end
        function zconv(X::Vector{Complex{$T}}, K::Vector{Complex{$T}})
            rn = length(X) + length(K) - 1
            result = Vector{Complex{$T}}(undef, rn)
            zconv!(result, X, K)
        end
    end

    # --- Complex matrix multiply ---

    # zmmul: complex matrix multiply C = A * B
    @eval begin
        function zmmul!(C::Matrix{Complex{$T}}, A::Matrix{Complex{$T}}, B::Matrix{Complex{$T}})
            m, p = size(A)
            p2, n = size(B)
            p == p2 || throw(DimensionMismatch("A columns ($p) ≠ B rows ($p2)"))
            size(C) == (m, n) || throw(DimensionMismatch("C must be $m×$n"))
            GC.@preserve A B C begin
                asplit = _split_view($DSPSplit, pointer(A))
                bsplit = _split_view($DSPSplit, pointer(B))
                osplit = _split_view($DSPSplit, pointer(C))
                # Same trick as mmul: swap for col-major → row-major
                ccall(($(string("vDSP_zmmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64, UInt64, UInt64),
                      bsplit, _CSTRIDE, asplit, _CSTRIDE, osplit, _CSTRIDE, UInt64(n), UInt64(m), UInt64(p))
            end
            return C
        end
        function zmmul(A::Matrix{Complex{$T}}, B::Matrix{Complex{$T}})
            m = size(A, 1)
            n = size(B, 2)
            C = Matrix{Complex{$T}}(undef, m, n)
            zmmul!(C, A, B)
        end
    end
end

@doc "Complex vector addition: `C = A + B`. Wraps [`vDSP_zvadd`](https://developer.apple.com/documentation/accelerate/vdsp_zvadd)." zvadd
@doc "Complex vector subtraction: `C = A - B`. Wraps [`vDSP_zvsub`](https://developer.apple.com/documentation/accelerate/vdsp_zvsub)." zvsub
@doc "Complex conjugate multiply: `C = conj(A) * B`. Wraps [`vDSP_zvcmul`](https://developer.apple.com/documentation/accelerate/vdsp_zvcmul)." zvcmul
@doc "Complex-real multiply: `C = A * B` (complex * real). Wraps [`vDSP_zrvmul`](https://developer.apple.com/documentation/accelerate/vdsp_zrvmul)." zrvmul
@doc "Complex-real divide: `C = A / B` (complex / real). Wraps [`vDSP_zrvdiv`](https://developer.apple.com/documentation/accelerate/vdsp_zrvdiv)." zrvdiv
@doc "Complex-real add: adds real vector to real part of complex. Wraps [`vDSP_zrvadd`](https://developer.apple.com/documentation/accelerate/vdsp_zrvadd)." zrvadd
@doc "Complex-real subtract: subtracts real from complex. Wraps [`vDSP_zrvsub`](https://developer.apple.com/documentation/accelerate/vdsp_zrvsub)." zrvsub
@doc "Complex conjugate multiply and add: `D = conj(A)*B + C`. Wraps [`vDSP_zvcma`](https://developer.apple.com/documentation/accelerate/vdsp_zvcma)." zvcma
@doc "Complex multiply and add: `D = A*B + C`. Wraps [`vDSP_zvma`](https://developer.apple.com/documentation/accelerate/vdsp_zvma)." zvma
@doc "Complex scalar multiply and add: `D = A*b + C`. Wraps [`vDSP_zvsma`](https://developer.apple.com/documentation/accelerate/vdsp_zvsma)." zvsma
@doc "Conjugate dot product: `sum(conj(A) .* B)`. Wraps [`vDSP_zidotpr`](https://developer.apple.com/documentation/accelerate/vdsp_zidotpr)." zidotpr
@doc "Complex-real dot product: `sum(A .* B)` where B is real. Wraps [`vDSP_zrdotpr`](https://developer.apple.com/documentation/accelerate/vdsp_zrdotpr)." zrdotpr
@doc "Fill complex vector with complex scalar. Wraps [`vDSP_zvfill`](https://developer.apple.com/documentation/accelerate/vdsp_zvfill)." zvfill!
@doc "Complex convolution. Wraps [`vDSP_zconv`](https://developer.apple.com/documentation/accelerate/vdsp_zconv)." zconv
@doc "Complex matrix multiply: `C = A * B`. Wraps [`vDSP_zmmul`](https://developer.apple.com/documentation/accelerate/vdsp_zmmul)." zmmul

# ============================================================
# Format conversion (interleaved ↔ split complex)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                             (Float64, "D", :DSPDoubleSplitComplex))
    @eval begin
        function ctoz(X::Vector{Complex{$T}})
            n = length(X)
            o_re = Vector{$T}(undef, n)
            o_im = Vector{$T}(undef, n)
            GC.@preserve o_re o_im begin
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_ctoz", suff)), libacc), Cvoid,
                      (Ptr{Complex{$T}}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      X, 2, osplit, 1, n)
            end
            return (o_re, o_im)
        end
        function ztoc(re::Vector{$T}, im::Vector{$T})
            n = length(re)
            result = Vector{Complex{$T}}(undef, n)
            GC.@preserve re im begin
                isplit = $DSPSplit(re, im)
                ccall(($(string("vDSP_ztoc", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{Complex{$T}}, Int64, UInt64),
                      isplit, 1, result, 2, n)
            end
            return result
        end
    end
end

@doc "Convert interleaved complex to split-complex (real, imag) vectors. Wraps [`vDSP_ctoz`](https://developer.apple.com/documentation/accelerate/vdsp_ctoz)." ctoz
@doc "Convert split-complex (real, imag) vectors to interleaved complex. Wraps [`vDSP_ztoc`](https://developer.apple.com/documentation/accelerate/vdsp_ztoc)." ztoc
