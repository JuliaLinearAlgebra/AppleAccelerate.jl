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

# ============================================================
# Complex → Complex unary operations (split-complex in/out)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))

    # vneg: negate complex vector
    @eval begin
        function vneg!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}})
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve x_re x_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvneg", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve x_re x_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvconj", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            x_re = real.(X); x_im = imag.(X)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve x_re x_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvmov", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, 1, osplit, 1, n)
            end
            return complex.(o_re, o_im)
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
# Complex → Complex binary operations (split-complex in/out)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))

    # vmul: element-wise complex multiply (conjugate flag = +1 for normal multiply)
    @eval begin
        function vmul!(result::Vector{Complex{$T}}, X::Vector{Complex{$T}}, Y::Vector{Complex{$T}})
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            y_re = real.(Y); y_im = imag.(Y)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve x_re x_im y_re y_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ysplit = $DSPSplit(y_re, y_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64, Cint),
                      xsplit, 1, ysplit, 1, osplit, 1, n, Cint(1))
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            y_re = real.(Y); y_im = imag.(Y)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve x_re x_im y_re y_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ysplit = $DSPSplit(y_re, y_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvdiv", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      ysplit, 1, xsplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            c_re = $T[real(c)]; c_im = $T[imag(c)]
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve x_re x_im c_re c_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                csplit = $DSPSplit(c_re, c_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvzsml", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Ref{$DSPSplit}, Int64, UInt64),
                      xsplit, 1, csplit, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
            return result
        end
        function vsmul(X::Vector{Complex{$T}}, c::Complex{$T})
            result = similar(X)
            vsmul!(result, X, c)
        end
    end
end

# ============================================================
# Complex → Real operations (split-complex in, real out)
# ============================================================
for (T, suff, DSPSplit) in ((Float32, "", :DSPSplitComplex),
                            (Float64, "D", :DSPDoubleSplitComplex))

    # vabs: complex absolute value (modulus)
    @eval begin
        function vabs!(result::Vector{$T}, X::Vector{Complex{$T}})
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            GC.@preserve x_re x_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ccall(($(string("vDSP_zvabs", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, 1, result, 1, n)
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            GC.@preserve x_re x_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ccall(($(string("vDSP_zvphas", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, 1, result, 1, n)
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            GC.@preserve x_re x_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ccall(($(string("vDSP_zvmags", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, 1, result, 1, n)
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            GC.@preserve x_re x_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ccall(($(string("vDSP_zvmgsa", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                      xsplit, 1, B, 1, result, 1, n)
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
            n = length(X)
            x_re = real.(X); x_im = imag.(X)
            y_re = real.(Y); y_im = imag.(Y)
            o_re = $T[0]; o_im = $T[0]
            GC.@preserve x_re x_im y_re y_im o_re o_im begin
                xsplit = $DSPSplit(x_re, x_im)
                ysplit = $DSPSplit(y_re, y_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zdotpr", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, UInt64),
                      xsplit, 1, ysplit, 1, osplit, n)
            end
            return complex(o_re[1], o_im[1])
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
            interleaved_in = reinterpret($T, X)
            interleaved_out = Vector{$T}(undef, 2 * n)
            ccall(($(string("vDSP_polar", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  interleaved_in, 2, interleaved_out, 2, n)
            magnitudes = interleaved_out[1:2:end]
            angles = interleaved_out[2:2:end]
            return (magnitudes, angles)
        end
    end

    # rect: (magnitude, angle) → Cartesian complex
    # vDSP_rect takes interleaved [mag, angle, mag, angle, ...] and writes [re, im, re, im, ...]
    @eval begin
        function rect(magnitudes::Vector{$T}, angles::Vector{$T})
            n = length(magnitudes)
            interleaved_in = Vector{$T}(undef, 2 * n)
            interleaved_in[1:2:end] .= magnitudes
            interleaved_in[2:2:end] .= angles
            interleaved_out = Vector{$T}(undef, 2 * n)
            ccall(($(string("vDSP_rect", suff)), libacc), Cvoid,
                  (Ptr{$T}, Int64, Ptr{$T}, Int64, UInt64),
                  interleaved_in, 2, interleaved_out, 2, n)
            return reinterpret(Complex{$T}, interleaved_out)
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im b_re b_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvadd", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, bsplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im b_re b_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvsub", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, bsplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im b_re b_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvcmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, bsplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zrvmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, B, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zrvdiv", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, B, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zrvadd", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, B, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zrvsub", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, B, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            c_re = real.(C); c_im = imag.(C)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im b_re b_im c_re c_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                csplit = $DSPSplit(c_re, c_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvcma", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, bsplit, 1, csplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            c_re = real.(C); c_im = imag.(C)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im b_re b_im c_re c_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                csplit = $DSPSplit(c_re, c_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvma", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, bsplit, 1, csplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = $T[real(b)]; b_im = $T[imag(b)]
            c_re = real.(C); c_im = imag.(C)
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve a_re a_im b_re b_im c_re c_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                csplit = $DSPSplit(c_re, c_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvsma", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64),
                      asplit, 1, bsplit, csplit, 1, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            o_re = $T[0]; o_im = $T[0]
            GC.@preserve a_re a_im b_re b_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                bsplit = $DSPSplit(b_re, b_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zidotpr", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, UInt64),
                      asplit, 1, bsplit, 1, osplit, n)
            end
            return complex(o_re[1], o_im[1])
        end
    end

    # zrdotpr: complex-real dot product — sum(A .* B) where B is real
    @eval begin
        function zrdotpr(A::Vector{Complex{$T}}, B::Vector{$T})
            n = length(A)
            a_re = real.(A); a_im = imag.(A)
            o_re = $T[0]; o_im = $T[0]
            GC.@preserve a_re a_im o_re o_im begin
                asplit = $DSPSplit(a_re, a_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zrdotpr", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ptr{$T}, Int64, Ref{$DSPSplit}, UInt64),
                      asplit, 1, B, 1, osplit, n)
            end
            return complex(o_re[1], o_im[1])
        end
    end

    # --- Complex fill ---

    # zvfill: fill complex vector with complex scalar
    @eval begin
        function zvfill!(result::Vector{Complex{$T}}, c::Complex{$T})
            n = length(result)
            c_re = $T[real(c)]; c_im = $T[imag(c)]
            o_re = Vector{$T}(undef, n); o_im = Vector{$T}(undef, n)
            GC.@preserve c_re c_im o_re o_im begin
                csplit = $DSPSplit(c_re, c_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zvfill", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Ref{$DSPSplit}, Int64, UInt64),
                      csplit, osplit, 1, n)
            end
            copyto!(result, complex.(o_re, o_im))
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
            x_re = real.(X); x_im = imag.(X)
            k_re = real.(K); k_im = imag.(K)
            o_re = Vector{$T}(undef, rn); o_im = Vector{$T}(undef, rn)
            # Pad X with kn-1 zeros like real conv
            xpad_re = [zeros($T, kn-1); x_re; zeros($T, kn)]
            xpad_im = [zeros($T, kn-1); x_im; zeros($T, kn)]
            GC.@preserve xpad_re xpad_im k_re k_im o_re o_im begin
                xsplit = $DSPSplit(xpad_re, xpad_im)
                ksplit = $DSPSplit(k_re, k_im)
                osplit = $DSPSplit(o_re, o_im)
                ccall(($(string("vDSP_zconv", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64, UInt64),
                      xsplit, 1, ksplit, 1, osplit, 1, rn, kn)
            end
            copyto!(result, complex.(o_re, o_im))
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
            a_re = real.(A); a_im = imag.(A)
            b_re = real.(B); b_im = imag.(B)
            o_re = Matrix{$T}(undef, m, n); o_im = Matrix{$T}(undef, m, n)
            GC.@preserve a_re a_im b_re b_im o_re o_im begin
                asplit = $DSPSplit(pointer(a_re), pointer(a_im))
                bsplit = $DSPSplit(pointer(b_re), pointer(b_im))
                osplit = $DSPSplit(pointer(o_re), pointer(o_im))
                # Same trick as mmul: swap for col-major → row-major
                ccall(($(string("vDSP_zmmul", suff)), libacc), Cvoid,
                      (Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, Ref{$DSPSplit}, Int64, UInt64, UInt64, UInt64),
                      bsplit, 1, asplit, 1, osplit, 1, UInt64(n), UInt64(m), UInt64(p))
            end
            copyto!(C, complex.(o_re, o_im))
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
