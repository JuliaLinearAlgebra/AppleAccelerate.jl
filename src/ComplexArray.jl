## ComplexArray.jl — Complex-valued vDSP operations using split-complex format ##

# ============================================================
# Split-complex structs (used here and by DSP.jl for FFT)
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
