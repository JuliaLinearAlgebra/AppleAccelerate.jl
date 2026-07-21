# Idiomatic numerical integration on top of the generated `LibAccelerate`
# Quadrature bindings. This is the hand-written layer: it owns the ergonomics
# (a Julia function `f`, keyword options, a tidy return value) while the struct
# layout, enum values, and ccall come from `src/lib/LibAccelerate.jl`.
#
# This file is the proof-of-concept for the Clang.jl-based architecture: note
# that it never names a C struct field offset or an enum integer — everything
# ABI-level is imported from the generated module.

using .LibAccelerate:
    quadrature_integrate,
    quadrature_integrate_function,
    quadrature_integrate_options,
    quadrature_status,
    QUADRATURE_SUCCESS,
    QUADRATURE_INTEGRATE_QNG,
    QUADRATURE_INTEGRATE_QAG,
    QUADRATURE_INTEGRATE_QAGS

# Mutable box passed to the C trampoline via the user-pointer. Holds the Julia
# integrand and a slot for any exception it throws, so the error can be carried
# back across the C QUADPACK frames and rethrown in `integrate` (throwing
# directly through the C frames is undefined behavior).
mutable struct _QuadratureBox
    f::Any
    err::Any
end

# Top-level (non-closure) callback for Accelerate's batched integrand evaluation. We recover
# the Julia integrand from the C `fun_arg` user-pointer rather than closing over it, because
# closure `@cfunction` is unsupported on macOS/aarch64 before Julia 1.12; this trampoline
# works on all supported platforms. `arg` points to a `Ref{Any}` boxing a `_QuadratureBox`.
function _quadrature_trampoline(arg::Ptr{Cvoid}, n::Csize_t,
                                x::Ptr{Cdouble}, y::Ptr{Cdouble})
    box = unsafe_pointer_to_objref(arg)[]::_QuadratureBox
    xs  = unsafe_wrap(Array, x, (Int(n),))
    out = unsafe_wrap(Array, y, (Int(n),))
    # If a previous batch already errored, stop evaluating; just fill sentinels.
    if box.err !== nothing
        @inbounds for i in 1:Int(n)
            out[i] = NaN
        end
        return nothing
    end
    f = box.f
    @inbounds for i in 1:Int(n)
        try
            out[i] = f(xs[i])
        catch e
            # Capture the first error and write a sentinel; do NOT let the
            # exception unwind through the C caller. `integrate` rethrows it.
            box.err = e
            for j in i:Int(n)
                out[j] = NaN
            end
            return nothing
        end
    end
    return nothing
end

"""
    integrate(f, a, b; integrator=:qags, abstol=1e-8, reltol=1e-8,
              max_intervals=200, qag_points=0)

Numerically integrate the real scalar function `f` over the interval `(a, b)`
using Apple's Accelerate Quadrature library (a C port of QUADPACK).

`f` is an ordinary Julia function `f(x::Float64) -> Float64`; Accelerate calls it
in batches, which is handled internally.

For `integrator=:qags`, one or both bounds may be infinite (`±Inf`).

For `integrator=:qag`, `qag_points` selects the Gauss-Kronrod rule and must be one
of `15`, `21`, `31`, `41`, `51`, `61` (or `0` for the library default). Any other
value throws an `ArgumentError` (the underlying library would otherwise silently
return `value = 0.0`).

Returns a named tuple `(value, abserr, status)` where `value` is the integral
estimate, `abserr` is the estimated absolute error, and `status` is the
`quadrature_status` enum (`QUADRATURE_SUCCESS` on success).

# Examples
```julia
integrate(x -> x^2, 0, 1).value          # ≈ 0.3333333
integrate(sin, 0, π).value               # ≈ 2.0
integrate(x -> exp(-x^2), -Inf, Inf).value  # ≈ √π
```
"""
function integrate(f, a::Real, b::Real;
                   integrator::Symbol = :qags,
                   abstol::Real = 1e-8,
                   reltol::Real = 1e-8,
                   max_intervals::Integer = 200,
                   qag_points::Integer = 0)

    intg = integrator === :qng  ? QUADRATURE_INTEGRATE_QNG  :
           integrator === :qag  ? QUADRATURE_INTEGRATE_QAG  :
           integrator === :qags ? QUADRATURE_INTEGRATE_QAGS :
           throw(ArgumentError("unknown integrator $(repr(integrator)); use :qng, :qag, or :qags"))

    # For QAG, only the six Gauss-Kronrod rules below are valid; any other nonzero
    # value makes the library silently return value = 0.0, so reject it up front.
    qag_points in (0, 15, 21, 31, 41, 51, 61) ||
        throw(ArgumentError("qag_points must be one of 0 (default), 15, 21, 31, 41, " *
                            "51, 61; got $(qag_points)"))

    # Only QAGS handles infinite/semi-infinite intervals; QNG and QAG operate on
    # finite intervals and would otherwise silently return a meaningless result.
    if (isinf(a) || isinf(b)) && intg != QUADRATURE_INTEGRATE_QAGS
        throw(ArgumentError("infinite bounds (±Inf) require integrator=:qags; " *
                            "got integrator=$(repr(integrator))"))
    end

    # Pass the integrand to the top-level trampoline via the C user-pointer. `fbox` boxes a
    # mutable `_QuadratureBox` in a Ref so we can take its address (pointer_from_objref needs a
    # mutable object); GC.@preserve keeps it alive across the ccall. The box also carries back
    # any exception the integrand throws so we can rethrow it without unwinding through C.
    cb = @cfunction(_quadrature_trampoline, Cvoid,
                    (Ptr{Cvoid}, Csize_t, Ptr{Cdouble}, Ptr{Cdouble}))
    fbox = Ref{Any}(_QuadratureBox(f, nothing))
    fun = Ref(quadrature_integrate_function(cb, pointer_from_objref(fbox)))
    opts = Ref(quadrature_integrate_options(
        intg, Float64(abstol), Float64(reltol),
        Csize_t(qag_points), Csize_t(max_intervals)))
    status = Ref(QUADRATURE_SUCCESS)
    abserr = Ref{Cdouble}(0.0)

    value = GC.@preserve fbox fun opts status abserr begin
        quadrature_integrate(fun, Float64(a), Float64(b), opts,
                             status, abserr, Csize_t(0), C_NULL)
    end

    # If the integrand threw during evaluation, surface that error now (rather
    # than returning a meaningless value/status).
    captured = (fbox[]::_QuadratureBox).err
    captured === nothing || throw(captured)

    return (value = value, abserr = abserr[], status = status[])
end
