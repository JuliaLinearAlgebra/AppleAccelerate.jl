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

"""
    integrate(f, a, b; integrator=:qags, abstol=1e-8, reltol=1e-8,
              max_intervals=200, qag_points=0)

Numerically integrate the real scalar function `f` over the interval `(a, b)`
using Apple's Accelerate Quadrature library (a C port of QUADPACK).

`f` is an ordinary Julia function `f(x::Float64) -> Float64`; Accelerate calls it
in batches, which is handled internally.

For `integrator=:qags`, one or both bounds may be infinite (`±Inf`).

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

    # Accelerate's callback fills y[i] = f(x[i]) for a batch of points. We close
    # over `f` directly, so the `arg` user-pointer is unused and we pass C_NULL.
    fill! = (_arg, n, x, y) -> begin
        xs  = unsafe_wrap(Array, x, (Int(n),))
        out = unsafe_wrap(Array, y, (Int(n),))
        @inbounds for i in 1:Int(n)
            out[i] = f(xs[i])
        end
        return nothing
    end
    cb = @cfunction($fill!, Cvoid, (Ptr{Cvoid}, Csize_t, Ptr{Cdouble}, Ptr{Cdouble}))

    fun = Ref(quadrature_integrate_function(
        Base.unsafe_convert(Ptr{Cvoid}, cb), C_NULL))
    opts = Ref(quadrature_integrate_options(
        intg, Float64(abstol), Float64(reltol),
        Csize_t(qag_points), Csize_t(max_intervals)))
    status = Ref(QUADRATURE_SUCCESS)
    abserr = Ref{Cdouble}(0.0)

    value = GC.@preserve cb fun opts status abserr begin
        quadrature_integrate(fun, Float64(a), Float64(b), opts,
                             status, abserr, Csize_t(0), C_NULL)
    end

    return (value = value, abserr = abserr[], status = status[])
end
