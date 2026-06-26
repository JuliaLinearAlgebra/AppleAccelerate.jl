# Numerical Integration

AppleAccelerate wraps Apple's
[Quadrature](https://developer.apple.com/documentation/accelerate/quadrature)
library (a C port of QUADPACK) for adaptive numerical integration of real scalar
functions.

```@setup quad
using AppleAccelerate
```

!!! note "Namespace"
    [`integrate`](@ref AppleAccelerate.integrate) is not exported. Access it via the
    `AppleAccelerate.` prefix.

## `integrate`

```julia
AppleAccelerate.integrate(f, a, b; integrator=:qags, abstol=1e-8, reltol=1e-8,
                          max_intervals=200, qag_points=0)
```

Integrate `f` over the interval `(a, b)`. `f` is an ordinary Julia function
`f(x::Float64) -> Float64`; Accelerate evaluates it in batches, which is handled
internally. The call returns a named tuple `(value, abserr, status)` where `value`
is the integral estimate, `abserr` is the estimated absolute error, and `status` is
the `quadrature_status` enum (`AppleAccelerate.LibAccelerate.QUADRATURE_SUCCESS` on
success).

```@example quad
r = AppleAccelerate.integrate(x -> x^2, 0, 1)
@assert isapprox(r.value, 1/3; atol = 1e-8)
r.value, r.abserr
```

```@example quad
AppleAccelerate.integrate(sin, 0, π).value   # ≈ 2
```

### Integrators

The `integrator` keyword selects the QUADPACK routine:

| `integrator` | Routine | Use case |
|--------------|---------|----------|
| `:qng`  | Non-adaptive Gauss-Kronrod | Smooth integrands; fastest, fixed rule |
| `:qag`  | Adaptive Gauss-Kronrod | General-purpose adaptive on a finite interval |
| `:qags` | Adaptive with extrapolation (default) | Endpoint singularities; supports infinite bounds |

For `:qag`, the `qag_points` keyword selects the Gauss-Kronrod rule and must be one
of `15`, `21`, `31`, `41`, `51`, `61` (or `0` for the library default). Any other
value throws an `ArgumentError` (the underlying library would otherwise silently
return `value = 0.0`).

```@example quad
AppleAccelerate.integrate(x -> exp(-x), 0, 5; integrator = :qag, qag_points = 21).value
```

### Infinite and reversed bounds

For `integrator = :qags`, one or both bounds may be infinite (`±Inf`):

```@example quad
r = AppleAccelerate.integrate(x -> exp(-x^2), -Inf, Inf)
@assert isapprox(r.value, sqrt(π); atol = 1e-6)
r.value   # ≈ √π
```

Reversed bounds (`a > b`) follow the usual sign convention,
`∫ₐᵇ f = -∫ᵦᵃ f`:

```@example quad
fwd = AppleAccelerate.integrate(sin, 0, π).value
rev = AppleAccelerate.integrate(sin, π, 0).value
@assert isapprox(fwd, -rev; atol = 1e-8)
nothing # hide
```

### Tolerances and subdivision

`abstol` and `reltol` set the requested absolute and relative error; `max_intervals`
bounds the number of subintervals the adaptive routines may create. If the integrand
throws, that exception is captured and rethrown from `integrate` (rather than
returning a meaningless value).

```@docs
AppleAccelerate.integrate
```
