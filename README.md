# OSXAccelerate.jl

[![Build Status](https://travis-ci.org/simonbyrne/OSXAccelerate.jl.svg?branch=master)](https://travis-ci.org/simonbyrne/OSXAccelerate.jl)

This provides a Julia interface to some of
[OS X's Accelerate framework](https://developer.apple.com/library/mac/documentation/Accelerate/Reference/AccelerateFWRef/). At
the moment, the package only provides an interface to the
[array-oriented functions](https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/index.html#//apple_ref/doc/uid/TP30000414-357225),
which provide a vectorised form of many common mathematical functions. In
general, the performance is significantly better than using standard libm
functions, though there does appear to be some reduced accuracy.

The following functions are supported:
 * *Rounding*: `ceil`, `floor`, `trunc`, `round`
 * *Logarithmic*: `exp`, `exp2`, `expm1`, `log`, `log1p`, `log2`, `log10`
 * *Trigonometric*: `sin`, `sinpi`, `cos`, `cospi`, `tan`, `tanpi`, `asin`, `acos`, `atan`, `atan2`
 * *Hyperbolic*: `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`
 * *Other*: `sqrt`, `copysign`, `exponent`, `abs`, `rem`

Note there are some slight differences from behaviour in Base:
 * No `DomainError`s are raised, instead `NaN` values are returned.
 * `round` breaks ties (values with a fractional part of 0.5) by choosing the
   nearest even value.
 * `exponent` returns a floating point value of the same type (instead of an `Int`).

Some additional functions that are also available:
* `rec(x)`: reciprocal (`1.0 ./ x`)
* `rsqrt(x)`: reciprocal square-root (`1.0 ./ sqrt(x)`)
* `pow(x,y)`: power (`x .^ y` in Base)
* `div(x,y)`: divide (`x ./ y` in Base)
* `sincos(x)`: returns `(sin(x), cos(x))`
* `cosisin(x)`: returns `cos(x)+im*sin(x)`

## Example

To avoid naming conflicts with Base, methods are not exported and so need to
be accessed via the namespace:
```julia
using OSXAccelerate
X = randn(1_000_000)
@time Y = exp(X) # standard libm function
@time Y = OSXAccelerate.exp(X) # OSXAccelerate array-oriented function
```

The `@replaceBase` macro replaces the relevant Base methods directly
```julia
OSXAccelerate.@replaceBase sin cos tan
OSXAccelerate.@replaceBase(.^, ./) # use parenthesised form for infix ops
@time sin(X) # will use OSXAccelerate methods for vectorised operations
```

Output arrays can be specified as first arguments of the functions suffixed
with `!`:
```julia
out = Array(Float64,1_000_000)
@time OSXAccelerate.exp!(out, X)
```

**Warning**: no dimension checks are performed on the `!` functions, so ensure
  your input and output arrays are of the same length.

Operations can be performed in-place by specifying the output array as the
input array (e.g. `OSXAccelerate.exp!(X,X)`), but this is not mentioned in the
OSXAccelerate docs, so do this at your own risk!
