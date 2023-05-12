# AppleAccelerate.jl

This provides a Julia interface to some of the
[macOS Accelerate framework](https://developer.apple.com/documentation/accelerate). At
the moment, the package consists mostly of an interface to the
[array-oriented functions](https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/index.html#//apple_ref/doc/uid/TP30000414-357225),
which provide a vectorised form of many common mathematical functions, however the package does provide access to several other vectorized operations and more are being added on a regular basis.   In
general, the performance is significantly better than using standard libm
functions, though there does appear to be some reduced accuracy.

## OS Requirements

MacOS 13.3 is required in order to run AppleAccelerate.jl. On older MacOS versions, your mileage may vary.

## Supported Functions

The following functions are supported:
 * *Rounding*: `ceil`, `floor`, `trunc`, `round`
 * *Logarithmic*: `exp`, `exp2`, `expm1`, `log`, `log1p`, `log2`, `log10`
 * *Trigonometric*: `sin`, `sinpi`, `cos`, `cospi`, `tan`, `tanpi`, `asin`, `acos`, `atan`, `atan2`, `cis`
 * *Hyperbolic*: `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`
 * *Convolution*: `conv`, `xcorr`
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
* `fdiv(x,y)`: divide (`x ./ y` in Base)
* `sincos(x)`: returns `(sin(x), cos(x))`

## Example

To avoid naming conflicts with Base, methods are not exported and so need to
be accessed via the namespace:
```julia
using AppleAccelerate
using BenchmarkTools
X = randn(1_000_000);
@btime exp.($X); # standard libm function
@btime AppleAccelerate.exp($X); # Accelerate array-oriented function
```

The `@replaceBase` macro replaces the relevant Base methods directly
```julia
@btime sin.($X); # standard libm function
AppleAccelerate.@replaceBase sin cos tan
@btime sin($X);  # will use AppleAccelerate methods for vectorised operations

X = randn(1_000_000);
Y = fill(3.0, 1_000_000);
@btime $X .^ $Y;
AppleAccelerate.@replaceBase(^, /) # use parenthesised form for infix ops
@btime $X ^ $Y;
```

Output arrays can be specified as first arguments of the functions suffixed
with `!`:
```julia
out = zeros(Float64, 1_000_000)
@btime AppleAccelerate.exp!($out, $X)
```

**Warning**: no dimension checks are performed on the `!` functions, so ensure
  your input and output arrays are of the same length.

Operations can be performed in-place by specifying the output array as the
input array (e.g. `AppleAccelerate.exp!(X,X)`). This is not mentioned in the
Accelerate docs, but [this comment](http://stackoverflow.com/a/28833191/392585) by one of the authors indicates that it is safe.
