# AppleAccelerate.jl

This provides a Julia interface to some of the
[macOS Accelerate framework](https://developer.apple.com/documentation/accelerate). At
the moment, this package provides:
1. Access to Accelerate BLAS and LAPACK using the [libblastrampoline](https://github.com/JuliaLinearAlgebra/libblastrampoline) framework,
2. An interface to the [array-oriented functions](https://developer.apple.com/documentation/accelerate/veclib),
which provide a vectorised form for many common mathematical functions

The performance is significantly better than using standard libm functions in some cases, though there does appear to be some reduced accuracy.

## OS Requirements

MacOS 13.4 is required in order to run AppleAccelerate.jl, especially for the libblastrampoline forwarding. On older MacOS versions, your mileage may vary.

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

## Setting the number of threads

Accelerate is multithreaded by default. You can set an upper limit on the number of threads through the `VECLIB_MAXIMUM_THREADS` environment variable; for example, for single-threaded execution, start Julia as `VECLIB_MAXIMUM_THREADS=1 julia`. Accelerate does not support the `BLAS.set_num_threads(nthreads)` and `BLAS.get_num_threads()` API used by other BLAS backends (`set_num_threads` is a no-op and `get_num_threads` returns a hardcoded default).

## Example

To use the Accelerate BLAS and LAPACK, simply load the library:
```julia
julia> peakflops(4096)
3.6024175318268243e11

julia> using AppleAccelerate

julia> peakflops(4096)
5.832806459434183e11
```

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
