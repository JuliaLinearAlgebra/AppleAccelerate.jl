# Accelerate.jl

[![Build Status](https://travis-ci.org/simonbyrne/Accelerate.jl.svg?branch=master)](https://travis-ci.org/simonbyrne/Accelerate.jl)

This provides a Julia interface to some of OS X's
[Accelerate framework](https://developer.apple.com/library/mac/documentation/Accelerate/Reference/AccelerateFWRef/). At
the moment, the package only provides an interface to the
[array-oriented functions](https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/index.html#//apple_ref/doc/uid/TP30000414-357225),
which provide a vectorised form of many common mathemetical functions. In
general, the performance is significantly better than using standard libm
functions, though there does appear to be some reduced accuracy.

## Example

```julia
using Accelerate
X = randn(1_000_000)
@time Y = exp(X) # standard libm function
@time Y = Accelerate.exp(X) # Accelerate array-oriented function
```
