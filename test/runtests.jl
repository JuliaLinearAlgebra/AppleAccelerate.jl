using Accelerate
using FactCheck

srand(7)

N = 1_000

facts("Rounding") do
    X = 100*randn(N)
    for f in [:floor,:ceil,:trunc,:round]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => fb(X) "Mismatch for $f"
    end
end

              
facts("Logarithmic/exponential") do
    X = 100*randn(N)
    for f in [:exp,:exp2,:expm1]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
    X = exp(10*randn(N))
    for f in [:log,:log2,:log10]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
    X = expm1(10*randn(N))
    for f in [:log1p]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
end

facts("Trigonometric") do
    X = 10*randn(N)
    for f in [:sin,:sinpi,:cos,:cospi,:tan,:atan] # tanpi not defined in Base
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    Y = 10*randn(N)
    for f in [:atan2]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X,Y) => roughly(fb(X,Y)) "Mismatch for $f"
    end

    X = 2*rand(N)-1
    for f in [:asin,:acos] # tanpi not defined in Base
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
end


facts("Hyperbolic") do
    X = 10*randn(N)
    for f in [:sinh,:cosh,:tanh,:asinh] 
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    X = exp(10*randn(N))+1
    for f in [:acosh] 
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    X = 2*rand(N)-1
    for f in [:atanh] 
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
end


facts("Misc.") do
    X = exp(10*randn(N))
    for f in [:sqrt]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    X = 10*randn(N)
    for f in [:exponent, :abs]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    Y = 10*randn(N)
    for f in [:copysign]
        @eval fb = $f
        @eval fa = Accelerate.$f
        @fact fa(X,Y) => roughly(fb(X,Y)) "Mismatch for $f"
    end
    @fact Accelerate.rem(X,Y) => [rem(X[i], Y[i]) for i=1:length(X)] # no vectorized rem
end

facts("Extra") do
    X = randn(N)
    Y = exp(randn(N))

    @fact Accelerate.rec(X) => roughly(1./X)
    @fact Accelerate.rsqrt(Y) => roughly(1./sqrt(Y))
    @fact Accelerate.pow(Y,X) => roughly(Y.^X)
    @fact Accelerate.div(X,Y) => roughly(X./Y)

    @fact [Accelerate.sincos(X)...] => roughly([sin(X),cos(X)])
    @fact Accelerate.cosisin(X) => roughly(cos(X)+im*sin(X))
end

FactCheck.exitstatus()
