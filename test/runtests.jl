using OSXAccelerate
using FactCheck

srand(7)

N = 1_000

facts("Rounding") do
    X = 100*randn(N)
    for f in [:floor,:ceil,:trunc,:round]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => fb(X) "Mismatch for $f"
    end
end

              
facts("Logarithmic/exponential") do
    X = 100*randn(N)
    for f in [:exp,:exp2,:expm1]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
    X = exp(10*randn(N))
    for f in [:log,:log2,:log10]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
    X = expm1(10*randn(N))
    for f in [:log1p]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
end

facts("Trigonometric") do
    X = 10*randn(N)
    for f in [:sin,:sinpi,:cos,:cospi,:tan,:atan] # tanpi not defined in Base
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    Y = 10*randn(N)
    for f in [:atan2]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X,Y) => roughly(fb(X,Y)) "Mismatch for $f"
    end

    X = 2*rand(N)-1
    for f in [:asin,:acos] # tanpi not defined in Base
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
end


facts("Hyperbolic") do
    X = 10*randn(N)
    for f in [:sinh,:cosh,:tanh,:asinh] 
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    X = exp(10*randn(N))+1
    for f in [:acosh] 
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    X = 2*rand(N)-1
    for f in [:atanh] 
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end
end


facts("Misc.") do
    X = exp(10*randn(N))
    for f in [:sqrt]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    X = 10*randn(N)
    for f in [:exponent, :abs]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X) => roughly(fb(X)) "Mismatch for $f"
    end

    Y = 10*randn(N)
    for f in [:copysign]
        @eval fb = $f
        @eval fa = OSXAccelerate.$f
        @fact fa(X,Y) => roughly(fb(X,Y)) "Mismatch for $f"
    end
    @fact OSXAccelerate.rem(X,Y) => [rem(X[i], Y[i]) for i=1:length(X)] # no vectorized rem
end

facts("Extra") do
    X = randn(N)
    Y = abs(randn(N))

    @fact OSXAccelerate.rec(X) => roughly(1./X)
    @fact OSXAccelerate.rsqrt(Y) => roughly(1./sqrt(Y))
    @fact OSXAccelerate.pow(Y,X) => roughly(Y.^X)
    @fact OSXAccelerate.div(X,Y) => roughly(X./Y)

    @fact [OSXAccelerate.sincos(X)...] => roughly([sin(X),cos(X)])
    @fact OSXAccelerate.cosisin(X) => roughly(cos(X)+im*sin(X))
end


facts("Replace Base") do
    X = randn(N)
    Y = abs(randn(N))

    OSXAccelerate.@replaceBase(sin,atan2,./,.^)
    @fact sin(X) => OSXAccelerate.sin(X)
    @fact atan2(X,Y) => OSXAccelerate.atan2(X,Y)
    @fact X ./ Y => OSXAccelerate.div(X,Y)
    @fact Y ./ X => OSXAccelerate.div(Y,X)
end

FactCheck.exitstatus()
