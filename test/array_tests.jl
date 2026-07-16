@testset "Array Operations" begin

for T in (Float32, Float64)
    @testset "Element-wise Operators::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)
        # Vector-vector
        @test (X .+ Y) ≈ AppleAccelerate.vadd(X, Y)
        @test (X .- Y) ≈ AppleAccelerate.vsub(X, Y)
        @test (X .* Y) ≈ AppleAccelerate.vmul(X, Y)
        @test (X ./ Y) ≈ AppleAccelerate.vdiv(X, Y)

        # Vector-vector non-allocating
        AppleAccelerate.vadd!(Z, X, Y)
        @test (X .+ Y) ≈ Z
        AppleAccelerate.vsub!(Z, X, Y)
        @test (X .- Y) ≈ Z
        AppleAccelerate.vmul!(Z, X, Y)
        @test (X .* Y) ≈ Z
        AppleAccelerate.vdiv!(Z, X, Y)
        @test (X ./ Y) ≈ Z

        # Vector-vector broadcasting
        @test (X .+ Y) ≈ AppleAccelerate.vadd.(X, Y)
        @test (X .- Y) ≈ AppleAccelerate.vsub.(X, Y)
        @test (X .* Y) ≈ AppleAccelerate.vmul.(X, Y)
        @test (X ./ Y) ≈ AppleAccelerate.vdiv.(X, Y)

        #Vector-scalar
        c::T         = randn()
        @test (X .+ c) ≈ AppleAccelerate.vsadd.(X, c)
        @test (X .- c) ≈ AppleAccelerate.vssub.(X, c)
        @test (c .- X) ≈ AppleAccelerate.svsub.(X, c)
        @test (X .* c) ≈ AppleAccelerate.vsmul.(X, c)
        @test (X ./ c) ≈ AppleAccelerate.vsdiv.(X, c)

        #Vector-scalar non-allocating
        AppleAccelerate.vsadd!(Y, X, c)
        @test (X .+ c) ≈ Y
        AppleAccelerate.vssub!(Y, X, c)
        @test (X .- c) ≈ Y
        AppleAccelerate.svsub!(Y, X, c)
        @test (c .- X) ≈ Y
        AppleAccelerate.vsmul!(Y, X, c)
        @test (X .* c) ≈ Y
        AppleAccelerate.vsdiv!(Y, X, c)
        @test (X ./ c) ≈ Y

        #Vector-scalar broadcasting
        @test (X .+ c) ≈ AppleAccelerate.vsadd.(X, c)
        @test (X .- c) ≈ AppleAccelerate.vssub.(X, c)
        @test (c .- X) ≈ AppleAccelerate.svsub.(X, c)
        @test (X .* c) ≈ AppleAccelerate.vsmul.(X, c)
        @test (X ./ c) ≈ AppleAccelerate.vsdiv.(X, c)

        @test (X .+ Y .+ Y) ≈ AppleAccelerate.vadd.(X, Y .+ Y)
    end
end

for T in (Float32, Float64)
    @testset "Rounding::$T" begin
        X::Array{T} = 100*randn(N)
        @testset "Testing $f::$T" for f in [:floor,:ceil,:trunc,:round]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb.(X)
            @test fa.(X) ≈ fb.(X)
        end
    end
end


for T in (Float32, Float64)
    @testset "Logarithmic::$T" begin
        X::Array{T} = exp.(10*randn(N))
        @testset "Testing $f::$T" for f in [:log,:log2,:log10, :log1p]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb.(X)
            @test fa.(X) ≈ fb.(X)
        end
    end
end


for T in (Float32, Float64)
    @testset "Exponential::$T" begin
        @testset "Testing $f::$T" for f in [:exp,:exp2,:expm1]
            X::Array{T} = 10*randn(N)
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb.(X)
            @test fa.(X) ≈ fb.(X)
        end
    end
end


for T in (Float32, Float64)
    X::Array{T} = 10*randn(N)
    @testset "Trigonometric::$T" begin
        @testset "Testing $f::$T" for f in [:sin,:sinpi,:cos,:cospi,:tan,:tanpi,:atan]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb.(X)
            @test fa.(X) ≈ fb.(X)
        end

        Y::Array{T} = 10*randn(N)
        @testset "Testing $f::$T" for f in [:atan]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X,Y) ≈ fb.(X,Y)
            @test fa.(X,Y) ≈ fb.(X,Y)
        end

        Z::Array{T} = 2*rand(N) .- 1
        @testset "Testing $f::$T" for f in [:asin,:acos]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Z) ≈ fb.(Z)
            @test fa.(Z) ≈ fb.(Z)
        end
    end
end


for T in (Float32, Float64)
    @testset "Hyperbolic::$T" begin
        X = 10*randn(N)
        @testset "Testing $f::$T" for f in [:sinh,:cosh,:tanh,:asinh]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb.(X)
            @test fa.(X) ≈ fb.(X)
        end

        Y = exp.(10*randn(N)) .+ 1
        @testset "Testing $f::$T" for f in [:acosh]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Y) ≈ fb.(Y)
            @test fa.(Y) ≈ fb.(Y)
        end

        Z = 2*rand(N) .- 1
        @testset "Testing $f::$T" for f in [:atanh]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Z) ≈ fb.(Z)
            @test fa.(Z) ≈ fb.(Z)
        end
    end
end

for T in (Float32, Float64)
    @testset "Array Properties::$T" begin
        X::Array{T} = randn(N)
        @testset "Testing $f::$T" for f in [:maximum, :minimum, :mean, :sum]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X)
        end

        @testset "Testing $f::$T" for f in [:findmax, :findmin]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X)[1] ≈ fb(X)[1]
            @test fa(X)[2] == fb(X)[2]
            # Public index return type must remain Int
            @test typeof(fa(X)[2]) == Int
        end

        @testset "Testing meanmag::$T" begin
            @test AppleAccelerate.meanmag(X) ≈ mean(abs, X)
        end

        @testset "Testing meansqr::$T" begin
            @test AppleAccelerate.meansqr(X) ≈ mean(X .* X)
        end

        @testset "Testing meanssqr::$T" begin
            @test AppleAccelerate.meanssqr(X) ≈ mean(X .* abs.(X))
        end

        @testset "Testing summag::$T" begin
            @test AppleAccelerate.summag(X) ≈ sum(abs, X)
        end

        @testset "Testing sumsqr::$T" begin
            @test AppleAccelerate.sumsqr(X) ≈ sum(abs2, X)
        end

        @testset "Testing sumssqr::$T" begin
            @test AppleAccelerate.sumssqr(X) ≈ sum(X .* abs.(X))
        end

    end
end

for T in (Float32, Float64)
    @testset "Misc::$T" begin
        X::Array{T} = exp.(10*randn(N))
        @testset "Testing $f::$T" for f in [:sqrt]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb.(X)
            @test fa.(X) ≈ fb.(X)
        end

        Y::Array{T} = 10*randn(N)
        @testset "Testing $f::$T" for f in [:exponent, :abs]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(Y) ≈ fb.(Y)
            @test fa.(Y) ≈ fb.(Y)
        end

        Z::Array{T} = 10*randn(N)
        @testset "Testing $f::$T" for f in [:copysign]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X,Y) ≈ fb.(X,Y)
            @test fa.(X,Y) ≈ fb.(X,Y)
        end

        @test AppleAccelerate.rem(X,Y) == rem.(X, Y)

    end
end


for T in (Float32, Float64)
    @testset "Extra::$T" begin
        X::Array{T} = randn(N)
        Y::Array{T} = abs.(randn(N))

        @test AppleAccelerate.rec(X) ≈ 1 ./ X
        @test AppleAccelerate.rsqrt(Y) ≈ 1 ./ sqrt.(Y)
        @test AppleAccelerate.pow(Y,X) ≈ Y.^X
        @test AppleAccelerate.div_float(X,Y) ≈ X./Y

        @test AppleAccelerate.sincos(X)[1] ≈ sin.(X)
        @test AppleAccelerate.sincos(X)[2] ≈ cos.(X)
        @test AppleAccelerate.cis(X) ≈ cis.(X)

    end
end

for T in (Float32, Float64)
    @testset "libvMisc/vForce::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)

        # cbrt: cube root (defined for negative inputs too)
        @test AppleAccelerate.cbrt(X) ≈ cbrt.(X)
        AppleAccelerate.cbrt!(Z, X)
        @test Z ≈ cbrt.(X)

        # remainder: IEEE-754 remainder == rem(x, y, RoundNearest)
        @test AppleAccelerate.remainder(X, Y) ≈ rem.(X, Y, RoundNearest)
        AppleAccelerate.remainder!(Z, X, Y)
        @test Z ≈ rem.(X, Y, RoundNearest)

        # nextafter: next representable value from X toward Y
        nextref = [X[i] < Y[i] ? nextfloat(X[i]) :
                   X[i] > Y[i] ? prevfloat(X[i]) : Y[i] for i in eachindex(X)]
        @test AppleAccelerate.nextafter(X, Y) == nextref
        AppleAccelerate.nextafter!(Z, X, Y)
        @test Z == nextref

        # pows: vector base raised to a scalar exponent
        Xpos::Vector{T} = abs.(randn(N))
        y::T = randn()
        @test AppleAccelerate.pows(Xpos, y) ≈ Xpos .^ y
        Zp::Vector{T} = similar(Xpos)
        AppleAccelerate.pows!(Zp, Xpos, y)
        @test Zp ≈ Xpos .^ y
    end
end

for T in (Float32, Float64)
    @testset "vDSP Unary::$T" begin
        X::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)

        @test AppleAccelerate.vneg(X) ≈ -X
        AppleAccelerate.vneg!(Z, X)
        @test Z ≈ -X

        @test AppleAccelerate.vnabs(X) ≈ -abs.(X)
        AppleAccelerate.vnabs!(Z, X)
        @test Z ≈ -abs.(X)

        @test AppleAccelerate.vsq(X) ≈ X .^ 2
        AppleAccelerate.vsq!(Z, X)
        @test Z ≈ X .^ 2

        @test AppleAccelerate.vssq(X) ≈ X .* abs.(X)
        AppleAccelerate.vssq!(Z, X)
        @test Z ≈ X .* abs.(X)

        @test AppleAccelerate.vfrac(X) ≈ X .- trunc.(X)
        AppleAccelerate.vfrac!(Z, X)
        @test Z ≈ X .- trunc.(X)

        @test AppleAccelerate.vabs(X) ≈ abs.(X)
        AppleAccelerate.vabs!(Z, X)
        @test Z ≈ abs.(X)

        Y = copy(X)
        AppleAccelerate.vreverse!(Y)
        @test Y ≈ reverse(X)
        @test AppleAccelerate.vreverse(X) ≈ reverse(X)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Two-Vector::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)

        @test AppleAccelerate.vmax(X, Y) ≈ max.(X, Y)
        AppleAccelerate.vmax!(Z, X, Y)
        @test Z ≈ max.(X, Y)

        @test AppleAccelerate.vmin(X, Y) ≈ min.(X, Y)
        AppleAccelerate.vmin!(Z, X, Y)
        @test Z ≈ min.(X, Y)

        @test AppleAccelerate.vmaxmg(X, Y) ≈ max.(abs.(X), abs.(Y))
        AppleAccelerate.vmaxmg!(Z, X, Y)
        @test Z ≈ max.(abs.(X), abs.(Y))

        @test AppleAccelerate.vminmg(X, Y) ≈ min.(abs.(X), abs.(Y))
        AppleAccelerate.vminmg!(Z, X, Y)
        @test Z ≈ min.(abs.(X), abs.(Y))

        @test AppleAccelerate.vdist(X, Y) ≈ hypot.(X, Y)
        AppleAccelerate.vdist!(Z, X, Y)
        @test Z ≈ hypot.(X, Y)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Scalar-Vector Divide::$T" begin
        X::Vector{T} = randn(N) .+ T(2)  # avoid division by zero
        c::T = T(3.5)
        Z::Vector{T} = similar(X)

        @test AppleAccelerate.svdiv(X, c) ≈ c ./ X
        AppleAccelerate.svdiv!(Z, X, c)
        @test Z ≈ c ./ X
    end
end

for T in (Float32, Float64)
    @testset "vDSP Compound Arithmetic::$T" begin
        A::Vector{T} = randn(N)
        B::Vector{T} = randn(N)
        C::Vector{T} = randn(N)
        D::Vector{T} = randn(N)
        R::Vector{T} = similar(A)

        # vam: (A+B)*C
        @test AppleAccelerate.vam(A, B, C) ≈ (A .+ B) .* C
        AppleAccelerate.vam!(R, A, B, C)
        @test R ≈ (A .+ B) .* C

        # vsbm: (A-B)*C
        @test AppleAccelerate.vsbm(A, B, C) ≈ (A .- B) .* C
        AppleAccelerate.vsbm!(R, A, B, C)
        @test R ≈ (A .- B) .* C

        # vaam: (A+B)*(C+D)
        @test AppleAccelerate.vaam(A, B, C, D) ≈ (A .+ B) .* (C .+ D)
        AppleAccelerate.vaam!(R, A, B, C, D)
        @test R ≈ (A .+ B) .* (C .+ D)

        # vsbsbm: (A-B)*(C-D)
        @test AppleAccelerate.vsbsbm(A, B, C, D) ≈ (A .- B) .* (C .- D)
        AppleAccelerate.vsbsbm!(R, A, B, C, D)
        @test R ≈ (A .- B) .* (C .- D)

        # vasbm: (A+B)*(C-D)
        @test AppleAccelerate.vasbm(A, B, C, D) ≈ (A .+ B) .* (C .- D)
        AppleAccelerate.vasbm!(R, A, B, C, D)
        @test R ≈ (A .+ B) .* (C .- D)

        # vpythg: sqrt((A-C)^2 + (B-D)^2)
        @test AppleAccelerate.vpythg(A, B, C, D) ≈ sqrt.((A .- C) .^ 2 .+ (B .- D) .^ 2)

        # vasm: (A+B)*c
        c::T = T(2.5)
        @test AppleAccelerate.vasm(A, B, c) ≈ (A .+ B) .* c

        # vsbsm: (A-B)*c
        @test AppleAccelerate.vsbsm(A, B, c) ≈ (A .- B) .* c

        # vsma: A*b + C
        b::T = T(1.5)
        @test AppleAccelerate.vsma(A, b, C) ≈ A .* b .+ C

        # vsmsa: A*b + c
        @test AppleAccelerate.vsmsa(A, b, c) ≈ A .* b .+ c

        # vaddsub: returns (A+B, B-A)
        (add_r, sub_r) = AppleAccelerate.vaddsub(A, B)
        @test add_r ≈ A .+ B
        @test sub_r ≈ B .- A
    end
end

for T in (Float32, Float64)
    @testset "vDSP Reductions::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)

        @test AppleAccelerate.dot(X, Y) ≈ LinearAlgebra.dot(X, Y)
        @test AppleAccelerate.distancesq(X, Y) ≈ sum((X .- Y) .^ 2)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Clipping::$T" begin
        X::Vector{T} = randn(N)

        low::T = T(-0.5)
        high::T = T(0.5)
        @test AppleAccelerate.vclip(X, low, high) ≈ clamp.(X, low, high)
        Z = similar(X)
        AppleAccelerate.vclip!(Z, X, low, high)
        @test Z ≈ clamp.(X, low, high)

        # vthr
        thr::T = T(0.0)
        @test AppleAccelerate.vthr(X, thr) ≈ max.(X, thr)

        # vthres
        expected_vthres = [x >= thr ? x : T(0) for x in X]
        @test AppleAccelerate.vthres(X, thr) ≈ expected_vthres
    end
end

@testset "vDSP Type Conversion" begin
    X32 = randn(Float32, N)
    X64 = randn(Float64, N)

    @test AppleAccelerate.vdouble(X32) ≈ Float64.(X32)
    @test AppleAccelerate.vsingle(X64) ≈ Float32.(X64)
end

for T in (Float32, Float64)
    @testset "vDSP Ramp::$T" begin
        start::T = T(1.0)
        step::T = T(0.5)
        n = 100
        ramp = AppleAccelerate.vramp(start, step, n)
        expected = [start + i * step for i in 0:(n-1)]
        @test ramp ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP Integration::$T" begin
        n = 101
        X::Vector{T} = ones(T, n)
        step::T = T(0.1)

        # vtrapz: trapezoidal integration of constant 1 over [0, 10] should give cumulative
        trapz = AppleAccelerate.vtrapz(X, step)
        # Result has same length as X; trapz[1] = 0, trapz[end] ≈ (n-1)*step
        @test trapz[1] ≈ T(0) atol=eps(T)
        @test trapz[end] ≈ T((n - 1) * step) atol=T(0.01)

        # vswsum: sliding window sum
        X2::Vector{T} = ones(T, 20)
        sw = AppleAccelerate.vswsum(X2, 5)
        @test length(sw) == 16
        @test all(x -> isapprox(x, T(5), atol=eps(T)), sw)

        # vswmax: sliding window max
        X3::Vector{T} = T.(collect(1:20))
        swm = AppleAccelerate.vswmax(X3, 5)
        @test length(swm) == 16
        @test swm ≈ T.(collect(5:20))
    end
end

for T in (Float32, Float64)
    @testset "vDSP Interpolation::$T" begin
        A::Vector{T} = randn(N)
        B::Vector{T} = randn(N)
        t::T = T(0.3)

        # vintb: A + t*(B-A)
        @test AppleAccelerate.vintb(A, B, t) ≈ A .+ t .* (B .- A)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Polynomial::$T" begin
        # Evaluate p(x) = 2x^2 + 3x + 1
        # vDSP_vpoly uses coefficients in reverse order: [highest degree first]
        # coeffs = [a_P, a_{P-1}, ..., a_1, a_0] for P-degree poly
        coeffs::Vector{T} = T[2, 3, 1]
        X::Vector{T} = T.(collect(-5.0:0.5:5.0))
        result = AppleAccelerate.vpoly(coeffs, X)
        expected = 2 .* X .^ 2 .+ 3 .* X .+ 1
        @test result ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP Normalization::$T" begin
        X::Vector{T} = randn(N)
        (result, m, s) = AppleAccelerate.vnormalize(X)
        @test m ≈ Statistics.mean(X) atol=T(1e-4)
        if s > 0
            @test result ≈ (X .- m) ./ s atol=T(1e-3)
        end
    end
end

for T in (Float32, Float64)
    @testset "vDSP Decibel Conversion::$T" begin
        X::Vector{T} = abs.(randn(N)) .+ T(0.01)  # positive values
        ref::T = T(1.0)

        # Power: 10*log10(X/ref)
        result_pow = AppleAccelerate.vdbcon(X, ref, true)
        expected_pow = 10 .* log10.(X ./ ref)
        @test result_pow ≈ expected_pow atol=T(1e-3)

        # Amplitude: 20*log10(X/ref)
        result_amp = AppleAccelerate.vdbcon(X, ref, false)
        expected_amp = 20 .* log10.(X ./ ref)
        @test result_amp ≈ expected_amp atol=T(1e-3)
    end
end

@testset "Broadcasted-Broadcasted dispatch" begin
    # Test that binary vDSP ops correctly handle chained broadcasting
    # where both arguments are Broadcasted objects (not yet materialized arrays).
    # This exercises the disambiguating methods for vadd/vsub/vmul/vdiv.
    for T in (Float32, Float64)
        X = randn(T, N)
        Y = randn(T, N)
        A = randn(T, N)
        B = randn(T, N)

        # vadd.(abs.(X), abs.(Y)) — both args are Broadcasted
        @test AppleAccelerate.vadd.(AppleAccelerate.abs.(X), AppleAccelerate.abs.(Y)) ≈ abs.(X) .+ abs.(Y)
        # vsub
        @test AppleAccelerate.vsub.(AppleAccelerate.abs.(X), AppleAccelerate.abs.(Y)) ≈ abs.(X) .- abs.(Y)
        # vmul
        @test AppleAccelerate.vmul.(AppleAccelerate.abs.(X), AppleAccelerate.abs.(Y)) ≈ abs.(X) .* abs.(Y)
        # vdiv
        @test AppleAccelerate.vdiv.(AppleAccelerate.abs.(X) .+ T(1), AppleAccelerate.abs.(Y) .+ T(1)) ≈ (abs.(X) .+ 1) ./ (abs.(Y) .+ 1)
    end
end

for T in (Float32, Float64)
    @testset "Complex Unary::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Z::Vector{Complex{T}} = similar(X)

        # vneg
        @test AppleAccelerate.vneg(X) ≈ -X
        AppleAccelerate.vneg!(Z, X)
        @test Z ≈ -X

        # vconj
        @test AppleAccelerate.vconj(X) ≈ conj.(X)
        AppleAccelerate.vconj!(Z, X)
        @test Z ≈ conj.(X)

        # vcopy
        @test AppleAccelerate.vcopy(X) ≈ X
    end
end

for T in (Float32, Float64)
    @testset "Complex Binary::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Y::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Z::Vector{Complex{T}} = similar(X)

        # vmul
        @test AppleAccelerate.vmul(X, Y) ≈ X .* Y
        AppleAccelerate.vmul!(Z, X, Y)
        @test Z ≈ X .* Y

        # vdiv
        @test AppleAccelerate.vdiv(X, Y) ≈ X ./ Y
        AppleAccelerate.vdiv!(Z, X, Y)
        @test Z ≈ X ./ Y

        # vsmul (complex scalar)
        c::Complex{T} = complex(T(2.5), T(-1.3))
        @test AppleAccelerate.vsmul(X, c) ≈ X .* c
        AppleAccelerate.vsmul!(Z, X, c)
        @test Z ≈ X .* c
    end
end

for T in (Float32, Float64)
    @testset "Complex → Real::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        R::Vector{T} = similar(X, T)

        # vabs
        @test AppleAccelerate.vabs(X) ≈ abs.(X)
        AppleAccelerate.vabs!(R, X)
        @test R ≈ abs.(X)

        # vphase
        @test AppleAccelerate.vphase(X) ≈ angle.(X)
        AppleAccelerate.vphase!(R, X)
        @test R ≈ angle.(X)

        # vmags
        @test AppleAccelerate.vmags(X) ≈ abs2.(X)
        AppleAccelerate.vmags!(R, X)
        @test R ≈ abs2.(X)

        # vmagsa
        B::Vector{T} = randn(T, N)
        @test AppleAccelerate.vmagsa(X, B) ≈ abs2.(X) .+ B
        AppleAccelerate.vmagsa!(R, X, B)
        @test R ≈ abs2.(X) .+ B
    end
end

for T in (Float32, Float64)
    @testset "Complex Dot Product::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Y::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # dot: conjugated dot product matching LinearAlgebra.dot — sum(conj(X) .* Y)
        @test AppleAccelerate.dot(X, Y) ≈ LinearAlgebra.dot(X, Y) ≈ sum(conj.(X) .* Y)
        # dotu: unconjugated (bilinear) dot product — sum(X .* Y)
        @test AppleAccelerate.dotu(X, Y) ≈ sum(X .* Y)
    end
end

for T in (Float32, Float64)
    @testset "Polar/Rect Conversion::$T" begin
        X::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # polar → rect roundtrip
        (mags, angs) = AppleAccelerate.polar(X)
        @test mags ≈ abs.(X)
        @test angs ≈ angle.(X)

        Y = AppleAccelerate.rect(mags, angs)
        @test Y ≈ X atol=T(1e-4)

        # rect from known values
        m::Vector{T} = ones(T, 10)
        a::Vector{T} = zeros(T, 10)
        result = AppleAccelerate.rect(m, a)
        @test real.(result) ≈ ones(T, 10) atol=T(1e-6)
        @test imag.(result) ≈ zeros(T, 10) atol=T(1e-6)
    end
end

for T in (Float32, Float64)
    @testset "vDSP Stereo Ramp Multiply (vrampmul2)::$T" begin
        n = 100
        I0 = ones(T, n)
        I1 = ones(T, n) .* T(2)
        start = T(1.0)
        step = T(0.5)

        O0, O1 = AppleAccelerate.vrampmul2(I0, I1, start, step)

        # O0[i] = (start + i*step) * I0[i], O1[i] = (start + i*step) * I1[i]
        ramp = [start + i * step for i in 0:(n-1)]
        @test O0 ≈ ramp .* I0
        @test O1 ≈ ramp .* I1

        # In-place variant
        R0 = similar(I0)
        R1 = similar(I1)
        AppleAccelerate.vrampmul2!(R0, R1, I0, I1, start, step)
        @test R0 ≈ O0
        @test R1 ≈ O1
    end
end

for T in (Float32, Float64)
    @testset "vDSP Vector Linear Average (vavlin)::$T" begin
        n = 50
        A = randn(T, n)
        C = randn(T, n)
        weight = T(3.0)

        # Operand order matches the mutating variant: vavlin(C, A, weight).
        result = AppleAccelerate.vavlin(C, A, weight)
        expected = (C .* weight .+ A) ./ (weight + 1)
        @test result ≈ expected

        # In-place variant
        C2 = copy(C)
        AppleAccelerate.vavlin!(C2, A, weight)
        @test C2 ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP vtmerg::$T" begin
        n = 100
        X = randn(T, n)
        Y = randn(T, n)

        # vtmerg: C[n] = A[n] + (B[n] - A[n]) * n/(N-1)  (0-based)
        result = AppleAccelerate.vtmerg(X, Y)
        expected = T[X[i] + (Y[i] - X[i]) * (i - 1) / (n - 1) for i in 1:n]
        @test result ≈ expected

        Z = similar(X)
        AppleAccelerate.vtmerg!(Z, X, Y)
        @test Z ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP venvlp::$T" begin
        n = 50
        A = randn(T, n)   # upper bound per element
        B = randn(T, n)   # lower bound per element
        C = randn(T, n)   # signal

        result = AppleAccelerate.venvlp(A, B, C)
        # venvlp: if C[n] < B[n] || A[n] < C[n] then D[n] = C[n] else D[n] = 0
        ref = T[(C[i] < B[i] || A[i] < C[i]) ? C[i] : T(0) for i in 1:n]
        @test result ≈ ref

        Z = similar(C)
        AppleAccelerate.venvlp!(Z, A, B, C)
        @test Z ≈ ref
    end
end

for T in (Float32, Float64)
    @testset "vDSP viclip::$T" begin
        X = T.(collect(-5.0:0.5:5.0))
        low = T(-2.0)
        high = T(2.0)

        # viclip: outside [low,high] pass through; inside: negative→low, non-negative→high
        result = AppleAccelerate.viclip(X, low, high)
        viclip_ref(x, lo, hi) = (x <= lo || hi <= x) ? x : (x < 0 ? lo : hi)
        expected = T[viclip_ref(x, low, high) for x in X]
        @test result ≈ expected

        Z = similar(X)
        AppleAccelerate.viclip!(Z, X, low, high)
        @test Z ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP vcmprs::$T" begin
        X = T[1, 2, 3, 4, 5, 6, 7, 8]
        gate = T[1, 0, 1, 0, 0, 1, 1, 0]

        result = AppleAccelerate.vcmprs(X, gate)
        n_kept = count(!iszero, gate)
        expected = X[gate .!= 0]
        @test result[1:n_kept] ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP vrampmul::$T" begin
        n = 100
        X = ones(T, n)
        start = T(1.0)
        step = T(0.5)

        result = AppleAccelerate.vrampmul(X, start, step)
        ramp = T[start + (i - 1) * step for i in 1:n]
        @test result ≈ X .* ramp

        Z = similar(X)
        AppleAccelerate.vrampmul!(Z, X, start, step)
        @test Z ≈ X .* ramp

        # Non-trivial input
        X2 = randn(T, n)
        result2 = AppleAccelerate.vrampmul(X2, start, step)
        @test result2 ≈ X2 .* ramp
    end
end

for T in (Float32, Float64)
    @testset "vDSP vrsum::$T" begin
        n = 50
        X = ones(T, n)
        scale = T(2.0)

        result = AppleAccelerate.vrsum(X, scale)
        # vrsum: result[0]=0; result[n] = result[n-1] + X[n]*scale  (0-based)
        ref = similar(X)
        ref[1] = T(0)
        for i in 2:n
            ref[i] = ref[i-1] + X[i] * scale
        end
        @test result ≈ ref

        Z = similar(X)
        AppleAccelerate.vrsum!(Z, X, scale)
        @test Z ≈ ref
    end
end

for T in (Float32, Float64)
    @testset "vDSP vsimps::$T" begin
        n = 101  # odd for Simpson's rule
        X = ones(T, n)
        step = T(0.1)

        result = AppleAccelerate.vsimps(X, step)
        # For constant f=1, Simpson integration gives cumulative area = n*step
        @test result[1] ≈ T(0) atol=eps(T)
        @test result[end] ≈ T((n - 1) * step) atol=T(0.01)

        Z = similar(X)
        AppleAccelerate.vsimps!(Z, X, step)
        @test Z ≈ result
    end
end

for T in (Float32, Float64)
    @testset "vDSP vlint::$T" begin
        # Ramp table: [0, 10, 20, 30, 40]
        table = T[0, 10, 20, 30, 40]
        # 0-based fractional indices
        indices = T[0.0, 0.5, 1.0, 2.5, 3.0]

        result = AppleAccelerate.vlint(table, indices)
        expected = T[0, 5, 10, 25, 30]
        @test result ≈ expected

        Z = Vector{T}(undef, length(indices))
        AppleAccelerate.vlint!(Z, table, indices)
        @test Z ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP vqint::$T" begin
        # Quadratic table: f(x) = x^2 → [0, 1, 4, 9, 16, 25]
        table = T[0, 1, 4, 9, 16, 25]
        # Fractional indices (must have p+2 < length(table))
        indices = T[0.0, 0.5, 1.0, 1.5, 2.0, 3.0]

        result = AppleAccelerate.vqint(table, indices)
        # For a quadratic table, quadratic interpolation should be exact
        expected = T[0, 0.25, 1, 2.25, 4, 9]
        @test result ≈ expected atol=T(0.01)

        Z = Vector{T}(undef, length(indices))
        AppleAccelerate.vqint!(Z, table, indices)
        @test Z ≈ expected atol=T(0.01)
    end
end

for T in (Float32, Float64)
    @testset "vDSP nzcros::$T" begin
        # Signal with known zero crossings at every consecutive pair
        X = T[1, -1, 1, -1, 1]
        (indices, cnt) = AppleAccelerate.nzcros(X)
        @test cnt == 4
        @test length(indices) == cnt

        # No crossings in a constant signal
        Y = ones(T, 10)
        (_, cnt2) = AppleAccelerate.nzcros(Y)
        @test cnt2 == 0
    end
end

# ============================================================
# Batch 1: Additional Compound Arithmetic
# ============================================================
for T in (Float32, Float64)
    @testset "vDSP Compound Arithmetic (Batch 1)::$T" begin
        A::Vector{T} = randn(N)
        B::Vector{T} = randn(N)
        C::Vector{T} = randn(N)
        D::Vector{T} = randn(N)
        R::Vector{T} = similar(A)

        # vma: A*B + C
        @test AppleAccelerate.vma(A, B, C) ≈ A .* B .+ C
        AppleAccelerate.vma!(R, A, B, C)
        @test R ≈ A .* B .+ C

        # vmsb: A*B - C
        @test AppleAccelerate.vmsb(A, B, C) ≈ A .* B .- C
        AppleAccelerate.vmsb!(R, A, B, C)
        @test R ≈ A .* B .- C

        # vmma: A*B + C*D
        @test AppleAccelerate.vmma(A, B, C, D) ≈ A .* B .+ C .* D
        AppleAccelerate.vmma!(R, A, B, C, D)
        @test R ≈ A .* B .+ C .* D

        # vmmsb: A*B - C*D
        @test AppleAccelerate.vmmsb(A, B, C, D) ≈ A .* B .- C .* D
        AppleAccelerate.vmmsb!(R, A, B, C, D)
        @test R ≈ A .* B .- C .* D

        # vmsa: A*B + c
        c::T = T(2.5)
        @test AppleAccelerate.vmsa(A, B, c) ≈ A .* B .+ c
        AppleAccelerate.vmsa!(R, A, B, c)
        @test R ≈ A .* B .+ c

        # vsmsb: A*b - C
        b::T = T(1.5)
        @test AppleAccelerate.vsmsb(A, b, C) ≈ A .* b .- C
        AppleAccelerate.vsmsb!(R, A, b, C)
        @test R ≈ A .* b .- C

        # vsmsma: A*b + C*d
        d::T = T(0.7)
        @test AppleAccelerate.vsmsma(A, b, C, d) ≈ A .* b .+ C .* d
        AppleAccelerate.vsmsma!(R, A, b, C, d)
        @test R ≈ A .* b .+ C .* d
    end
end

# ============================================================
# Batch 2: Extra Reductions
# ============================================================
for T in (Float32, Float64)
    @testset "vDSP Extra Reductions::$T" begin
        X::Vector{T} = randn(N)

        # rmsqv
        @test AppleAccelerate.rmsqv(X) ≈ sqrt(sum(X .^ 2) / length(X))

        # sve_svesq
        (s, ssq) = AppleAccelerate.sve_svesq(X)
        @test s ≈ sum(X)
        @test ssq ≈ sum(X .^ 2)

        # maxmgv / minmgv
        @test AppleAccelerate.maxmgv(X) ≈ maximum(abs.(X))
        @test AppleAccelerate.minmgv(X) ≈ minimum(abs.(X))

        # maxmgvi / minmgvi
        (val, idx) = AppleAccelerate.maxmgvi(X)
        absX = abs.(X)
        @test val ≈ maximum(absX)
        @test absX[idx] ≈ val

        (val2, idx2) = AppleAccelerate.minmgvi(X)
        @test val2 ≈ minimum(absX)
        @test absX[idx2] ≈ val2
    end
end

# ============================================================
# Batch 3: Vector Utility
# ============================================================
for T in (Float32, Float64)
    @testset "vDSP Vector Fill/Clear/Swap::$T" begin
        X = randn(T, N)
        Y = randn(T, N)
        X_orig = copy(X)
        Y_orig = copy(Y)

        # vclr!
        Z = randn(T, N)
        AppleAccelerate.vclr!(Z)
        @test all(iszero, Z)

        # vfill!
        AppleAccelerate.vfill!(Z, T(3.14))
        @test all(==(T(3.14)), Z)

        # vswap!
        A = copy(X_orig)
        B = copy(Y_orig)
        AppleAccelerate.vswap!(A, B)
        @test A ≈ Y_orig
        @test B ≈ X_orig
    end
end

for T in (Float32, Float64)
    @testset "vDSP Gather/Index::$T" begin
        A = T.(collect(1:10))
        B_idx = UInt[1, 3, 5, 7, 9]
        @test AppleAccelerate.vgathr(A, B_idx) ≈ T[1, 3, 5, 7, 9]

        # vindex uses 0-based float indices
        B_float = T[0.0, 2.0, 4.0, 6.0, 8.0]
        @test AppleAccelerate.vindex(A, B_float) ≈ T[1, 3, 5, 7, 9]
    end
end

for T in (Float32, Float64)
    @testset "vDSP Generation::$T" begin
        # vgen: linear ramp between two values
        result = AppleAccelerate.vgen(T(0), T(10), 11)
        expected = T[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        @test result ≈ expected
    end
end

for T in (Float32, Float64)
    @testset "vDSP Clipping Variants::$T" begin
        X = T.(collect(-5.0:0.5:5.0))

        # vclipc: clip with count
        low = T(-2.0)
        high = T(2.0)
        (clipped, nlow, nhigh) = AppleAccelerate.vclipc(X, low, high)
        @test clipped ≈ clamp.(X, low, high)
        @test nlow == count(x -> x < low, X)
        @test nhigh == count(x -> x > high, X)

        # vlim
        b = T(0.0)
        c = T(1.0)
        result_vlim = AppleAccelerate.vlim(X, b, c)
        expected_vlim = T[(x >= b) ? c : -c for x in X]
        @test result_vlim ≈ expected_vlim
    end
end

for T in (Float32, Float64)
    @testset "vDSP Sort::$T" begin
        X = randn(T, 100)

        # vsort! ascending
        Y = copy(X)
        AppleAccelerate.vsort!(Y, true)
        @test Y ≈ sort(X)

        # vsort! descending
        Y = copy(X)
        AppleAccelerate.vsort!(Y, false)
        @test Y ≈ sort(X, rev=true)

        # vsorti
        perm = AppleAccelerate.vsorti(X, true)
        @test X[perm] ≈ sort(X)

        # vsorti! directly, with the documented pre-filled 0-based identity
        # buffer. The returned indices are 0-based; +1 gives a Julia permutation.
        idx = UInt.(0:length(X)-1)
        ret = AppleAccelerate.vsorti!(idx, X, true)
        @test ret === idx                                 # mutates & returns same buffer
        @test sort(Int.(idx)) == collect(0:length(X)-1)   # still a valid permutation
        @test X[Int.(idx) .+ 1] ≈ sort(X)                 # and it actually sorts X

        # descending
        idx_desc = UInt.(0:length(X)-1)
        AppleAccelerate.vsorti!(idx_desc, X, false)
        @test sort(Int.(idx_desc)) == collect(0:length(X)-1)
        @test X[Int.(idx_desc) .+ 1] ≈ sort(X, rev=true)

        # vsorti! and the allocating vsorti must agree
        @test (Int.(idx) .+ 1) == AppleAccelerate.vsorti(X, true)
    end
end

@testset "vDSP Integer Operations" begin
    A = Int32[1, -2, 3, -4, 5]
    B = Int32[10, 20, 30, 40, 50]

    # vaddi
    @test AppleAccelerate.vaddi(A, B) == A .+ B

    # vabsi
    @test AppleAccelerate.vabsi(A) == abs.(A)

    # vfilli!
    C = Vector{Int32}(undef, 5)
    AppleAccelerate.vfilli!(C, Int32(42))
    @test all(==(Int32(42)), C)

    # veqvi (XNOR)
    @test AppleAccelerate.veqvi(A, B) == .~(A .⊻ B)
end

# ============================================================
# Batch 4: Matrix Operations
# ============================================================
for T in (Float32, Float64)
    @testset "vDSP Matrix Operations::$T" begin
        m, p, n = 4, 3, 5
        A = randn(T, m, p)
        B = randn(T, p, n)

        # mmul
        C = AppleAccelerate.mmul(A, B)
        @test C ≈ A * B atol=T(1e-4)

        C2 = Matrix{T}(undef, m, n)
        AppleAccelerate.mmul!(C2, A, B)
        @test C2 ≈ A * B atol=T(1e-4)

        # mtrans
        D = randn(T, 3, 4)
        Dt = AppleAccelerate.mtrans(D)
        @test Dt ≈ transpose(D) atol=T(1e-6)

        # mmov
        E = randn(T, 3, 4)
        F = AppleAccelerate.mmov(E)
        @test F ≈ E
    end
end

# ============================================================
# Batch 6: Type Conversion (int ↔ float)
# ============================================================
for T in (Float32, Float64)
    @testset "vDSP Type Conversion int↔float::$T" begin
        # float → int (truncating)
        X = T[1.7, -2.3, 3.9, -4.1, 0.5]
        @test AppleAccelerate.vfix8(X) == Int8.(trunc.(X))
        @test AppleAccelerate.vfix16(X) == Int16.(trunc.(X))
        @test AppleAccelerate.vfix32(X) == Int32.(trunc.(X))

        # float → unsigned int (truncating)
        Xpos = T[1.7, 2.3, 3.9, 4.1, 0.5]
        @test AppleAccelerate.vfixu8(Xpos) == UInt8.(trunc.(Xpos))
        @test AppleAccelerate.vfixu16(Xpos) == UInt16.(trunc.(Xpos))
        @test AppleAccelerate.vfixu32(Xpos) == UInt32.(trunc.(Xpos))

        # float → int (rounding)
        @test AppleAccelerate.vfixr8(X) == Int8.(round.(X))
        @test AppleAccelerate.vfixr16(X) == Int16.(round.(X))
        @test AppleAccelerate.vfixr32(X) == Int32.(round.(X))

        # float → unsigned int (rounding)
        @test AppleAccelerate.vfixru8(Xpos) == UInt8.(round.(Xpos))
        @test AppleAccelerate.vfixru16(Xpos) == UInt16.(round.(Xpos))
        @test AppleAccelerate.vfixru32(Xpos) == UInt32.(round.(Xpos))

        # int → float
        Ai8 = Int8[1, -2, 3, -4, 5]
        @test AppleAccelerate.vflt8(Ai8, T) ≈ T.(Ai8)

        Ai16 = Int16[100, -200, 300, -400, 500]
        @test AppleAccelerate.vflt16(Ai16, T) ≈ T.(Ai16)

        Ai32 = Int32[1000, -2000, 3000, -4000, 5000]
        @test AppleAccelerate.vflt32(Ai32, T) ≈ T.(Ai32)

        # unsigned int → float
        Au8 = UInt8[1, 2, 3, 4, 5]
        @test AppleAccelerate.vfltu8(Au8, T) ≈ T.(Au8)

        Au16 = UInt16[100, 200, 300, 400, 500]
        @test AppleAccelerate.vfltu16(Au16, T) ≈ T.(Au16)

        Au32 = UInt32[1000, 2000, 3000, 4000, 5000]
        @test AppleAccelerate.vfltu32(Au32, T) ≈ T.(Au32)
    end
end

# ============================================================
# Batch 8: Image Convolution
# ============================================================
# Reference zero-padded 2D convolution (correlation, as vDSP computes), used to
# validate the image-convolution wrappers on non-square inputs and asymmetric
# filters where any row/column mix-up would show up.
function naive_imgfir(A::Matrix, F::Matrix)
    nr, nc = size(A)
    fr, fc = size(F)
    pr, pc = fr ÷ 2, fc ÷ 2
    C = zeros(eltype(A), nr, nc)
    for i in 1:nr, j in 1:nc
        s = zero(eltype(A))
        for u in 1:fr, v in 1:fc
            ai = i + (u - 1 - pr)
            aj = j + (v - 1 - pc)
            (1 <= ai <= nr && 1 <= aj <= nc) || continue
            s += A[ai, aj] * F[u, v]
        end
        C[i, j] = s
    end
    return C
end

for T in (Float32, Float64)
    @testset "vDSP Image Convolution::$T" begin
        # Non-square input so a row/column swap can't pass unnoticed.
        A = randn(T, 11, 9)

        # f3x3 with an asymmetric kernel
        F3 = T[1 2 3; 4 5 6; 7 8 9]
        C = AppleAccelerate.f3x3(A, F3)
        ref3 = naive_imgfir(A, F3)
        @test C[2:10, 2:8] ≈ ref3[2:10, 2:8] atol=T(1e-4)

        # f5x5 with an asymmetric kernel
        F5 = T.(reshape(1:25, 5, 5))
        C5 = AppleAccelerate.f5x5(A, F5)
        ref5 = naive_imgfir(A, F5)
        @test C5[3:9, 3:7] ≈ ref5[3:9, 3:7] atol=T(1e-3)

        # imgfir with a non-square, asymmetric kernel (odd dims so the
        # centering is unambiguous while fr != fc still exercises the swap)
        F = T.(reshape(1:15, 3, 5))   # 3×5
        C_imgfir = AppleAccelerate.imgfir(A, F)
        ref = naive_imgfir(A, F)
        @test C_imgfir[2:10, 3:7] ≈ ref[2:10, 3:7] atol=T(1e-3)

        # imgfir with a 3×3 identity should reproduce the interior of A
        Fid = zeros(T, 3, 3); Fid[2, 2] = T(1)
        @test AppleAccelerate.imgfir(A, Fid)[2:10, 2:8] ≈ A[2:10, 2:8]
    end
end

# ============================================================
# Format conversion tests (ctoz / ztoc)
# ============================================================
for T in (Float32, Float64)
    @testset "ctoz/ztoc::$T" begin
        X = complex.(randn(T, N), randn(T, N))
        (re, im) = AppleAccelerate.ctoz(X)
        @test re ≈ real.(X)
        @test im ≈ imag.(X)

        Y = AppleAccelerate.ztoc(re, im)
        @test Y ≈ X
    end
end

# ============================================================
# svsub correctness + edge cases
# ============================================================
for T in (Float32, Float64)
    @testset "svsub correctness::$T" begin
        X::Vector{T} = randn(N)
        c::T = randn()

        # Allocating svsub computes c - X
        @test AppleAccelerate.svsub(X, c) ≈ (c .- X)

        # Mutating svsub! matches allocating result
        out = similar(X)
        AppleAccelerate.svsub!(out, X, c)
        @test out ≈ (c .- X)
        @test out ≈ AppleAccelerate.svsub(X, c)

        # In-place aliasing (out === X) works for a mutating op
        Z = copy(X)
        AppleAccelerate.vadd!(Z, Z, X)
        @test Z ≈ (X .+ X)

        # Single-element array
        x1 = T[3]
        @test AppleAccelerate.svsub(x1, T(5)) ≈ T[2]
        @test AppleAccelerate.vadd(x1, x1) ≈ T[6]
    end
end

# Regression: mutating multi-vector ops used to size the loop from one argument
# and pass the others as bare pointers, so mismatched lengths read out of bounds
# with no error. They must now throw DimensionMismatch.
for T in (Float32, Float64)
    @testset "length validation throws DimensionMismatch::$T" begin
        a5 = randn(T, 5); b5 = randn(T, 5); c5 = randn(T, 5); d5 = randn(T, 5)
        out5 = similar(a5)
        a3 = randn(T, 3)
        # Element-wise op (vadd!): mismatched input.
        @test_throws DimensionMismatch AppleAccelerate.vadd!(out5, a5, a3)
        # Element-wise op: mismatched output.
        @test_throws DimensionMismatch AppleAccelerate.vadd!(similar(a3), a5, b5)
        # Compound op (vma!): result = A*B + C.
        @test_throws DimensionMismatch AppleAccelerate.vma!(out5, a5, b5, a3)
        # 4-vector compound op (vmma!).
        @test_throws DimensionMismatch AppleAccelerate.vmma!(out5, a5, b5, c5, a3)
        # Reduction (dot).
        @test_throws DimensionMismatch AppleAccelerate.dot(a5, a3)
    end
end

# ============================================================
# Reference reconciliation: cross-validate every operation
# against its Base / Statistics / LinearAlgebra equivalent,
# including mutating (!) variants and ops not otherwise
# numerically reconciled above. (Coverage is already 100%;
# this testset hardens the confidence goal.)
# ============================================================
for T in (Float32, Float64)
    @testset "Reference reconciliation::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        Z::Vector{T} = similar(X)
        Xpos::Vector{T} = abs.(randn(N)) .+ T(0.5)
        Ypos::Vector{T} = abs.(randn(N)) .+ T(0.5)

        # --- 1-arg elementwise mutating variants vs Base ---
        for (af, bf, src) in (
                (AppleAccelerate.sin,  sin,  X), (AppleAccelerate.cos, cos, X),
                (AppleAccelerate.tan,  tan,  X), (AppleAccelerate.exp, exp, X),
                (AppleAccelerate.expm1, expm1, X), (AppleAccelerate.log, log, Xpos),
                (AppleAccelerate.log1p, log1p, Xpos), (AppleAccelerate.sqrt, sqrt, Xpos),
                (AppleAccelerate.cbrt, cbrt, X), (AppleAccelerate.abs, abs, X),
                (AppleAccelerate.atan, atan, X), (AppleAccelerate.floor, floor, X),
                (AppleAccelerate.ceil, ceil, X), (AppleAccelerate.round, round, X),
                (AppleAccelerate.trunc, trunc, X))
            af!  = getfield(AppleAccelerate, Symbol(string(nameof(af), "!")))
            @test af(src) ≈ bf.(src)
            af!(Z, src)
            @test Z ≈ bf.(src)
        end

        # asin/acos need |x|<=1
        W = T(2) .* rand(T, N) .- T(1)
        @test AppleAccelerate.asin(W) ≈ asin.(W)
        AppleAccelerate.asin!(Z, W); @test Z ≈ asin.(W)
        @test AppleAccelerate.acos(W) ≈ acos.(W)
        AppleAccelerate.acos!(Z, W); @test Z ≈ acos.(W)

        # --- 2-arg elementwise mutating variants vs Base ---
        @test AppleAccelerate.copysign(Xpos, Y) ≈ copysign.(Xpos, Y)
        AppleAccelerate.copysign!(Z, Xpos, Y); @test Z ≈ copysign.(Xpos, Y)
        @test AppleAccelerate.rem(X, Ypos) ≈ rem.(X, Ypos)
        AppleAccelerate.rem!(Z, X, Ypos); @test Z ≈ rem.(X, Ypos)
        @test AppleAccelerate.atan(X, Y) ≈ atan.(X, Y)
        AppleAccelerate.atan!(Z, X, Y); @test Z ≈ atan.(X, Y)
        @test AppleAccelerate.div_float(X, Ypos) ≈ X ./ Ypos
        AppleAccelerate.div_float!(Z, X, Ypos); @test Z ≈ X ./ Ypos

        # --- pow allocating + mutating vs Base .^ ---
        @test AppleAccelerate.pow(Xpos, Y) ≈ Xpos .^ Y
        AppleAccelerate.pow!(Z, Xpos, Y); @test Z ≈ Xpos .^ Y

        # --- reductions vs Base/Statistics/LinearAlgebra ---
        @test AppleAccelerate.maximum(X) ≈ maximum(X)
        @test AppleAccelerate.minimum(X) ≈ minimum(X)
        @test AppleAccelerate.sum(X)     ≈ sum(X)
        @test AppleAccelerate.mean(X)    ≈ mean(X)
        @test AppleAccelerate.findmax(X) == findmax(X)
        @test AppleAccelerate.findmin(X) == findmin(X)
        @test AppleAccelerate.dot(X, Y)  ≈ LinearAlgebra.dot(X, Y)
        @test AppleAccelerate.rmsqv(X)   ≈ sqrt(mean(abs2, X))

        # vnormalize stddev/mean vs Statistics (population, corrected=false)
        (norm, m, s) = AppleAccelerate.vnormalize(X)
        @test m ≈ mean(X) atol=T(1e-4)
        @test s ≈ std(X, corrected=false) atol=T(1e-3)

        # --- vector-scalar mutating variants vs Base ---
        c::T = randn()
        AppleAccelerate.vsadd!(Z, X, c); @test Z ≈ X .+ c
        AppleAccelerate.vssub!(Z, X, c); @test Z ≈ X .- c
        AppleAccelerate.svsub!(Z, X, c); @test Z ≈ c .- X
        AppleAccelerate.vsmul!(Z, X, c); @test Z ≈ X .* c
        cnz::T = c == 0 ? one(T) : c
        AppleAccelerate.vsdiv!(Z, X, cnz); @test Z ≈ X ./ cnz

        # --- threshold/limit/table ops vs explicit reference ---
        @test AppleAccelerate.vthrsc(X, T(0), T(1)) ≈ T[x >= 0 ? T(1) : T(-1) for x in X]
        # vtabi identity lookup reproduces the table
        tab = T.(collect(0:10:40))
        idxA = T.(collect(0:4))
        @test AppleAccelerate.vtabi(idxA, T(1), T(0), tab) ≈ tab
        # vgenp piecewise-linear matches manual linear interpolation
        bpos = T[0, 2, 4]; bval = T[0, 10, 30]
        gp = AppleAccelerate.vgenp(bval, bpos, 5)
        gpref = T[0, 5, 10, 20, 30]
        @test gp ≈ gpref
    end
end

# Validation / error paths reachable from public API
for T in (Float32, Float64)
    @testset "Validation paths::$T" begin
        a3 = randn(T, 3); b5 = randn(T, 5)
        # 2-arg elementwise allocating shape check
        @test_throws DimensionMismatch AppleAccelerate.copysign(a3, b5)
        @test_throws DimensionMismatch AppleAccelerate.atan(a3, b5)
        @test_throws DimensionMismatch AppleAccelerate.pow(a3, b5)
        @test_throws DimensionMismatch AppleAccelerate.remainder(a3, b5)
        # reductions
        @test_throws DimensionMismatch AppleAccelerate.distancesq(a3, b5)
        # two-vector vDSP op
        @test_throws DimensionMismatch AppleAccelerate.vmax!(similar(a3), a3, b5)
        # compound ops with scalar
        c = T(1)
        @test_throws DimensionMismatch AppleAccelerate.vasm!(similar(a3), a3, b5, c)
        @test_throws DimensionMismatch AppleAccelerate.vmsa!(similar(a3), a3, b5, c)
        @test_throws DimensionMismatch AppleAccelerate.vsma!(similar(a3), a3, c, b5)
        @test_throws DimensionMismatch AppleAccelerate.vintb!(similar(a3), a3, b5, c)
        # matrix shape checks
        @test_throws DimensionMismatch AppleAccelerate.mmul(randn(T, 3, 2), randn(T, 3, 2))
        @test_throws DimensionMismatch AppleAccelerate.f3x3(randn(T, 5, 5), randn(T, 2, 2))
        @test_throws DimensionMismatch AppleAccelerate.f5x5(randn(T, 7, 7), randn(T, 3, 3))
    end
end

# Out-of-bounds-write guards on single-input / interpolation / conversion
# mutating `!` methods: undersized `result`/`C` buffers must throw before the
# C routine writes past the end of the buffer.
for T in (Float32, Float64)
    @testset "Bounds-check guards::$T" begin
        X = T[1, 2, 3, 4, 5]
        small = Vector{T}(undef, 2)   # deliberately too small
        ok = similar(X)

        # --- sliding window (A2): window range + result length ---
        # valid call still works
        @test AppleAccelerate.vswsum(X, 2) ≈ [X[i] + X[i+1] for i in 1:length(X)-1]
        @test AppleAccelerate.vswmax(X, 2) ≈ [max(X[i], X[i+1]) for i in 1:length(X)-1]
        # invalid window (allocating + mutating, sum + max)
        @test_throws ArgumentError AppleAccelerate.vswsum(X, 0)
        @test_throws ArgumentError AppleAccelerate.vswsum(X, length(X) + 1)
        @test_throws ArgumentError AppleAccelerate.vswmax(X, 0)
        @test_throws ArgumentError AppleAccelerate.vswmax(X, length(X) + 1)
        @test_throws ArgumentError AppleAccelerate.vswsum!(Vector{T}(undef, 4), X, 0)
        @test_throws ArgumentError AppleAccelerate.vswmax!(Vector{T}(undef, 4), X, 0)
        # undersized result for a valid window (n_out = 4)
        @test_throws DimensionMismatch AppleAccelerate.vswsum!(small, X, 2)
        @test_throws DimensionMismatch AppleAccelerate.vswmax!(small, X, 2)

        # --- clipping / thresholding (A1) ---
        lo = T(2); hi = T(4); thr = T(3)
        @test AppleAccelerate.vclip!(similar(X), X, lo, hi) ≈ clamp.(X, lo, hi)
        @test_throws DimensionMismatch AppleAccelerate.vclip!(small, X, lo, hi)
        @test_throws DimensionMismatch AppleAccelerate.viclip!(small, X, lo, hi)
        @test_throws DimensionMismatch AppleAccelerate.vthr!(small, X, thr)
        @test_throws DimensionMismatch AppleAccelerate.vthres!(small, X, thr)
        @test_throws DimensionMismatch AppleAccelerate.vlim!(small, X, T(0), T(1))
        @test_throws DimensionMismatch AppleAccelerate.vthrsc!(small, X, T(0), T(1))
        @test_throws DimensionMismatch AppleAccelerate.vclipc!(small, X, lo, hi)

        # --- interpolation (A1): result vs length(indices) ---
        table = T[10, 20, 30, 40]
        idx = T[0.0, 0.5, 1.5, 2.5]  # 4 fractional indices
        @test length(AppleAccelerate.vlint(table, idx)) == length(idx)
        @test_throws DimensionMismatch AppleAccelerate.vlint!(small, table, idx)
        @test_throws DimensionMismatch AppleAccelerate.vqint!(small, table, idx)

        # --- gather / index (A1 + A3): result length + index bounds ---
        Bidx = UInt[1, 3, 5]
        @test AppleAccelerate.vgathr(X, Bidx) ≈ T[1, 3, 5]
        @test_throws DimensionMismatch AppleAccelerate.vgathr!(Vector{T}(undef, 1), X, Bidx)
        @test_throws DimensionMismatch AppleAccelerate.vindex!(Vector{T}(undef, 1), X, T[0, 2, 4])
        # A3: out-of-range 1-based index (0 and length+1) must throw
        @test_throws BoundsError AppleAccelerate.vgathr(X, UInt[1, UInt(length(X) + 1)])
        @test_throws BoundsError AppleAccelerate.vgathr(X, UInt[0])

        # --- table lookup (A1): D vs length(A) ---
        tab = T[1, 2, 3, 4]
        idxA = T[0, 1, 2, 3]
        @test AppleAccelerate.vtabi(idxA, T(1), T(0), tab) ≈ tab
        @test_throws DimensionMismatch AppleAccelerate.vtabi!(Vector{T}(undef, 1), idxA, T(1), T(0), tab)

        # --- piecewise generation (A1): C vs n ---
        @test_throws DimensionMismatch AppleAccelerate.vgenp!(Vector{T}(undef, 2), T[1, 2], T[0, 10], 5)

        # --- type conversion (A1): C vs length(A) ---
        Xf = T[1.2, 2.7, 3.1, 4.9]
        @test AppleAccelerate.vfix32(Xf) == Int32.(trunc.(Xf))
        @test_throws DimensionMismatch AppleAccelerate.vfix32!(Vector{Int32}(undef, 1), Xf)
        @test_throws DimensionMismatch AppleAccelerate.vfixu8!(Vector{UInt8}(undef, 1), Xf)
        @test_throws DimensionMismatch AppleAccelerate.vfixr16!(Vector{Int16}(undef, 1), Xf)
        @test_throws DimensionMismatch AppleAccelerate.vfixru32!(Vector{UInt32}(undef, 1), Xf)
        Ai = Int32[1, 2, 3, 4]
        @test AppleAccelerate.vflt32(Ai, T) ≈ T.(Ai)
        @test_throws DimensionMismatch AppleAccelerate.vflt32!(Vector{T}(undef, 1), Ai)
        @test_throws DimensionMismatch AppleAccelerate.vfltu8!(Vector{T}(undef, 1), UInt8[1, 2, 3, 4])

        # --- unary ops: undersized result must throw ---
        @test_throws DimensionMismatch AppleAccelerate.vneg!(small, X)
        @test_throws DimensionMismatch AppleAccelerate.vnabs!(small, X)
        @test_throws DimensionMismatch AppleAccelerate.vsq!(small, X)
        @test_throws DimensionMismatch AppleAccelerate.vssq!(small, X)
        @test_throws DimensionMismatch AppleAccelerate.vfrac!(small, X)
        @test_throws DimensionMismatch AppleAccelerate.vabs!(small, X)

        # --- vector-scalar ops: result and X must match (count is taken
        #     from result, so an oversized result would read past X) ---
        big = Vector{T}(undef, length(X) + 3)
        c = T(2)
        @test_throws DimensionMismatch AppleAccelerate.vsadd!(small, X, c)
        @test_throws DimensionMismatch AppleAccelerate.vsadd!(big, X, c)
        @test_throws DimensionMismatch AppleAccelerate.vsdiv!(big, X, c)
        @test_throws DimensionMismatch AppleAccelerate.vsmul!(big, X, c)
        @test_throws DimensionMismatch AppleAccelerate.vssub!(big, X, c)
        @test_throws DimensionMismatch AppleAccelerate.svsub!(big, X, c)

        # --- scalar / vector divide ---
        @test_throws DimensionMismatch AppleAccelerate.svdiv!(small, X, c)

        # --- tapered merge: all three lengths must match ---
        @test_throws DimensionMismatch AppleAccelerate.vtmerg!(small, X, X)
        @test_throws DimensionMismatch AppleAccelerate.vtmerg!(ok, X, T[1, 2])

        # --- ramp multiply ---
        @test_throws DimensionMismatch AppleAccelerate.vrampmul!(small, X, T(0), T(1))

        # --- integration / running ops ---
        @test_throws DimensionMismatch AppleAccelerate.vrsum!(small, X, T(1))
        @test_throws DimensionMismatch AppleAccelerate.vsimps!(small, X, T(1))
        @test_throws DimensionMismatch AppleAccelerate.vtrapz!(small, X, T(1))

        # --- polynomial evaluation ---
        @test_throws DimensionMismatch AppleAccelerate.vpoly!(small, T[1, 0], X)

        # --- normalization ---
        @test_throws DimensionMismatch AppleAccelerate.vnormalize!(small, X)
    end
end

# Int32 mutating ops: mismatched buffers must throw before the C routine
# reads or writes out of bounds.
@testset "Bounds-check guards::Int32" begin
    Ai = Int32[1, -2, 3, -4, 5]
    Bi = Int32[1, 2, 3, 4, 5]
    smalli = Vector{Int32}(undef, 2)
    @test_throws DimensionMismatch AppleAccelerate.vaddi!(smalli, Ai, Bi)
    @test_throws DimensionMismatch AppleAccelerate.vaddi!(similar(Ai), Ai, Int32[1, 2])
    @test_throws DimensionMismatch AppleAccelerate.vabsi!(smalli, Ai)
    @test_throws DimensionMismatch AppleAccelerate.veqvi!(smalli, Ai, Bi)
    @test_throws DimensionMismatch AppleAccelerate.veqvi!(similar(Ai), Ai, Int32[1, 2])
end

end # @testset "Array Operations"
