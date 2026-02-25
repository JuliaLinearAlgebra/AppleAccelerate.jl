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
        @testset "Testing $f::$T" for f in [:sin,:sinpi,:cos,:cospi,:tan,:atan] # tanpi not defined in Base
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
            @test fa(X)[2] ≈ fb(X)[2]
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

        # dot: unconjugated dot product — sum(X .* Y)
        @test AppleAccelerate.dot(X, Y) ≈ sum(X .* Y)
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

end # @testset "Array Operations"
