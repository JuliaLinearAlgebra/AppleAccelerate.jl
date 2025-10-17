using LinearAlgebra
using AppleAccelerate
using DSP, Test, Random, Statistics

if !Sys.isapple()
    @info("AppleAccelerate.jl will be tested only on macOS. Exiting.")
    exit(0)
end

Random.seed!(7)
N = 1_000

@testset "AppleAccelerate.jl" begin
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


@testset "DCT::Float32" begin
    r=rand(Float32,2^16)
    d1=DSP.dct(r)
    plan_accel = AppleAccelerate.plan_dct(length(r), 2)
    d2=AppleAccelerate.dct(r, plan_accel)
    @test norm(d1[2]/d2[2]*d2[2:end]-d1[2:end])≤1000eps(Float32)
end


for T in (Float32,  Float64)
    @testset "Convolution & Correlation::$T" begin
        X::Vector{T} = randn(N)
        Y::Vector{T} = randn(N)
        @testset "Testing $f::$T" for f in [:conv, :xcorr]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X, Y) ≈ fb(X, Y)
        end

        @testset "Testing $f::$T" for f in [:xcorr]
            @eval fb = $f
            @eval fa = AppleAccelerate.$f
            @test fa(X) ≈ fb(X, copy(X))
        end
    end
end

for T in (Float64, )
    @testset "Biquadratic Flitering::$T" begin
        @testset "Single Section::$T" begin
            X::Vector{T} = randn(10)
            d::Vector{T} = zeros(4)
            c::Vector{T} = [x%0.5 for x in randn(5)]
            fdsp = DSP.Biquad(c[1], c[2], c[3], c[4], c[5])
            fa = AppleAccelerate.biquadcreate(c, 1)
            @test DSP.filt(fdsp, X) ≈ AppleAccelerate.biquad(X, d, length(X), fa)
        end
    end
end


for T in (Float32, Float64)
    @testset "Window Functions::$T" begin
        N = 64
        @testset "Testing blackman::$T" begin
            Wa = AppleAccelerate.blackman(N, T)
            Wb = 0.42 .- (0.5cos.(2pi.*(0:(N-1))./N)) .+ (0.08cos.(4pi.*(0:(N-1))./N))
            @test Wa ≈ Wb
        end

        @testset "Testing hamming::$T" begin
            Wa = AppleAccelerate.hamming(N, T)
            Wb = 0.54 .- 0.46cos.(2pi.*(0:(N-1))./N)
            @test Wa ≈ Wb
        end

        @testset "Testing hanning::$T" begin
            Wa = AppleAccelerate.hanning(N, T)
            Wb = 0.5(1 .- cos.(2pi.*(0:(N-1))./N))
            @test Wa ≈ Wb
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

#=
AppleAccelerate.@replaceBase(sin, atan, /)

@testset "Replace Base::$T" for T in (Float32, Float64)
X::Array{T} = randn(N)
Y::Array{T} = abs.(randn(N))

@test Base.sin.(X) == AppleAccelerate.sin(X)
@test Base.atan.(X, Y) == AppleAccelerate.atan(X, Y)
@test X ./ Y  == AppleAccelerate.div_float(X, Y)
end
=#
end

if AppleAccelerate.get_macos_version() < v"13.4"
    @info("AppleAccelerate.jl needs macOS >= 13.4 for BLAS forwarding. Not testing forwarding capabilities.")
    exit(0)
end

# Set up a debugging fallback function that prints out a stacktrace if the LinearAlgebra
# tests end up calling a function that we don't have forwarded.
accelerate_header = """
import LinearAlgebra
function debug_missing_function()
    println("Missing BLAS/LAPACK function!")
    display(stacktrace())
end
LinearAlgebra.BLAS.lbt_set_default_func(@cfunction(debug_missing_function, Cvoid, ()))
import AppleAccelerate
"""
eval(Meta.parseall(accelerate_header))

config = BLAS.get_config()
@info("Running with $(length(config.loaded_libs)) libraries loaded:")
display(config.loaded_libs)

using Test
@testset "Accelerate Forwarding Sanity Tests" begin
    ver = AppleAccelerate.get_macos_version()
    if ver >= v"13.4"
        @test LinearAlgebra.peakflops() > 0
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :lp64).libname, "Accelerate")
        @test endswith(BLAS.lbt_find_backing_library("dgemm_", :ilp64).libname, "Accelerate")

        # Accelerate has `_rook` symbols:
        @test endswith(BLAS.lbt_find_backing_library("dsytrf_rook_", :ilp64).libname, "Accelerate")
    end
end

@testset "CBLAS dot test" begin
    a = ComplexF64[
        1 + 1im,
        2 - 2im,
        3 + 3im
    ]
    @test BLAS.dotc(a, a) ≈ ComplexF64(28)
    @test BLAS.dotu(a, a) ≈ ComplexF64(12im)

    a = Float32[1, 2, 3]
    @test BLAS.dot(a, a) ≈ 14f0
end

@testset "BLAS threading tests" begin
    if AppleAccelerate.get_macos_version() >= v"15"
        AppleAccelerate.set_num_threads(1)
        @test AppleAccelerate.get_num_threads() == 1
        AppleAccelerate.set_num_threads(4)
        @test AppleAccelerate.get_num_threads() == Sys.CPU_THREADS
    else
        @test AppleAccelerate.get_num_threads() == -1
    end
end

linalg_stdlib_test_path = joinpath(dirname(pathof(LinearAlgebra)), "..", "test")

# Don't run blas.jl tests since the "strided interface blas" tests are currently failing
#=
@testset verbose=true "LinearAlgebra.jl BLAS tests" begin
    joinpath(linalg_stdlib_test_path, "blas.jl") |> include
end
=#

@testset verbose=true "LinearAlgebra.jl LAPACK tests" begin
    joinpath(linalg_stdlib_test_path, "lapack.jl") |> include
end
