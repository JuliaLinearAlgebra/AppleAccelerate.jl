@testset "Complex Array Operations (Batch 5)" begin

for T in (Float32, Float64)
    @testset "Complex-Complex Binary::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # zvadd
        @test AppleAccelerate.zvadd(A, B) ≈ A .+ B
        R = similar(A)
        AppleAccelerate.zvadd!(R, A, B)
        @test R ≈ A .+ B

        # zvsub
        @test AppleAccelerate.zvsub(A, B) ≈ A .- B
        AppleAccelerate.zvsub!(R, A, B)
        @test R ≈ A .- B

        # zvcmul: conj(A) * B
        @test AppleAccelerate.zvcmul(A, B) ≈ conj.(A) .* B
        AppleAccelerate.zvcmul!(R, A, B)
        @test R ≈ conj.(A) .* B
    end
end

for T in (Float32, Float64)
    @testset "Complex-Real Binary::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{T} = randn(T, N)

        # zrvmul: complex * real
        @test AppleAccelerate.zrvmul(A, B) ≈ A .* B
        R = similar(A)
        AppleAccelerate.zrvmul!(R, A, B)
        @test R ≈ A .* B

        # zrvadd: complex + real (adds to real part)
        expected_add = complex.(real.(A) .+ B, imag.(A))
        @test AppleAccelerate.zrvadd(A, B) ≈ expected_add
        AppleAccelerate.zrvadd!(R, A, B)
        @test R ≈ expected_add

        # zrvsub: complex - real (subtracts from real part)
        expected_sub = complex.(real.(A) .- B, imag.(A))
        @test AppleAccelerate.zrvsub(A, B) ≈ expected_sub
        AppleAccelerate.zrvsub!(R, A, B)
        @test R ≈ expected_sub

        # zrvdiv: complex / real (elementwise)
        Bnz::Vector{T} = B .+ T(2)  # ensure nonzero divisor
        expected_div = A ./ Bnz
        @test AppleAccelerate.zrvdiv(A, Bnz) ≈ expected_div
        AppleAccelerate.zrvdiv!(R, A, Bnz)
        @test R ≈ expected_div
    end
end

for T in (Float32, Float64)
    @testset "Complex Compound::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        C::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # zvcma: conj(A)*B + C
        @test AppleAccelerate.zvcma(A, B, C) ≈ conj.(A) .* B .+ C
        R = similar(A)
        AppleAccelerate.zvcma!(R, A, B, C)
        @test R ≈ conj.(A) .* B .+ C

        # zvma: A*B + C
        @test AppleAccelerate.zvma(A, B, C) ≈ A .* B .+ C
        AppleAccelerate.zvma!(R, A, B, C)
        @test R ≈ A .* B .+ C

        # zvsma: A*b + C (b is complex scalar)
        b::Complex{T} = complex(T(2.5), T(-1.3))
        @test AppleAccelerate.zvsma(A, b, C) ≈ A .* b .+ C
        AppleAccelerate.zvsma!(R, A, b, C)
        @test R ≈ A .* b .+ C
    end
end

for T in (Float32, Float64)
    @testset "Complex Dot Products::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))

        # zidotpr: conjugate dot — sum(conj(A) .* B)
        @test AppleAccelerate.zidotpr(A, B) ≈ sum(conj.(A) .* B)

        # zrdotpr: complex-real dot
        Br::Vector{T} = randn(T, N)
        @test AppleAccelerate.zrdotpr(A, Br) ≈ sum(A .* Br)
    end
end

for T in (Float32, Float64)
    @testset "Complex Fill::$T" begin
        c = complex(T(3.0), T(-1.5))
        R = Vector{Complex{T}}(undef, 100)
        AppleAccelerate.zvfill!(R, c)
        @test all(==(c), R)
    end
end

for T in (Float32, Float64)
    @testset "Complex Convolution::$T" begin
        X = complex.(randn(T, 50), randn(T, 50))
        K = complex.(randn(T, 5), randn(T, 5))
        result = AppleAccelerate.zconv(X, K)
        @test length(result) == length(X) + length(K) - 1
        @test eltype(result) == Complex{T}
        # Numerical correctness: vDSP_zconv applies the kernel without reversal
        # (correlation-style), matching DSP.conv with a reversed kernel.
        @test result ≈ DSP.conv(X, reverse(K)) atol=T(1e-3)
    end
end

for T in (Float32, Float64)
    @testset "Complex DimensionMismatch::$T" begin
        A::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        B::Vector{Complex{T}} = complex.(randn(T, N), randn(T, N))
        Rbad = similar(A, N - 1)
        @test_throws DimensionMismatch AppleAccelerate.zvadd!(Rbad, A, B)
        @test_throws DimensionMismatch AppleAccelerate.zvsub!(Rbad, A, B)
    end
end

for T in (Float32, Float64)
    @testset "Complex Matrix Multiply::$T" begin
        m, p, n = 4, 3, 5
        A = complex.(randn(T, m, p), randn(T, m, p))
        B = complex.(randn(T, p, n), randn(T, p, n))

        C = AppleAccelerate.zmmul(A, B)
        @test C ≈ A * B atol=T(1e-3)

        C2 = Matrix{Complex{T}}(undef, m, n)
        AppleAccelerate.zmmul!(C2, A, B)
        @test C2 ≈ A * B atol=T(1e-3)
    end
end

end # @testset
