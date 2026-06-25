using AppleAccelerate
using AppleAccelerate: AASparseMatrix, AAFactorization, refactor!, solve, factor!,
                       SparseFactorizationCholesky
using SparseArrays
using LinearAlgebra
using Test

# Helpers to build well-conditioned test matrices with Int64-indexed CSC.
_spd(::Type{T}, n) where {T<:Real} = begin
    M = sprandn(n, n, 0.4)
    SparseMatrixCSC{T,Int64}(sparse(T.(M * M') + n * I))
end
_spd(::Type{T}, n) where {T<:Complex} = begin
    M = sprandn(T, n, n, 0.4)
    SparseMatrixCSC{T,Int64}(sparse(M * M' + n * I))
end
_gen(::Type{T}, n) where {T} = SparseMatrixCSC{T,Int64}(sparse(T.(sprandn(n, n, 0.5)) + n * I))

@testset "SparseArrays extension" begin
    @testset "extension is loaded" begin
        # The extension provides SparseMatrixCSC(::AASparseMatrix); its presence
        # confirms the ext compiled & loaded.
        @test hasmethod(SparseMatrixCSC, Tuple{AASparseMatrix})
    end

    @testset "backslash $(T)" for T in (Float64, Float32, ComplexF64, ComplexF32)
        n = 8
        A = _spd(T, n)
        b = T <: Complex ? randn(T, n) : randn(T, n)
        x = A \ b
        ref = Matrix(A) \ b
        @test x ≈ ref rtol=sqrt(eps(real(T))) * 50
        # multiple RHS
        B = T <: Complex ? randn(T, n, 3) : randn(T, n, 3)
        X = A \ B
        @test X ≈ Matrix(A) \ B rtol=sqrt(eps(real(T))) * 50
    end

    @testset "factorization entry points $(T)" for T in (Float64, Float32)
        n = 8
        A = _spd(T, n)
        b = randn(T, n)
        ref = Matrix(A) \ b
        for f in (cholesky(A), ldlt(A), factorize(A))
            @test solve(f, b) ≈ ref rtol=sqrt(eps(T)) * 50
        end
        # qr & lu on a general (nonsymmetric) matrix
        G = _gen(T, n)
        gref = Matrix(G) \ b
        @test solve(qr(G), b) ≈ gref rtol=sqrt(eps(T)) * 50
        if something(AppleAccelerate.get_macos_version(), v"0") >= v"15.5"
            @test solve(lu(G), b) ≈ gref rtol=sqrt(eps(T)) * 50
        end
    end

    @testset "qr least squares (rectangular)" begin
        m, n = 10, 6
        A = SparseMatrixCSC{Float64,Int64}(sparse(sprandn(m, n, 0.5) + sparse(1:n, 1:n, fill(3.0, n), m, n)))
        b = randn(m)
        x = A \ b
        @test x ≈ Matrix(A) \ b rtol=1e-8
    end

    @testset "SparseMatrixCSC <- AASparseMatrix roundtrip" begin
        # general
        G = _gen(Float64, 7)
        @test SparseMatrixCSC(AASparseMatrix(G)) == G
        # symmetric (stored as one triangle internally)
        S = _spd(Float64, 7)
        @test Matrix(SparseMatrixCSC(AASparseMatrix(S))) ≈ Matrix(S)
        # complex Hermitian
        H = _spd(ComplexF64, 5)
        @test Matrix(SparseMatrixCSC(AASparseMatrix(H))) ≈ Matrix(H)
        # transpose view
        Gt = transpose(AASparseMatrix(G))
        @test Matrix(SparseMatrixCSC(Gt)) ≈ Matrix(G)'
    end
end

@testset "refactor! (numeric refactorization)" begin
    @testset "Cholesky refactor $(T)" for T in (Float64, Float32)
        n = 8
        A = _spd(T, n)
        b = randn(T, n)
        f = AAFactorization(A)
        factor!(f, SparseFactorizationCholesky)
        @test solve(f, b) ≈ Matrix(A) \ b rtol=sqrt(eps(T)) * 50

        # New matrix, same sparsity pattern, scaled values.
        A2 = copy(A)
        A2.nzval .= A.nzval .* T(2)
        out = refactor!(f, A2)
        @test out === f
        @test solve(f, b) ≈ Matrix(A2) \ b rtol=sqrt(eps(T)) * 50
    end

    @testset "QR refactor (general)" begin
        n = 8
        G = _gen(Float64, n)
        b = randn(n)
        f = AAFactorization(G)
        factor!(f, AppleAccelerate.SparseFactorizationQR)
        @test solve(f, b) ≈ Matrix(G) \ b rtol=1e-8
        G2 = copy(G); G2.nzval .= G.nzval .* 3.0
        refactor!(f, G2)
        @test solve(f, b) ≈ Matrix(G2) \ b rtol=1e-8
    end

    @testset "refactor! input guards" begin
        A = _spd(Float64, 6)
        b = randn(6)
        # Not yet factored -> error.
        f0 = AAFactorization(A)
        @test_throws ArgumentError refactor!(f0, A)
        # Factor, then mismatched dimensions.
        f = AAFactorization(A); factor!(f)
        Bbig = _spd(Float64, 8)
        @test_throws DimensionMismatch refactor!(f, Bbig)
        # Same size, different nnz -> error.
        Adense = SparseMatrixCSC{Float64,Int64}(sparse(Matrix(A) .+ 1.0))
        if nnz(Adense) != nnz(A)
            @test_throws ArgumentError refactor!(f, Adense)
        end
    end
end
