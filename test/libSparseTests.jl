using AppleAccelerate
using SparseArrays
using LinearAlgebra
@testset "libSparse" begin
    @testset "wrappers" begin
        @testset "SparseMultiply and SparseMultiplyAdd" begin
            # less copy-paste heavy way via meta programming features?
            for T in (Float32, Float64)
                @eval begin
                    dense = rand($T, 3,3)
                    denseV = rand($T, 3)
                    dense2 = zeros($T, 3,3)
                    denseV2 = zeros($T, 3)
                    sparseM = sprand($T, 3, 3, 0.5)
                    sparse_data = sparseM.nzval
                    col_inds =  Clong.(sparseM.colptr .+ -1)
                    row_inds = Cint.(sparseM.rowval .+ -1)
                    GC.@preserve col_inds row_inds sparse_data begin
                        s = AppleAccelerate.SparseMatrixStructure(3, 3,
                            pointer(col_inds), pointer(row_inds),
                            AppleAccelerate.ATT_ORDINARY, 1
                        )
                        sparse_matrix = AppleAccelerate.SparseMatrix{$T}(s, pointer(sparse_data))
                        AppleAccelerate.SparseMultiply(sparse_matrix, dense, dense2)
                        @test dense2 ≈ sparseM * dense
                        AppleAccelerate.SparseMultiply(sparse_matrix, denseV, denseV2)
                        @test denseV2 ≈ sparseM * denseV
                        scalar = rand($T)
                        AppleAccelerate.SparseMultiply(scalar, sparse_matrix, dense, dense2)
                        @test dense2 ≈ scalar * sparseM * dense
                        AppleAccelerate.SparseMultiply(scalar, sparse_matrix, denseV, denseV2)
                        @test denseV2 ≈ scalar * sparseM * denseV
                        fill!(dense2, 1.0)
                        fill!(denseV2, 1.0)
                        AppleAccelerate.SparseMultiplyAdd(sparse_matrix, dense, dense2)
                        @test dense2≈ ones(size(dense2))+sparseM * dense
                        AppleAccelerate.SparseMultiplyAdd(sparse_matrix, denseV, denseV2)
                        @test denseV2 ≈ ones(size(denseV2))+sparseM * denseV
                        fill!(dense2, 1.0)
                        fill!(denseV2, 1.0)
                        scalar = rand($T)
                        AppleAccelerate.SparseMultiplyAdd(scalar, sparse_matrix, dense, dense2)
                        @test dense2 ≈ ones(size(dense2))+ scalar * sparseM * dense
                        AppleAccelerate.SparseMultiplyAdd(scalar, sparse_matrix, denseV, denseV2)
                        @test denseV2 ≈ ones(size(denseV2))+ scalar * sparseM * denseV
                    end
                end
            end
        end
        @testset "cconvert and unsafe_convert dense" begin
            dense = rand(3,3)
            denseV = rand(3)
            dense2 = zeros(3, 3)
            denseV2 = zeros(3)
            sparseM = sprand(3, 3, 0.5)
            sparse_data = sparseM.nzval
            col_inds =  Clong.(sparseM.colptr .+ -1)
            row_inds = Cint.(sparseM.rowval .+ -1)
            # for lack of a way to call cconvert directly, I'll use the matrix multiply routines.
            GC.@preserve sparse_data col_inds row_inds begin
                s = AppleAccelerate.SparseMatrixStructure(3, 3,
                        pointer(col_inds), pointer(row_inds),
                        AppleAccelerate.ATT_ORDINARY, 1)
                sparse_matrix = AppleAccelerate.SparseMatrix{Cdouble}(s, pointer(sparse_data))
                AppleAccelerate.SparseMultiply(sparse_matrix, dense, dense2)
                @test dense2 ≈ sparseM * dense
                AppleAccelerate.SparseMultiply(sparse_matrix, denseV, denseV2)
                @test denseV2 ≈ sparseM * denseV
            end
        end
        @testset "transpose" begin
            for T in (Float32, Float64)
                @eval begin
                    sparseM = sprand($T, 3, 4, 0.5)
                    sparse_data = sparseM.nzval
                    col_inds =  Clong.(sparseM.colptr .+ -1)
                    row_inds = Cint.(sparseM.rowval .+ -1)
                    GC.@preserve sparse_data col_inds row_inds begin
                        s = AppleAccelerate.SparseMatrixStructure(3, 4,
                                pointer(col_inds), pointer(row_inds),
                                AppleAccelerate.ATT_ORDINARY, 1)
                        sparse_matrix = AppleAccelerate.SparseMatrix{$T}(s,
                                                pointer(sparse_data))
                        ATT_TRANSP = AppleAccelerate.ATT_TRANSPOSE
                        res = AppleAccelerate.SparseGetTranspose(sparse_matrix)
                        @test res.data == sparse_matrix.data &&
                                    (res.structure.attributes & ATT_TRANSP != zero(typeof(ATT_TRANSP)))
                        original = AppleAccelerate.SparseGetTranspose(res)
                        @test original.data == sparse_matrix.data &&
                                (original.structure.attributes & ATT_TRANSP == zero(typeof(ATT_TRANSP)))
                    end
                end
            end
        end
        @testset "SparseFactor" begin
            for T in (Float32, Float64)
                @eval begin
                    sparseM = sprand($T, 4, 4, 0.5)
                    while det(sparseM) == 0
                        global sparseM = sprand($T, 4, 4, 0.5)
                    end
                    sparse_data = copy(sparseM.nzval)
                    col_inds = Clong.(sparseM.colptr .+ -1)
                    row_inds = Cint.(sparseM.rowval .+ -1)
                    expected = rand($T, 4, 4)
                    B = Array(sparseM) * expected
                    GC.@preserve sparse_data col_inds row_inds begin
                        s = AppleAccelerate.SparseMatrixStructure(4, 4,
                            pointer(col_inds), pointer(row_inds),
                            AppleAccelerate.ATT_ORDINARY, 1)
                        sparse_matrix = AppleAccelerate.SparseMatrix{$T}(s,
                                                pointer(sparse_data))
                        qrType = AppleAccelerate.SparseFactorizationQR
                        sf = AppleAccelerate.SparseFactor(qrType, s)
                        @test sf.status == AppleAccelerate.SparseStatusOk
                        f = AppleAccelerate.SparseFactor(qrType, sparse_matrix)
                        AppleAccelerate.SparseSolve(f, B)
                        @test B ≈ expected
                        AppleAccelerate.SparseCleanup(f)
                        AppleAccelerate.SparseCleanup(sf)
                    end
                end
            end
        end
    end
    @testset "AASparseMatrix" begin
        @testset "arithmetic" begin
            N = 10
            sparseM = sprand(N, N, 0.3)
            alpha = rand(Float64)
            aaM = AASparseMatrix(sparseM)
            X, x = rand(N, N), rand(N)
            @test aaM*x ≈ sparseM*x
            @test aaM*X ≈ sparseM*X
            @test alpha*aaM*x ≈ alpha*sparseM*x
            @test alpha*aaM*X ≈ alpha*sparseM*X
            Y, y = rand(N,N), rand(N)
            FMA, fma = sparseM*X + Y, sparseM*x + y
            muladd!(aaM, x, y)
            muladd!(aaM, X, Y)
            @test Y ≈ FMA
            @test y ≈ fma
            Y, y = rand(N,N), rand(N)
            FMA2, fma2 = alpha*sparseM*X + Y, alpha*sparseM*x + y
            muladd!(alpha, aaM, x, y)
            muladd!(alpha, aaM, X, Y)
            @test Y ≈ FMA2
            @test y ≈ fma2
        end
        @testset "getindex" begin
            N = 3
            M = 4
            sparseM = sprand(N, M, 0.5)
            aaM = AASparseMatrix(sparseM)
            @test size(aaM) == size(sparseM)
            for i in 1:N
                for j in 1:M
                    @test aaM[i,j] == sparseM[i,j]
                end
            end
            denseM = Array(sparseM)
            for i in 1:(N*M)
                @test denseM[i] == aaM[i]
            end
        end
        @testset "special matrices" begin
            N = 3
            sparseM = sprand(N, N, 0.5)
            while istriu(sparseM) || istril(sparseM)
                sparseM = sprand(N, N, 0.5)
            end
            ordinary = AASparseMatrix(sparseM)
            @test !issymmetric(ordinary) && !istril(ordinary) && !istriu(ordinary)
            symmetricM = sparse(sparseM + sparseM')
            @assert issymmetric(symmetricM) && !istril(symmetricM) && !istriu(symmetricM)
            symmetric = AASparseMatrix(symmetricM)
            @test issymmetric(symmetric) && !istril(symmetric) && !istriu(symmetric)
            upper_tri, lower_tri = sparse(I, N, N), sparse(I, N, N)
            upper_tri[1, 3] = 1.0
            lower_tri[3, 1] = 1.0
            @test !issymmetric(upper_tri) && !istril(upper_tri) && istriu(upper_tri)
            @test !issymmetric(lower_tri) && istril(lower_tri) && !istriu(lower_tri)
        end
    end

    @testset "AAFactorization" begin
        @testset "QR factor/solve/ldiv" begin
            for T in (Float32, Float64)
                @eval begin
                    N = 3
                    jlA = sprand($T, N, N, 0.5)
                    while det(jlA) ≈ 0
                        global jlA = sprand($T, N, N, 0.5)
                    end
                    test_fact = AAFactorization(jlA)
                    B = rand($T, N, N)
                    if issymmetric(jlA)
                        jlA[1,2] += 0.1
                    end
                    @test solve(test_fact, B) ≈ Array(jlA) \ B
                    b = rand($T, N)
                    @test solve(test_fact, b) ≈ Array(jlA) \ b
                    @test test_fact \ B ≈ Array(jlA) \ B
                    @test test_fact \ b ≈ Array(jlA) \ b
                end
            end
        end
        @testset "Error handling" begin
            N = 5
            err1 = nothing
            singular = spzeros(Float64, N, N)
            singular[1, 2] = 1.0
            try
                f = AAFactorization(singular)
                factor!(f)
            catch err1
            end
            @test sprint(showerror, err1) == "The matrix is singular.\n"

            err2 = nothing
            temp = sprand(N, N, 0.5)
            # add big negative value to lower-right corner of positive definite.
            nonPosDef = sparse(temp*temp' + diagm(rand(N)) +
                            diagm(vcat(zeros(Float64, N-1), [-1*float(N+1)])))
            try
                f = AAFactorization(nonPosDef)
                factor!(f, AppleAccelerate.SparseFactorizationCholesky)
            catch err2
            end
            @test err2 !== nothing

            err3 = nothing
            nonSymmetric = sparse(temp*temp' + diagm(rand(N)) + singular)
            try
                f = AAFactorization(nonSymmetric)
                @assert size(f.matrixObj)[1] == size(f.matrixObj)[2]
                factor!(f, AppleAccelerate.SparseFactorizationCholesky)
            catch err3
            end
            # they say "non-square:" I think they really mean "non-symmetric."
            @test sprint(showerror, err3) == "Cannot perform symmetric" *
                    " matrix factorization of non-square matrix.\n"
        end

        @testset "Cholesky" begin
            # create a random symmetric positive-definite matrix A.
            N = 4
            temp = sprand(N, N, 0.3)
            A = sparse(temp*temp' + diagm(rand(N)))
            sym_fact = AAFactorization(A)
            factor!(sym_fact)
            x, X = rand(N), rand(N, 4)
            b, B = A*x, A*X
            @test solve(sym_fact, b) ≈ x
            @test solve(sym_fact, B) ≈ X
        end

        @testset "non-square in-place solve" begin
            # using high density to avoid structurally singular matrices.
            tallMatrix = sprand(6,3,0.9)
            aa_fact = AAFactorization(tallMatrix)
            x, X = rand(3), rand(3, 3)
            b, B = tallMatrix * x, tallMatrix * X
            @test isapprox(solve!(aa_fact, b), x; 0.001)
            # solve!(aa_fact, B)
            # @test isapprox(B, X; 0.001)
            shortMatrix = sprand(3,4,0.9)
            aa_fact2 = AAFactorization(shortMatrix)
            x, X = rand(4), rand(4,4)
            b, B = shortMatrix * x, shortMatrix * X
            bx, BX = zeros(4), zeros(4,4)
            bx[1:3], BX[1:3, :] = b, B
            # Seems like julia and apple make different choices
            # when there's multiple solutions. which is odd because I thought
            # they both did the minimal norm solution...
            # @test isapprox(solve!(aa_fact2, copy(bx)), shortMatrix\b)
            @test isapprox(shortMatrix * solve!(aa_fact2, bx), b)

            # solve!(aa_fact2, BX)
            # @test BX ≈ X
        end
    end
    @testset "factorize" begin
        jlM = sprand(4,4, .9)
        aaM = AASparseMatrix(jlM)
        f = factorize(aaM)
        b = rand(4)
        @test f \ b ≈ jlM \ b
    end
end
