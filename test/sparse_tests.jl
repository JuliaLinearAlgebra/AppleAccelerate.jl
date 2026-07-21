using AppleAccelerate
using SparseArrays
using LinearAlgebra

import AppleAccelerate: AAFactorization, AASparseMatrix, factor!, muladd!, refactor!, solve, solve!

@testset "Sparse Linear Algebra" begin
    @testset "attribute bitfields" begin
        AA = AppleAccelerate
        # Independent C bitfields must not share bits with each other.
        @test AA.ATT_TRANSPOSE           & AA.ATT_KIND_MASK_COMPLEX == 0
        @test AA.ATT_LOWER_TRIANGLE      & AA.ATT_KIND_MASK_COMPLEX == 0
        @test AA.ATT_CONJUGATE_TRANSPOSE & AA.ATT_KIND_MASK_COMPLEX == 0
        @test AA.ATT_TRANSPOSE      & AA.ATT_LOWER_TRIANGLE      == 0
        @test AA.ATT_TRANSPOSE      & AA.ATT_CONJUGATE_TRANSPOSE == 0
        @test AA.ATT_LOWER_TRIANGLE & AA.ATT_CONJUGATE_TRANSPOSE == 0

        # OR-combined constants must let each field be extracted back via its mask.
        @test AA.ATT_TRI_LOWER & AA.ATT_KIND_MASK     == AA.ATT_TRIANGULAR
        @test AA.ATT_TRI_LOWER & AA.ATT_TRIANGLE_MASK == AA.ATT_LOWER_TRIANGLE
        @test AA.ATT_TRI_UPPER & AA.ATT_KIND_MASK     == AA.ATT_TRIANGULAR
        @test AA.ATT_TRI_UPPER & AA.ATT_TRIANGLE_MASK == AA.ATT_UPPER_TRIANGLE
    end

    @testset "wrappers" begin
        @testset "SparseMultiply and SparseMultiplyAdd" begin
            # less copy-paste heavy way via meta programming features?
            for T in (Float32, Float64)
                @eval begin
                    # Non-square (3 rows × 4 cols): the input operands have 4
                    # rows (the column count) and the results have 3.
                    dense = rand($T, 4,3)
                    denseV = rand($T, 4)
                    dense2 = zeros($T, 3,3)
                    denseV2 = zeros($T, 3)
                    sparseM = sprand($T, 3, 4, 0.5)
                    sparse_data = sparseM.nzval
                    col_inds =  Clong.(sparseM.colptr .+ -1)
                    row_inds = Cint.(sparseM.rowval .+ -1)
                    GC.@preserve col_inds row_inds sparse_data begin
                        s = AppleAccelerate.SparseMatrixStructure(3, 4,
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
            dense = rand(4,3)
            denseV = rand(4)
            dense2 = zeros(3, 3)
            denseV2 = zeros(3)
            sparseM = sprand(3, 4, 0.5)
            sparse_data = sparseM.nzval
            col_inds =  Clong.(sparseM.colptr .+ -1)
            row_inds = Cint.(sparseM.rowval .+ -1)
            # for lack of a way to call cconvert directly, I'll use the matrix multiply routines.
            GC.@preserve sparse_data col_inds row_inds begin
                s = AppleAccelerate.SparseMatrixStructure(3, 4,
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

        @testset "_libsparse_throw status mapping" begin
            # Code coverage reasons: no easy was to trigger the SparseInternalError.
            AA = AppleAccelerate
            @test AA._libsparse_throw(AA.SparseStatusOk, "x") === nothing
            @test_throws SingularException AA._libsparse_throw(
                AA.SparseMatrixIsSingular, "factor")
            @test_throws ArgumentError AA._libsparse_throw(
                AA.SparseParameterError, "factor")
            @test_throws ErrorException AA._libsparse_throw(
                AA.SparseStatusFailed, "factor")
            @test_throws ErrorException AA._libsparse_throw(
                AA.SparseInternalError, "factor")
            # Unknown / unmapped status falls through to the generic branch.
            @test_throws ErrorException AA._libsparse_throw(
                AA.SparseYetToBeFactored, "factor")
        end

        @testset "SparseFactorNoErrors" begin
            # SparseFactor(kind, M, true) should route to SparseFactorNoErrors,
            # which calls libSparse without our error-handling options wrapper.
            jlA = sprand(Float64, 5, 5, 0.5) + 5.0I
            aa = AASparseMatrix(jlA)
            qr = AppleAccelerate.SparseFactorizationQR
            f = AppleAccelerate.SparseFactor(qr, aa.matrix, true)
            @test f.status == AppleAccelerate.SparseStatusOk
            # Use it to solve so we know it produced a working factorization.
            x = rand(5); b = jlA * x
            xb = copy(b)
            AppleAccelerate.SparseSolve(f, xb)
            @test xb ≈ x
            AppleAccelerate.SparseCleanup(f)
            # TBD is rejected at the SparseFactorNoErrors entry point.
            @test_throws ArgumentError AppleAccelerate.SparseFactorNoErrors(
                AppleAccelerate.SparseFactorizationTBD, aa.matrix)
        end
    end
    @testset "AASparseMatrix" begin
        @testset "arithmetic" begin
            # Non-square (R rows, C cols) so a row/column mix-up in the
            # multiply/muladd wrappers can't slip through.
            R, C = 10, 7
            sparseM = sprand(R, C, 0.3)
            alpha = rand(Float64)
            aaM = AASparseMatrix(sparseM)
            X, x = rand(C, 4), rand(C)
            @test aaM*x ≈ sparseM*x
            @test aaM*X ≈ sparseM*X
            @test alpha*aaM*x ≈ alpha*sparseM*x
            @test alpha*aaM*X ≈ alpha*sparseM*X
            Y, y = rand(R, 4), rand(R)
            FMA, fma = sparseM*X + Y, sparseM*x + y
            muladd!(aaM, x, y)
            muladd!(aaM, X, Y)
            @test Y ≈ FMA
            @test y ≈ fma
            Y, y = rand(R, 4), rand(R)
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
        @testset "transpose" begin
            N = 3
            M = 4
            sparseM = sprand(N, M, 0.5)
            aaM = AASparseMatrix(sparseM)
            aaMt = transpose(aaM)
            # SparseGetTranspose sets the transpose flag; the underlying structure
            # retains the original dimensions, so verify the flag is toggled.
            ATT_TRANSP = AppleAccelerate.ATT_TRANSPOSE
            @test (aaMt.matrix.structure.attributes & ATT_TRANSP) != zero(typeof(ATT_TRANSP))
            # Double transpose clears the flag.
            aaMtt = transpose(aaMt)
            @test (aaMtt.matrix.structure.attributes & ATT_TRANSP) == zero(typeof(ATT_TRANSP))
            # Dimensions reverse under transpose.
            @test size(aaM)  == (N, M)
            @test size(aaMt) == (M, N)
            # Adjoint is tranpose.
            aaMa = adjoint(aaM)
            @test (aaMa.matrix.structure.attributes & ATT_TRANSP) != zero(typeof(ATT_TRANSP))
            @test size(aaMa) == (M, N)
            x = rand(N)
            @test aaMa * x ≈ transpose(Array(sparseM)) * x

            # transpose(M)[i,j] equals M[j, i]
            for I in 1:M, J in 1:N
                @test aaMt[I, J] == sparseM[J, I]
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
            # verify that the julia matrix is has the desired properties.
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
                    N = 100
                    # Generate a well-conditioned sparse matrix by adding a
                    # diagonal shift, avoiding flaky failures from
                    # ill-conditioned random matrices.
                    jlA = sprand($T, N, N, 0.1) + $T(N) * I
                    test_fact = AAFactorization(jlA)
                    B = rand($T, N, N)
                    @test solve(test_fact, B) ≈ Array(jlA) \ B
                    b = rand($T, N)
                    @test solve(test_fact, b) ≈ Array(jlA) \ b
                    @test test_fact \ B ≈ Array(jlA) \ B
                    @test test_fact \ b ≈ Array(jlA) \ b

                    # Overdetermined (tall) system: QR solve should match the
                    # dense least-squares solution. The rhs has N rows and the
                    # solution has M = N ÷ 2 rows.
                    M = N ÷ 2
                    tallA = sprand($T, N, M, 0.3) + $T(N) * [I; spzeros($T, N - M, M)]
                    tall_fact = AAFactorization(tallA)
                    Btall = rand($T, N, 3)
                    @test solve(tall_fact, Btall) ≈ Array(tallA) \ Btall
                    btall = rand($T, N)
                    @test solve(tall_fact, btall) ≈ Array(tallA) \ btall
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
            @test err1 isa SingularException ||
                  occursin("singular", lowercase(sprint(showerror, err1)))

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
            # Cholesky of a non-positive-definite matrix must raise, and the
            # error should actually describe a factorization/property failure
            # rather than being any arbitrary exception.
            @test err2 isa Exception
            @test occursin(r"properties|singular|factoriz|parameter|positive",
                           lowercase(sprint(showerror, err2)))

            err3 = nothing
            nonSymmetric = sparse(temp*temp' + diagm(rand(N)) + singular)
            try
                f = AAFactorization(nonSymmetric)
                # Verify the matrix is square
                @test size(f.matrixObj)[1] == size(f.matrixObj)[2]
                factor!(f, AppleAccelerate.SparseFactorizationCholesky)
            catch err3
            end
            @test startswith(sprint(showerror, err3),
                "Cannot perform symmetric matrix factorization"
            )
        end

        @testset "DimensionMismatch errors" begin
            N = 4
            M = 5
            jlA = sprand(Float64, N, M, 0.5)
            while rank(Array(jlA)) < min(N, M)
                jlA = sprand(Float64, N, M, 0.5)
            end
            f = AAFactorization(jlA)
            # For an N×M matrix the right-hand side must have N rows (one entry
            # per equation) and the solution must have M rows (one per unknown).

            # Test solve with wrong dimension vector (needs N rows, not M)
            wrong_b = rand(Float64, M)
            @test_throws DimensionMismatch solve(f, wrong_b)

            # Test solve with wrong dimension matrix (needs N rows, not M)
            wrong_B = rand(Float64, M, 3)
            @test_throws DimensionMismatch solve(f, wrong_B)

            # Test ldiv! with a correctly-sized rhs but wrong-sized output
            good_b = rand(Float64, N)
            wrong_x = rand(Float64, N)  # should have M rows
            @test_throws DimensionMismatch LinearAlgebra.ldiv!(wrong_x, f, good_b)
        end

        @testset "ldiv!(x, f, b) writes solution" begin
            # Exercises the successful body of the 3-arg ldiv!(x, f, b), checked
            # against the dense solve. Vector and matrix RHS, square system.
            N = 6
            A = sprandn(N, N, 0.5) + N * I
            f = AAFactorization(A)
            x0 = rand(N); b = A * x0
            x = similar(x0)
            ret = LinearAlgebra.ldiv!(x, f, b)
            @test ret === x
            @test isapprox(x, Matrix(A) \ b; rtol = 1e-8)

            X0 = rand(N, 3); B = A * X0
            X = similar(X0)
            @test isapprox(LinearAlgebra.ldiv!(X, f, B), Matrix(A) \ B; rtol = 1e-8)
        end

        @testset "ArgumentError for in-place solve" begin
            N = 4
            M = 5
            # Test that in-place solve with matrix throws ArgumentError for non-square factorization
            jlA = sprand(Float64, N, M, 0.9)
            while rank(Array(jlA)) < min(N, M)
                jlA = sprand(Float64, N, M, 0.9)
            end
            f = AAFactorization(jlA)

            xb_matrix = rand(Float64, M, 3)
            @test_throws ArgumentError solve!(f, xb_matrix)
        end

        @testset "lu! on unfactored F throws" begin
            A = sprandn(6, 6, 0.4) + 10I
            # refactor!-style spellings require F to already hold a factorization.
            @test_throws ArgumentError lu!(AAFactorization(A), A)
        end

        @testset "solve!(f, B::Matrix) square" begin
            N = 5
            A = sprandn(N, N, 0.5) + N * I
            f = AAFactorization(A)
            X = rand(N, 3)
            B = A * X
            solve!(f, B)             # square system: B overwritten with solution
            @test isapprox(B, X; rtol = 1e-8)
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

        @testset "LDLT" begin
            N = 4
            # symmetric indefinite matrix (not positive definite)
            temp = sprand(N, N, 0.9)
            A = sparse(temp + temp')
            while det(A) ≈ 0 || !issymmetric(A)
                temp = sprand(N, N, 0.9)
                A = sparse(temp + temp')
            end
            ldlt_fact = AAFactorization(A)
            factor!(ldlt_fact, AppleAccelerate.SparseFactorizationLDLTTPP)
            x = rand(N)
            b = A * x
            @test solve(ldlt_fact, b) ≈ x
        end

        something(AppleAccelerate.get_macos_version(), v"0.0.0") >= v"15.5" && @testset "LU" begin
            for T in (Float32, Float64)
                @eval begin
                    N = 50
                    # Square, non-symmetric, well-conditioned.
                    jlA = sprand($T, N, N, 0.1) + $T(N) * I
                    @test !issymmetric(jlA)
                    lu_fact = AAFactorization(jlA)
                    # Default for square non-symmetric should be LU.
                    factor!(lu_fact)
                    # Default LU currently resolves to LUTPP internally.
                    @test lu_fact._factorization.symbolicFactorization.type in
                        (AppleAccelerate.SparseFactorizationLU,
                         AppleAccelerate.SparseFactorizationLUUnpivoted,
                         AppleAccelerate.SparseFactorizationLUSPP,
                         AppleAccelerate.SparseFactorizationLUTPP)
                    x = rand($T, N)
                    X = rand($T, N, 3)
                    b, B = jlA * x, jlA * X
                    @test solve(lu_fact, b) ≈ x
                    @test solve(lu_fact, B) ≈ X

                    # Also exercise the explicit LU variants.
                    for kind in (AppleAccelerate.SparseFactorizationLUUnpivoted,
                                 AppleAccelerate.SparseFactorizationLUSPP,
                                 AppleAccelerate.SparseFactorizationLUTPP)
                        f2 = AAFactorization(jlA)
                        factor!(f2, kind)
                        @test solve(f2, b) ≈ x
                    end
                end
            end
        end

        @testset "non-square in-place solve" begin
            # using high density to avoid structurally singular matrices.
            tallMatrix = sprand(6,3,0.9)
            aa_fact = AAFactorization(tallMatrix)
            x, X = rand(3), rand(3, 3)
            b, B = tallMatrix * x, tallMatrix * X
            @test isapprox(solve!(aa_fact, b), x; atol=1e-3)
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
    # Complex sparse support is macOS 15.5+ only.
    something(AppleAccelerate.get_macos_version(), v"0.0.0") >= v"15.5" && @testset "Complex" begin
        @testset "multiply $T" for T in (ComplexF32, ComplexF64)
            N = 30
            # Diagonal shift keeps the matrix well-conditioned for the factor/solve tests.
            jlA = sprand(T, N, N, 0.1) + T(N) * I
            aaA = AASparseMatrix(jlA)
            x = rand(T, N)
            X = rand(T, N, 3)
            # Sparse * dense matches the dense product.
            @test aaA * x ≈ Array(jlA) * x
            @test aaA * X ≈ Array(jlA) * X
            # Scaled multiply (exercises the Cf/Cd-mangled scalar variant).
            α = T(2) + T(im) * T(3)
            @test α * aaA * x ≈ Array(jlA) * (α * x)  # α*A*x = A*(α*x)
            # muladd!: y += A * x.
            y = rand(T, N); y0 = copy(y)
            muladd!(aaA, x, y)
            @test y ≈ y0 + Array(jlA) * x
        end

        @testset "adjoint vs transpose $T" for T in (ComplexF32, ComplexF64)
            N = 8
            jlA = sprand(T, N, N, 0.5) + T(N) * I
            aaA = AASparseMatrix(jlA)
            x = rand(T, N)
            # transpose and adjoint are distinct for complex.
            @test transpose(aaA) * x ≈ transpose(Array(jlA)) * x
            @test adjoint(aaA)  * x ≈ adjoint(Array(jlA))  * x
            @test transpose(aaA) * x ≉ adjoint(aaA) * x
            # Round-trip: (A')' has the same action as A.
            @test adjoint(adjoint(aaA)) * x ≈ aaA * x
            # getindex on adjoint should swap indices AND conjugate.
            adjA = adjoint(aaA)
            for I in 1:N, J in 1:N
                @test adjA[I, J] == conj(jlA[J, I])
            end
            # adjoint(adjoint(M)) lands in the "no meaning" (F, T) attribute state.
            # libSparse treats that as identity, and so must getindex — i.e., no
            # accidental conjugation here. Regression guard for the `ct alone`
            # check that would have silently broken this.
            adjAdj = adjoint(adjoint(aaA))
            for I in 1:N, J in 1:N
                @test adjAdj[I, J] == jlA[I, J]
            end
            # libSparse traps on transpose-of-adjoint or adjoint-of-transpose
            # Our wrappers convert those to clean Julia errors.
            @test_throws ArgumentError transpose(adjoint(aaA))
            @test_throws ArgumentError adjoint(transpose(aaA))
            # (F,T) state from adjoint(adjoint) is libSparse-identity, so
            # transposing/re-adjointing it is legal and acts on the original.
            @test transpose(adjoint(adjoint(aaA))) * x ≈ transpose(Array(jlA)) * x
            @test adjoint(adjoint(adjoint(aaA)))   * x ≈ adjoint(Array(jlA))   * x
        end

        @testset "complex symmetric (not Hermitian)" begin
            # Older macOS libSparse rejects LDLT for complex-symmetric matrices
            # ("Cannot perform Hermitian matrix factorization..."), even though
            # newer macOS implements it as LDLᵀ correctly. For cross-version
            # portability we don't auto-tag — these get stored as ordinary and
            # factored via LU on full storage. See `factor!` docs for details.
            N = 8
            T = ComplexF64
            M = sprand(T, N, N, 0.3)
            A = sparse(M + transpose(M)) + T(N) * I
            @assert issymmetric(A) && !ishermitian(A)
            aaA = AASparseMatrix(A)
            @test !ishermitian(aaA) && !issymmetric(aaA)  # stored as ordinary
            f = AAFactorization(aaA)
            x = rand(T, N)
            @test solve(f, A * x) ≈ x
        end

        @testset "Hermitian detection and factor!" begin
            N = 10
            T = ComplexF64
            # Build a Hermitian positive-definite matrix.
            temp = sprand(T, N, N, 0.3)
            H = sparse(temp * adjoint(temp) + N * I)
            @test ishermitian(H) && !issymmetric(H)
            aaH = AASparseMatrix(H)
            @test ishermitian(aaH)
            @test !issymmetric(aaH)
            # factor! default for a Hermitian matrix should pick Cholesky.
            f = AAFactorization(aaH)
            factor!(f)
            @test f._factorization.symbolicFactorization.type ==
                  AppleAccelerate.SparseFactorizationCholesky
            x = rand(T, N)
            @test solve(f, H * x) ≈ x
        end

        @testset "factor/solve square non-Hermitian $T" for T in (ComplexF32, ComplexF64)
            N = 20
            jlA = sprand(T, N, N, 0.2) + T(N) * I
            f = AAFactorization(jlA)
            # Default for square non-Hermitian complex is LU.
            x = rand(T, N)
            B = rand(T, N, 2)
            @test solve(f, jlA * x) ≈ x
            @test solve(f, jlA * B) ≈ B
            @test f._factorization.symbolicFactorization.type in
                (AppleAccelerate.SparseFactorizationLU,
                 AppleAccelerate.SparseFactorizationLUUnpivoted,
                 AppleAccelerate.SparseFactorizationLUSPP,
                 AppleAccelerate.SparseFactorizationLUTPP)
        end
    end

    @testset "factorize" begin
        # Use 100x100 to avoid singular matrices from small random sparse matrices
        jlM = sprandn(100, 100, 0.3) + 10I
        aaM = AASparseMatrix(jlM)
        f = factorize(aaM)
        b = rand(100)
        @test f \ b ≈ jlM \ b
    end

    # Regression: unit-triangular auto-tagging used to OR in ATT_UNIT_TRIANGULAR
    # while leaving the diagonal stored, so libSparse double-counted the diagonal
    # (implicit unit diagonal + stored diagonal) in both multiply and solve.
    @testset "unit-diagonal triangular not double-counted" begin
        for T in (Float32, Float64)
            # Lower-triangular with an all-ones diagonal.
            A = sparse(T[1 0 0; 5 1 0; 2 3 1])
            aa = AASparseMatrix(A)
            x = T[2, 3, 4]
            @test aa * x ≈ A * x
            # Specifically the case from the bug report.
            A2 = sparse(T[1 0; 5 1])
            aa2 = AASparseMatrix(A2)
            @test aa2 * T[2, 3] ≈ A2 * T[2, 3] ≈ T[2, 13]
            # Solve must also be correct (not corrupted by an implicit diagonal).
            f = AAFactorization(aa)
            @test solve(f, A * x) ≈ x
        end
    end

    # Regression: getindex used to read raw CSC and ignore the symmetric/Hermitian
    # attribute, returning 0 for the un-stored triangle.
    @testset "getindex mirrors symmetric/Hermitian un-stored triangle" begin
        # Real symmetric.
        S = sparse([1.0 2.0 0.0; 2.0 3.0 4.0; 0.0 4.0 5.0])
        @test issymmetric(S)
        aaS = AASparseMatrix(S)
        @test Array(aaS) == Array(S)
        @test aaS[1, 2] == 2.0
        @test aaS[2, 3] == 4.0

        # Complex Hermitian: upper triangle is conj of lower.
        H = sparse(ComplexF64[2 (1+2im) 0; (1-2im) 3 (0-1im); 0 (0+1im) 4])
        @test ishermitian(H)
        aaH = AASparseMatrix(H)
        @test Array(aaH) == Array(H)
        @test aaH[1, 2] == conj(H[2, 1]) == (1 + 2im)
    end

    @testset "refactor!" begin
        @testset "real $T" for T in (Float32, Float64)
            # Diagonally dominant so QR/LU refactor is well-conditioned.
            A1 = sprandn(T, 40, 40, 0.2) + T(10) * I
            f = AAFactorization(AASparseMatrix(A1))
            factor!(f)  # compute symbolic + numeric once
            # Same pattern, different values (scaling preserves the stored pattern).
            A2 = copy(A1)
            A2.nzval .*= T(1.5)
            b = rand(T, 40)
            refactor!(f, AASparseMatrix(A2))
            xref = solve(f, b)
            # Must match a fresh factorization of A2.
            xfresh = solve(AAFactorization(AASparseMatrix(A2)), b)
            @test isapprox(xref, xfresh; rtol = sqrt(eps(T)))
            @test isapprox(Matrix(A2) * xref, b; rtol = 10 * sqrt(eps(T)))
        end

        @testset "via SparseMatrixCSC" begin
            A1 = sprandn(30, 30, 0.2) + 10I
            f = AAFactorization(AASparseMatrix(A1))
            solve(f, rand(30))  # forces factorization
            A2 = copy(A1); A2.nzval .*= 1.25
            b = rand(30)
            refactor!(f, A2)  # CSC overload
            @test isapprox(solve(f, b), Matrix(A2) \ b; rtol = 1e-8)
        end

        @testset "errors before factorization" begin
            A = sprandn(5, 5, 0.4) + 10I
            f = AAFactorization(AASparseMatrix(A))
            # not yet factored
            @test_throws ArgumentError refactor!(f, AASparseMatrix(A))
            factor!(f)
            # dimension mismatch
            @test_throws DimensionMismatch refactor!(f, AASparseMatrix(sprandn(6, 6, 0.4) + 10I))
            # different nnz (sparsity pattern)
            B = copy(A); B[1, 1] == 0 && (B[1, 1] = 1.0)
            different = copy(A)
            # add a stored zero somewhere not already present
            idx = findfirst(==(0.0), Matrix(A))
            if idx !== nothing
                different[idx] = 1.0
                @test_throws ArgumentError refactor!(f, AASparseMatrix(different))
            end
        end

        @testset "lu!/cholesky!/ldlt! spellings" begin
            b = rand(20)
            # LU: nonsymmetric square -> factor! defaults to LU.
            A1 = sprandn(20, 20, 0.3) + 10I
            f = lu(AASparseMatrix(A1))
            A2 = copy(A1); A2.nzval .*= 1.3
            lu!(f, AASparseMatrix(A2))
            @test isapprox(solve(f, b), Matrix(A2) \ b; rtol = 1e-8)
            @test_throws ArgumentError cholesky!(f, AASparseMatrix(A2))  # wrong kind

            # Cholesky / LDLᵀ on a symmetric positive-definite matrix.
            M = sprandn(20, 20, 0.3); S = SparseMatrixCSC(M * M' + 20I)
            S2 = copy(S); S2.nzval .*= 1.1
            fc = cholesky(AASparseMatrix(S))
            cholesky!(fc, S2)                       # CSC overload
            @test isapprox(solve(fc, b), Matrix(S2) \ b; rtol = 1e-7)
            @test_throws ArgumentError lu!(fc, S2)  # wrong kind

            fl = AAFactorization(AASparseMatrix(S))
            factor!(fl, AppleAccelerate.SparseFactorizationLDLT)
            ldlt!(fl, AASparseMatrix(S2))
            @test isapprox(solve(fl, b), Matrix(S2) \ b; rtol = 1e-7)
        end

        # Complex sparse support (and the Hermitian refactor kernels) require
        # macOS 15.5+. Mirrors the real refactor coverage above. Regression for
        # the Hermitian-refactor crash: libSparse's *Symmetric_Complex* refactor
        # rejects Hermitian factorizations, so Cholesky/LDLT refactor of a complex
        # Hermitian matrix must route through the dedicated Hermitian kernels.
        something(AppleAccelerate.get_macos_version(), v"0.0.0") >= v"15.5" &&
        @testset "complex $T" for T in (ComplexF64, ComplexF32)
            rtol = T == ComplexF32 ? 1e-2 : 1e-7
            N = 12

            # Build a Hermitian positive-definite matrix and a same-pattern variant.
            herm(B) = sparse(Hermitian(B + B' + N * I))
            B1 = sprandn(T, N, N, 0.3)
            H1 = herm(B1)
            @test ishermitian(H1) && !issymmetric(H1)
            # Same sparsity pattern, scaled values (scaling preserves Hermitian-ness
            # and the stored pattern).
            H2 = copy(H1); H2.nzval .*= T(1.3)
            @test ishermitian(H2)
            b = rand(T, N)

            @testset "Hermitian Cholesky refactor" begin
                f = cholesky(AASparseMatrix(H1))
                cholesky!(f, AASparseMatrix(H2))
                xref = solve(f, b)
                xfresh = solve(cholesky(AASparseMatrix(H2)), b)
                @test isapprox(xref, xfresh; rtol = rtol)
                @test isapprox(Matrix(H2) * xref, b; rtol = 10 * rtol)
            end

            @testset "complex LU refactor" begin
                # Square, non-Hermitian, diagonally dominant.
                A1 = sprandn(T, N, N, 0.3) + T(N) * I
                @test !ishermitian(A1)
                A2 = copy(A1); A2.nzval .*= T(1.4)
                f = lu(AASparseMatrix(A1))
                refactor!(f, AASparseMatrix(A2))
                xref = solve(f, b)
                xfresh = solve(lu(AASparseMatrix(A2)), b)
                @test isapprox(xref, xfresh; rtol = rtol)
                @test isapprox(Matrix(A2) * xref, b; rtol = 10 * rtol)
            end

            @testset "complex QR refactor" begin
                A1 = sprandn(T, N, N, 0.3) + T(N) * I
                A2 = copy(A1); A2.nzval .*= T(1.2)
                f = qr(AASparseMatrix(A1))
                refactor!(f, AASparseMatrix(A2))
                xref = solve(f, b)
                xfresh = solve(qr(AASparseMatrix(A2)), b)
                @test isapprox(xref, xfresh; rtol = rtol)
                @test isapprox(Matrix(A2) * xref, b; rtol = 10 * rtol)
            end
        end

        # Hermitian LDLT refactor must route through the dedicated Hermitian
        # kernel (the Symmetric_Complex kernel rejects Hermitian factorizations).
        # This is verified in a *fresh subprocess*: Apple's libSparse has a
        # process-global defect where any prior complex-Hermitian factorization
        # (e.g. the Hermitian Cholesky factorize exercised earlier in this suite)
        # poisons a later complex-Hermitian LDLT refactor, crashing inside
        # libSparse's own `_SparseGetWorkspaceRequired`. Isolating it in its own
        # process avoids that libSparse bug while still covering our routing fix.
        something(AppleAccelerate.get_macos_version(), v"0.0.0") >= v"15.5" &&
        @testset "Hermitian LDLT refactor (isolated subprocess)" begin
            code = """
            using AppleAccelerate, LinearAlgebra, SparseArrays
            import AppleAccelerate: AAFactorization, AASparseMatrix, factor!, solve
            T = ComplexF64; N = 12
            B1 = sprandn(T, N, N, 0.3)
            H1 = sparse(Hermitian(B1 + B1' + N*I))
            H2 = copy(H1); H2.nzval .*= T(1.3)
            b = rand(T, N)
            f = AAFactorization(AASparseMatrix(H1))
            factor!(f, AppleAccelerate.SparseFactorizationLDLT)
            ldlt!(f, AASparseMatrix(H2))                 # Hermitian LDLT refactor
            xref = solve(f, b)
            ok = isapprox(Matrix(H2) * xref, b; rtol = 1e-6)
            exit(ok ? 0 : 1)
            """
            proj = Base.active_project()
            p = run(ignorestatus(`$(Base.julia_cmd()) --project=$(proj) -e $code`))
            @test p.exitcode == 0
        end
    end

    @testset "SparseMatrixCSC(::AASparseMatrix) round-trip" begin
        @testset "ordinary $T" for T in (Float64, ComplexF64)
            M = T <: Complex ? sprandn(T, 7, 7, 0.4) : sprandn(7, 7, 0.4)
            aa = AASparseMatrix(M)
            @test SparseMatrixCSC(aa) == M
        end

        @testset "transpose / adjoint" begin
            M = sprandn(ComplexF64, 6, 6, 0.4)
            aa = AASparseMatrix(M)
            @test SparseMatrixCSC(transpose(aa)) == sparse(transpose(M))
            @test SparseMatrixCSC(adjoint(aa)) == sparse(adjoint(M))
        end

        @testset "symmetric" begin
            S = sparse(Symmetric(sprandn(8, 8, 0.3)))
            @test issymmetric(S)
            @test SparseMatrixCSC(AASparseMatrix(S)) == S
        end

        @testset "Hermitian" begin
            B = sprandn(ComplexF64, 8, 8, 0.3)
            H = sparse(Hermitian(B + 8I))
            @test ishermitian(H)
            @test SparseMatrixCSC(AASparseMatrix(H)) == H
        end

        @testset "triangular $T" for T in (Float64, ComplexF64)
            B = T <: Complex ? sprandn(T, 6, 6, 0.5) : sprandn(6, 6, 0.5)
            L = sparse(tril(B))
            U = sparse(triu(B))
            @test SparseMatrixCSC(AASparseMatrix(L)) == L
            @test SparseMatrixCSC(AASparseMatrix(U)) == U
        end
    end

    @testset "factorization entry points on AASparseMatrix" begin
        @testset "qr (real $T)" for T in (Float32, Float64)
            M = sprandn(T, 20, 20, 0.3) + T(5) * I
            b = rand(T, 20)
            f = qr(AASparseMatrix(M))
            @test f isa AAFactorization
            @test isapprox(solve(f, b), Matrix(M) \ b; rtol = 10 * sqrt(eps(T)))
            # generic \ via Factorization supertype also works
            @test isapprox(f \ b, Matrix(M) \ b; rtol = 10 * sqrt(eps(T)))
        end

        @testset "cholesky (SPD)" begin
            B = sprandn(15, 15, 0.3)
            A = sparse(B * B' + 15I)  # symmetric positive definite
            b = rand(15)
            f = cholesky(AASparseMatrix(A))
            @test isapprox(solve(f, b), Matrix(A) \ b; rtol = 1e-7)
        end

        @testset "ldlt (symmetric)" begin
            B = sprandn(12, 12, 0.3)
            A = sparse(B + B' + 12I)  # symmetric
            b = rand(12)
            f = ldlt(AASparseMatrix(A))
            @test isapprox(solve(f, b), Matrix(A) \ b; rtol = 1e-7)
        end

        if something(AppleAccelerate.get_macos_version(), v"0.0.0") >= v"15.5"
            @testset "lu (square, macOS 15.5+)" begin
                M = sprandn(18, 18, 0.3) + 10I
                b = rand(18)
                f = lu(AASparseMatrix(M))
                @test isapprox(solve(f, b), Matrix(M) \ b; rtol = 1e-7)
            end
        end

        @testset "does not pirate Julia's sparse lu" begin
            # Loading AppleAccelerate must not override Julia's UMFPACK lu for CSC.
            m = which(LinearAlgebra.lu, Tuple{SparseMatrixCSC{Float64,Int}})
            @test parentmodule(m) !== AppleAccelerate
        end
    end

    # Regression for the workspace-taking SparseConvertFromCoord dispatch: the
    # Float32 method must call the `float` C symbol and Float64 the `double` one.
    # A swap (the bug that was fixed) would make the Float32 path read 8 bytes per
    # element from a Float32 array — garbage values and an out-of-bounds read — so
    # a functional coordinate→CSC round-trip catches it for both element types.
    @testset "SparseConvertFromCoord workspace dispatch::$T" for T in (Float32, Float64)
        AA = AppleAccelerate
        Mref = T[10 0 30; 0 20 0; 40 0 50]
        rows, cols = size(Mref)
        # Coordinate (COO) triplets, 0-based indices as Accelerate expects.
        rowidx = Cint[0, 2, 1, 0, 2]
        colidx = Cint[0, 0, 1, 2, 2]
        vals   = T[10, 40, 20, 30, 50]
        nnz = length(vals)
        blockSize = 1

        # Caller-owned storage/workspace, sized per Accelerate's Solve.h contract.
        storage_bytes = 48 + (cols + 1) * sizeof(Clong) + nnz * sizeof(Cint) +
                        nnz * blockSize * blockSize * sizeof(T)
        storage   = Vector{UInt8}(undef, storage_bytes)
        workspace = Vector{UInt8}(undef, rows * sizeof(Cint))

        dense = GC.@preserve rowidx colidx vals storage workspace begin
            sm = AA.SparseConvertFromCoord(
                Cint(rows), Cint(cols), Clong(nnz), Cuchar(blockSize),
                AA.ATT_ORDINARY,
                pointer(rowidx), pointer(colidx), pointer(vals),
                Ptr{Cvoid}(pointer(storage)), Ptr{Cvoid}(pointer(workspace)))
            # The returned CSC arrays point into `storage`; keep it alive while reading.
            colStarts = unsafe_wrap(Array, sm.structure.columnStarts, cols + 1)
            rowInds   = unsafe_wrap(Array, sm.structure.rowIndices, nnz)
            dataOut   = unsafe_wrap(Array, sm.data, nnz)
            D = zeros(T, rows, cols)
            for c in 1:cols
                for k in (colStarts[c] + 1):colStarts[c + 1]   # colStarts are 0-based
                    D[rowInds[k] + 1, c] += dataOut[k]
                end
            end
            D
        end
        @test dense ≈ Mref
    end

    # COO -> CSC constructor built on the workspace SparseConvertFromCoord.
    @testset "AASparseMatrix from COO triplets::$T" for T in (Float32, Float64)
        # Random reference via SparseArrays.sparse, then rebuild from its triplets.
        ref = sprandn(T, 12, 9, 0.3)
        I, J, Vv = findnz(ref)
        A = AASparseMatrix(I, J, Vv, 12, 9)
        @test Matrix(A) ≈ Matrix(ref)
        # Functional check: matches a dense multiply.
        x = randn(T, 9)
        @test A * x ≈ Matrix(ref) * x rtol = sqrt(eps(T))

        # Duplicate coordinates are summed (like SparseArrays.sparse).
        Ad = AASparseMatrix(Int[1, 1, 2], Int[1, 1, 2], T[2, 3, 5], 2, 2)
        @test Matrix(Ad) ≈ T[5 0; 0 5]

        # Input validation.
        @test_throws DimensionMismatch AASparseMatrix(Int[1, 2], Int[1], T[1, 2], 2, 2)
        @test_throws ArgumentError AASparseMatrix(Int[1, 3], Int[1, 1], T[1, 2], 2, 2)  # row OOR
        @test_throws ArgumentError AASparseMatrix(Int[1, 1], Int[1, 3], T[1, 2], 2, 2)  # col OOR
    end
end
