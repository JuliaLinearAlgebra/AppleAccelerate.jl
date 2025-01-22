@doc """Factorization object.

Create via `f = AAFactorization(A::SparseMatrixCSC{T, Int64})`. Calls to `solve`,
`ldiv`, and their in-place versions require explicitly passing in the
factorization object as the first argument. On construction, the struct stores a
placeholder yet-to-be-factored object: the factorization is computed upon the first call
to `solve`, or by explicitly calling `factor!`. If the matrix is symmetric, it defaults to
a Cholesky factorization; otherwise, it defaults to QR.
"""
mutable struct AAFactorization{T<:vTypes} <: LinearAlgebra.Factorization{T}
    matrixObj::AASparseMatrix{T}
    _factorization::SparseOpaqueFactorization{T}
end

# returns an AAFactorization containing A and a dummy "yet-to-be-factored" factorization.
function AAFactorization(A::AASparseMatrix{T}) where T<:vTypes
    obj = AAFactorization(A, SparseOpaqueFactorization(T))
    function cleanup(aa_fact)
        # If it's yet-to-be-factored, then there's nothing to release
        if !(aa_fact._factorization.status in (SparseYetToBeFactored,
                                                            SparseStatusReleased))
            SparseCleanup(aa_fact._factorization)
        end
    end
    return finalizer(cleanup, obj)
end

LinearAlgebra.factorize(A::AASparseMatrix{T}) where T<:vTypes = AAFactorization(A)

# julia's LinearAlgebra module doesn't provide similar constructors.
AAFactorization(M::SparseMatrixCSC{T, Int64}) where T<:vTypes =
                                        AAFactorization(AASparseMatrix(M))

# easiest way to make this follow the defaults and naming conventions of LinearAlgebra?
# TODO: add tests for the different kinds of factorizations, beyond QR.
function factor!(aa_fact::AAFactorization{T},
            kind::SparseFactorization_t = SparseFactorizationTBD) where T<:vTypes
    if aa_fact._factorization.status == SparseYetToBeFactored
        if kind == SparseFactorizationTBD
            # so far I'm only dealing with ordinary and symmetric
            kind = issymmetric(aa_fact.matrixObj) ? SparseFactorizationCholesky :
                            SparseFactorizationQR
        end
        aa_fact._factorization = SparseFactor(kind, aa_fact.matrixObj.matrix)
        if aa_fact._factorization.status == SparseMatrixIsSingular
            # throw a SingularException? Factoring a singular matrix usually does not
            # make it this far: the call to SparseFactor throws an error.
            throw(ErrorException("The matrix is singular."))
        elseif aa_fact._factorization.status == SparseStatusFailed
            throw(ErrorException("Factorization failed: check that the matrix"
                        * " has the correct properties for the factorization."))
        elseif aa_fact._factorization.status != SparseStatusOk
            throw(ErrorException("Something went wrong internally. Error type: "
                                * String(aa_fact._factorization.status)))
        end
    end
end

function solve(aa_fact::AAFactorization{T}, b::StridedVecOrMat{T}) where T<:vTypes
    @assert size(aa_fact.matrixObj)[2] == size(b, 1)
    factor!(aa_fact)
    x = Array{T}(undef, size(aa_fact.matrixObj)[2], size(b)[2:end]...)
    SparseSolve(aa_fact._factorization, b, x)
    return x
end

function solve!(aa_fact::AAFactorization{T}, xb::StridedVecOrMat{T}) where T<:vTypes
    @assert (xb isa StridedVector) ||
            (size(aa_fact.matrixObj)[1] == size(aa_fact.matrixObj)[2]) "Can't in-place " *
            "solve: x and b are different sizes and Julia cannot resize a matrix."
    factor!(aa_fact)
    SparseSolve(aa_fact._factorization, xb)
    return xb # written in imitation of KLU.jl, which also returns
end

LinearAlgebra.ldiv!(aa_fact::AAFactorization{T}, xb::StridedVecOrMat{T}) where T<:vTypes =
        solve!(aa_fact, xb)

function LinearAlgebra.ldiv!(x::StridedVecOrMat{T},
                            aa_fact::AAFactorization{T},
                            b::StridedVecOrMat{T}) where T<:vTypes
    @assert size(aa_fact.matrixObj)[2] == size(b, 1)
    factor!(aa_fact)
    SparseSolve(aa_fact._factorization, b, x)
    return x
end



