module AppleAccelerateLinearAlgebraSparseArraysExt

@static if Sys.isapple()

using AppleAccelerate
using AppleAccelerate: vTypes, AASparseMatrix,
                       SparseFactorization_t,
                       SparseFactorizationCholesky, SparseFactorizationQR,
                       SparseFactorizationLU, SparseFactorizationLDLT,
                       factor!, solve, refactor!
using LinearAlgebra
using SparseArrays
using SparseArrays: SparseMatrixCSC

# This extension provides the idiomatic LinearAlgebra entry points
# (`lu`/`cholesky`/`qr`/`ldlt`/`factorize`/`\`) and the `AAFactorization`
# constructor that accept a standard Julia `SparseMatrixCSC`. They require BOTH
# the `SparseMatrixCSC` type (SparseArrays) AND the `AAFactorization` machinery
# (defined in the LinearAlgebra extension), hence the dual trigger.
#
# `AppleAccelerate.AAFactorization` is installed into the AppleAccelerate module
# by the LinearAlgebra extension's `__init__`. Reference it through the parent
# module so this extension does not depend on the other extension's module path.

# AAFactorization construction from a CSC (the core lacks the type; the
# SparseArrays ext lacks the LinearAlgebra-defined AAFactorization).
function AppleAccelerate.AAFactorization(M::SparseMatrixCSC{T, Int64}) where T<:vTypes
    return AppleAccelerate.AAFactorization(AASparseMatrix(M))
end

# Build (but don't yet compute) an AAFactorization from a CSC, then force a
# specific factorization kind.
function _aa_factorize(A::SparseMatrixCSC{T,Int64}, kind::SparseFactorization_t) where {T<:vTypes}
    f = AppleAccelerate.AAFactorization(A)
    factor!(f, kind)
    return f
end

"""
    factorize(A::SparseMatrixCSC) -> AAFactorization

Return an `AAFactorization` of `A` using Apple's Sparse solvers, with the
factorization computed lazily on first solve (Cholesky for symmetric/Hermitian,
LU/QR otherwise — see [`AppleAccelerate.factor!`](@ref)).
"""
LinearAlgebra.factorize(A::SparseMatrixCSC{T,Int64}) where {T<:vTypes} =
    AppleAccelerate.AAFactorization(A)

"""
    cholesky(A::SparseMatrixCSC) -> AAFactorization

Cholesky factorization via Apple's Sparse solvers. `A` must be symmetric
(real) or Hermitian (complex) and positive definite.
"""
LinearAlgebra.cholesky(A::SparseMatrixCSC{T,Int64}) where {T<:vTypes} =
    _aa_factorize(A, SparseFactorizationCholesky)

"""
    ldlt(A::SparseMatrixCSC) -> AAFactorization

LDLᵀ factorization via Apple's Sparse solvers, for symmetric/Hermitian
(indefinite) matrices.
"""
LinearAlgebra.ldlt(A::SparseMatrixCSC{T,Int64}) where {T<:vTypes} =
    _aa_factorize(A, SparseFactorizationLDLT)

"""
    qr(A::SparseMatrixCSC) -> AAFactorization

QR factorization via Apple's Sparse solvers. Works for rectangular matrices and
solves least-squares systems through [`\\`](@ref).
"""
LinearAlgebra.qr(A::SparseMatrixCSC{T,Int64}) where {T<:vTypes} =
    _aa_factorize(A, SparseFactorizationQR)

"""
    lu(A::SparseMatrixCSC) -> AAFactorization

LU factorization via Apple's Sparse solvers (requires macOS 15.5+ for the
libSparse LU path).
"""
LinearAlgebra.lu(A::SparseMatrixCSC{T,Int64}) where {T<:vTypes} =
    _aa_factorize(A, SparseFactorizationLU)

"""
    A \\ b

Solve the sparse linear system `A x = b` (or least-squares for rectangular `A`)
using Apple's Sparse solvers, for `A::SparseMatrixCSC` and a dense `b`. Builds an
`AAFactorization` internally and discards it; for repeated solves with the same
`A`, build the factorization once with [`factorize`](@ref)/[`cholesky`](@ref)/…
and reuse it.
"""
Base.:\(A::SparseMatrixCSC{T,Int64}, b::StridedVecOrMat{T}) where {T<:vTypes} =
    solve(AppleAccelerate.AAFactorization(A), b)

# refactor! accepting a CSC (needs both AAFactorization and SparseMatrixCSC).
AppleAccelerate.refactor!(aa_fact, A::SparseMatrixCSC{T, Int64}) where T<:vTypes =
    refactor!(aa_fact, AASparseMatrix(A))

end # @static if Sys.isapple()

end # module
