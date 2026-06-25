module AppleAccelerateSparseArraysExt

@static if Sys.isapple()

using AppleAccelerate
using AppleAccelerate: AASparseMatrix, AAFactorization, vTypes, factor!, solve,
                       SparseFactorization_t,
                       SparseFactorizationCholesky, SparseFactorizationQR,
                       SparseFactorizationLU, SparseFactorizationLDLT,
                       ATT_TRANSPOSE, ATT_CONJUGATE_TRANSPOSE,
                       ATT_TRIANGLE_MASK, ATT_LOWER_TRIANGLE
using LinearAlgebra
using SparseArrays
using SparseArrays: SparseMatrixCSC

# ---------------------------------------------------------------------------
# Background / design note
# ---------------------------------------------------------------------------
# `SparseArrays` is a *hard* dependency of AppleAccelerate (the core
# `src/sparse.jl` builds `AASparseMatrix` directly from `SparseMatrixCSC` and
# uses `tril`/`diag`/`ishermitian` on it). Demoting it to a weakdep would break
# the core package, so we keep it as a regular dependency and layer this
# extension *on top* purely to add idiomatic LinearAlgebra-style entry points
# (`\`, `cholesky`, `lu`, `qr`, `ldlt`, `factorize`) that accept a standard
# Julia `SparseMatrixCSC`, plus a `SparseMatrixCSC <- AASparseMatrix`
# round-trip. Wiring this through the extension mechanism keeps the
# SparseArrays-typed glue cleanly separated and mirrors the layout of the other
# package extensions (e.g. AbstractFFTs).

# Apple's Sparse solvers index with Clong column pointers / Cint row indices and
# only support 64-bit-indexed CSC on the Julia side (`SparseMatrixCSC{T,Int64}`),
# matching the core `AASparseMatrix` constructor.
const _CSC = SparseMatrixCSC{<:vTypes, Int64}

# ---------------------------------------------------------------------------
# Conversion: AASparseMatrix -> SparseMatrixCSC
# ---------------------------------------------------------------------------

"""
    SparseMatrixCSC(A::AASparseMatrix)

Materialize an Accelerate `AASparseMatrix` back into a standard Julia
`SparseMatrixCSC`. The transpose/adjoint attribute bits are honored, so the
result equals the logical matrix `A` represents (`size`, `getindex`).
"""
function SparseArrays.SparseMatrixCSC(A::AASparseMatrix{T}) where {T<:vTypes}
    attrs = A.matrix.structure.attributes
    transposed = (attrs & ATT_TRANSPOSE) != 0
    conj_t     = (attrs & ATT_CONJUGATE_TRANSPOSE) != 0
    # Rebuild the raw *stored* CSC from the underlying buffers. These describe
    # the matrix in libSparse's own (untransposed, lower/upper-as-stored) layout.
    rawrows = Int(A.matrix.structure.rowCount)
    rawcols = Int(A.matrix.structure.columnCount)
    raw = SparseMatrixCSC{T,Int64}(rawrows, rawcols,
                                   Int64.(A._colptr) .+ 1,
                                   Int64.(A._rowval) .+ 1,
                                   copy(A._nzval))

    sym  = LinearAlgebra.issymmetric(A)
    herm = T <: Complex && LinearAlgebra.ishermitian(A)
    if sym || herm
        # Stored as one triangle (plus diagonal); reflect to the full matrix.
        # For symmetric/Hermitian kind, istril/istriu return false (they only
        # report the *triangular* kind), so read the triangle bit directly.
        lower = (attrs & ATT_TRIANGLE_MASK) == ATT_LOWER_TRIANGLE
        tri = lower ? tril(raw) : triu(raw)
        offdiag = lower ? tril(raw, -1) : triu(raw, 1)
        full = tri + (herm ? adjoint(offdiag) : transpose(offdiag))
        return SparseMatrixCSC{T,Int64}(full)
    end

    # Ordinary / triangular matrices: apply the view bits to the raw CSC.
    out = raw
    if transposed
        out = conj_t ? SparseMatrixCSC{T,Int64}(copy(adjoint(out))) :
                       SparseMatrixCSC{T,Int64}(copy(transpose(out)))
    end
    return out
end

# ---------------------------------------------------------------------------
# Factorization entry points accepting a SparseMatrixCSC
# ---------------------------------------------------------------------------

# Build (but don't yet compute) an AAFactorization from a CSC, then force a
# specific factorization kind.
function _aa_factorize(A::SparseMatrixCSC{T,Int64}, kind::SparseFactorization_t) where {T<:vTypes}
    f = AAFactorization(A)
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
    AAFactorization(A)

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

# ---------------------------------------------------------------------------
# Backslash: solve A x = b for a Julia sparse A
# ---------------------------------------------------------------------------

"""
    A \\ b

Solve the sparse linear system `A x = b` (or least-squares for rectangular `A`)
using Apple's Sparse solvers, for `A::SparseMatrixCSC` and a dense `b`. Builds an
`AAFactorization` internally and discards it; for repeated solves with the same
`A`, build the factorization once with [`factorize`](@ref)/[`cholesky`](@ref)/…
and reuse it.
"""
Base.:\(A::SparseMatrixCSC{T,Int64}, b::StridedVecOrMat{T}) where {T<:vTypes} =
    solve(AAFactorization(A), b)

end # @static if Sys.isapple()

end # module
