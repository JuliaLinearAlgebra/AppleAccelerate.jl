module AppleAccelerateSparseArraysExt

@static if Sys.isapple()

using AppleAccelerate
using AppleAccelerate: AASparseMatrix, vTypes,
                       att_type, ATT_ORDINARY, ATT_HERMITIAN, ATT_SYMMETRIC,
                       ATT_LOWER_TRIANGLE, ATT_TRI_LOWER, ATT_TRI_UPPER,
                       ATT_UNIT_TRIANGULAR, ATT_TRANSPOSE, ATT_CONJUGATE_TRANSPOSE,
                       ATT_TRIANGLE_MASK,
                       _is_symmetric_attr, _is_hermitian_attr
using SparseArrays
using SparseArrays: SparseMatrixCSC

# SparseArrays depends on (and loads) LinearAlgebra, but LinearAlgebra is only a
# *weak* dependency of AppleAccelerate, so this extension cannot `using
# LinearAlgebra` directly. We reach the (always-loaded) LinearAlgebra module
# through SparseArrays' own binding for the few predicates we need (ishermitian,
# istril, istriu, diag). tril/triu are re-exported by SparseArrays itself.
const LA = SparseArrays.LinearAlgebra

# Apple's Sparse solvers index with Clong column pointers / Cint row indices and
# only support 64-bit-indexed CSC on the Julia side (`SparseMatrixCSC{T,Int64}`).
const _CSC = SparseMatrixCSC{<:vTypes, Int64}

# ---------------------------------------------------------------------------
# Construction: SparseMatrixCSC -> AASparseMatrix
# ---------------------------------------------------------------------------

function AppleAccelerate.AASparseMatrix(sparseM::SparseMatrixCSC{T, Int64},
                        attributes::att_type = ATT_ORDINARY) where T<:vTypes
    if attributes == ATT_ORDINARY
        if LA.ishermitian(sparseM)
            kind = T <: Complex ? ATT_HERMITIAN : ATT_SYMMETRIC
            return AASparseMatrix(LA.tril(sparseM), kind | ATT_LOWER_TRIANGLE)
        elseif LA.istril(sparseM) || LA.istriu(sparseM)
            attributes = LA.istril(sparseM) ? ATT_TRI_LOWER : ATT_TRI_UPPER
        end
    end
    # NOTE: we deliberately do NOT auto-tag unit-triangular matrices. libSparse
    # treats ATT_UNIT_TRIANGULAR as an *implicit* unit diagonal, but we keep the
    # stored diagonal in the CSC data, so tagging it would double-count the
    # diagonal in both multiply and solve. The plain ATT_TRIANGULAR path is
    # correct for stored diagonals.
    c = Clong.(sparseM.colptr .+ -1)
    r = Cint.(sparseM.rowval .+ -1)
    vals = copy(sparseM.nzval)
    return AASparseMatrix(size(sparseM)..., c, r, vals, attributes)
end

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

    sym  = _is_symmetric_attr(A)
    herm = T <: Complex && _is_hermitian_attr(A)
    if sym || herm
        # Stored as one triangle (plus diagonal); reflect to the full matrix.
        # For symmetric/Hermitian kind, istril/istriu return false (they only
        # report the *triangular* kind), so read the triangle bit directly.
        lower = (attrs & ATT_TRIANGLE_MASK) == ATT_LOWER_TRIANGLE
        tri = lower ? LA.tril(raw) : LA.triu(raw)
        offdiag = lower ? LA.tril(raw, -1) : LA.triu(raw, 1)
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

end # @static if Sys.isapple()

end # module
