
@doc """Matrix wrapper, containing the Apple sparse matrix struct
and the pointed-to data. Construct from a `SparseMatrixCSC`.

Multiplication (`*`) and multiply-add (`muladd!`) with both
`Vector` and `Matrix` objects are working.
`transpose` creates a new matrix structure with the opposite
transpose flag, that references the same CSC data.
"""
mutable struct AASparseMatrix{T<:vTypes} <: AbstractMatrix{T}
    matrix::SparseMatrix{T}
    _colptr::Vector{Clong}
    _rowval::Vector{Cint}
    _nzval::Vector{T}
end

# I use StridedVector here because it allows for views/references,
# so you can do shallow copies: same pointed-to data. Better way?
"""Constructor for advanced usage: col and row here are 0-indexed CSC data.
Could allow for shared  `_colptr`, `_rowval`, `_nzval` between multiple
structs via views or references. Currently unused."""
function AASparseMatrix(n::Int, m::Int,
            col::StridedVector{Clong}, row::StridedVector{Cint},
            data::StridedVector{T},
            attributes::att_type = ATT_ORDINARY) where T<:vTypes
    (stride(col, 1) != 1 || stride(row, 1) != 1 || stride(data, 1) != 1) &&
        throw(ArgumentError("col, row, and data must have stride 1"))
    # I'm assuming here that pointer(row) == pointer(_row_inds),
    # ie that col, row, and data are passed by reference, not by value.
    s = SparseMatrixStructure(n, m, pointer(col),
                    pointer(row), attributes, 1)
    m = SparseMatrix(s, pointer(data))
    return AASparseMatrix(m, col, row, data)
end

function AASparseMatrix(sparseM::SparseMatrixCSC{T, Int64},
                        attributes::att_type = ATT_ORDINARY) where T<:vTypes
    if issymmetric(sparseM) && attributes == ATT_ORDINARY
        return AASparseMatrix(tril(sparseM), ATT_SYMMETRIC | ATT_LOWER_TRIANGLE)
    elseif (istril(sparseM) || istriu(sparseM)) && attributes == ATT_ORDINARY
        attributes = istril(sparseM) ? ATT_TRI_LOWER : ATT_TRI_UPPER
    end
    if attributes in (ATT_TRI_LOWER, ATT_TRI_UPPER) &&
                    all(diag(sparseM) .== one(eltype(sparseM)))
        attributes |= ATT_UNIT_TRIANGULAR
    end
    c = Clong.(sparseM.colptr .+ -1)
    r = Cint.(sparseM.rowval .+ -1)
    vals = copy(sparseM.nzval)
    return AASparseMatrix(size(sparseM)..., c, r, vals, attributes)
end

Base.size(M::AASparseMatrix) = (M.matrix.structure.rowCount,
                                    M.matrix.structure.columnCount)
Base.eltype(M::AASparseMatrix) = eltype(M._nzval)
LinearAlgebra.issymmetric(M::AASparseMatrix) = (M.matrix.structure.attributes &
                                                ATT_KIND_MASK) == ATT_SYMMETRIC
istri(M::AASparseMatrix) = (M.matrix.structure.attributes
                                    & ATT_KIND_MASK) == ATT_TRIANGULAR
LinearAlgebra.istriu(M::AASparseMatrix) = istri(M) && (MM.matrix.structure.attributes
                                        & ATT_TRIANGLE_MASK == ATT_LOWER_TRIANGLE)
LinearAlgebra.istril(M::AASparseMatrix) = istri(M) && (MM.matrix.structure.attributes
                                        & ATT_TRIANGLE_MASK == ATT_UPPER_TRIANGLE)

function Base.getindex(M::AASparseMatrix, i::Int, j::Int)
    ((size(M)[1] >= i >= 1) && (size(M)[2] >= j >= 1)) || throw(BoundsError(M, (i, j)))
    (startCol, endCol) = (M._colptr[j], M._colptr[j+1]-1) .+ 1
    rowsInCol = @view M._rowval[startCol:endCol]
    ind = searchsortedfirst(rowsInCol, i-1)
    if ind <= length(rowsInCol) && rowsInCol[ind] == i-1
        return M._nzval[startCol+ind-1]
    end
    return zero(eltype(M))
end

function Base.getindex(M::AASparseMatrix, i::Int)
    1 <= i <= size(M)[1]*size(M)[2] || throw(BoundsError(M, i))
    return M[(i-1) % size(M)[1] + 1, div(i-1, size(M)[1]) + 1]
end
# Creates a new structure, referring to the same data,
# but with the transpose flag (in attributes) flipped.
# TODO: untested
Base.transpose(M::AASparseMatrix) = AASparseMatrix(SparseGetTranspose(M.matrix),
                        M._colptr, M._rowval, M._nzval)

function Base.:(*)(A::AASparseMatrix{T}, x::StridedVecOrMat{T}) where T<:vTypes
    size(x)[1] == size(A)[2] || throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(A)[2]) and $(size(x, 1))"))
    y = Array{T}(undef, size(A)[1], size(x)[2:end]...)
    SparseMultiply(A.matrix, x, y)
    return y
end

function Base.:(*)(alpha::T, A::AASparseMatrix{T},
                            x::StridedVecOrMat{T}) where T<:vTypes
    size(x)[1] == size(A)[2] || throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(A)[2]) and $(size(x, 1))"))
    y = Array{T}(undef, size(A)[1], size(x)[2:end]...)
    SparseMultiply(alpha, A.matrix, x, y)
    return y
end

"""
Computes y += A*x in place. Note that this modifies its LAST argument.
"""
function muladd!(A::AASparseMatrix{T}, x::StridedVecOrMat{T},
                    y::StridedVecOrMat{T}) where T<:vTypes
    (size(x) == size(y) && size(x)[1] == size(A)[2]) || throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(A)[2]) and $(size(x, 1))"))
    SparseMultiplyAdd(A.matrix, x, y)
end

"""
Computes y += alpha*A*x in place. Note that this modifies its LAST argument.
"""
function muladd!(alpha::T, A::AASparseMatrix{T},
                x::StridedVecOrMat{T}, y::StridedVecOrMat{T}) where T<:vTypes
    (size(x) == size(y) && size(x)[1] == size(A)[2]) || throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(A)[2]) and $(size(x, 1))"))
    SparseMultiplyAdd(alpha, A.matrix, x, y)
end
