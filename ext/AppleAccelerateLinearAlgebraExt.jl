module AppleAccelerateLinearAlgebraExt

@static if Sys.isapple()

using Libdl
using LinearAlgebra
using AppleAccelerate
using AppleAccelerate: libacc, vTypes,
    AASparseMatrix, AAFactorization, SparseOpaqueFactorization, SparseMatrix,
    SparseFactorization_t, SparseFactorizationTBD,
    SparseFactorizationCholesky, SparseFactorizationLU, SparseFactorizationQR,
    SparseStatus_t, SparseStatusOk, SparseMatrixIsSingular, SparseParameterError,
    SparseStatusFailed, SparseInternalError, SparseYetToBeFactored,
    SparseStatusReleased,
    SparseFactor, SparseSolve, SparseCleanup,
    SparseNumericFactorOptions, _sparse_refactor!, _refactor_family,
    _is_symmetric_attr, _is_hermitian_attr, _is_triu_attr, _is_tril_attr,
    get_macos_version, _macos_version

import AppleAccelerate: factor!, solve, solve!, refactor!

# ===========================================================================
# BLAS/LAPACK forwarding to Accelerate (moved here from the core __init__)
# ===========================================================================

function AppleAccelerate.forward_accelerate(interface::Symbol;
                            new_lapack::Bool = interface == :ilp64,
                            clear::Bool = false,
                            verbose::Bool = false)
    kwargs = Dict{Symbol,String}()
    if new_lapack
        if interface == :ilp64
            kwargs[:suffix_hint] = "\x1a\$NEWLAPACK\$ILP64"
        else
            kwargs[:suffix_hint] = "\x1a\$NEWLAPACK"
        end
    else
        if interface == :ilp64
            throw(ArgumentError("ILP64 accelerate requires new_lapack"))
        end
    end
    BLAS.lbt_forward(libacc; clear, verbose, kwargs...)
end

function AppleAccelerate.load_accelerate(; clear::Bool = false,
                           verbose::Bool = false,
                           load_ilp64::Bool = true)
    libacc_hdl = Libdl.dlopen_e(libacc)
    if libacc_hdl == C_NULL
        return
    end

    # Check to see if we can load ILP64 symbols
    if load_ilp64 && Libdl.dlsym_e(libacc_hdl, "dgemm\$NEWLAPACK\$ILP64") == C_NULL
        error("Unable to load ILP64 interface from '$(libacc)'; you are running macOS $(get_macos_version()), you need v13.4+")
    end

    # First, load :lp64 symbols, optionally clearing the current LBT forwarding tables
    AppleAccelerate.forward_accelerate(:lp64; new_lapack=true, clear, verbose)
    if load_ilp64
        AppleAccelerate.forward_accelerate(:ilp64; new_lapack=true, verbose)
    end
end

# ===========================================================================
# AASparseMatrix predicates (delegate to core attribute-bit helpers)
# ===========================================================================

LinearAlgebra.issymmetric(M::AASparseMatrix) = _is_symmetric_attr(M)
LinearAlgebra.ishermitian(M::AASparseMatrix) = _is_hermitian_attr(M)
LinearAlgebra.istriu(M::AASparseMatrix) = _is_triu_attr(M)
LinearAlgebra.istril(M::AASparseMatrix) = _is_tril_attr(M)

# ===========================================================================
# The factorize/solve subsystem on the core `AAFactorization` type
# ===========================================================================
#
# `AAFactorization` is defined in core (src/sparse.jl) as a plain struct so that
# the dual LinearAlgebra+SparseArrays extension can add a `SparseMatrixCSC`
# constructor for it. Core deliberately does NOT subtype
# `LinearAlgebra.Factorization` (no LinearAlgebra dependency). This extension
# supplies the full `Factorization`-style API, including the `Base.:\` method
# that the `Factorization` supertype would otherwise have provided.

LinearAlgebra.factorize(A::AASparseMatrix{T}) where T<:vTypes = AAFactorization(A)

# easiest way to make this follow the defaults and naming conventions of LinearAlgebra?
"""
    factor!(f::AAFactorization, [type::SparseFactorization_t])

Explicitly compute the factorization. If `type` is not specified, the default
is chosen from the matrix's kind attribute:

  - Hermitian (real symmetric or complex Hermitian) → Cholesky
  - square, non-Hermitian → LU (macOS 15.5+; falls back to QR on older systems)
  - rectangular → QR

Called automatically by [`solve`](@ref) if the factorization has not yet been
computed.

!!! note "Complex symmetric (not Hermitian)"
    When the user leaves the factorization type unspecified and passes
    a complex symmetric (but not Hermitian) matrix, we default to LU factorization.

    On older versions of MacOS, Apple's `SparseFactorizationLDLT` errors on complex
    symmetric (non-Hermitian) matrices, reporting "Cannot
    perform Hermitian matrix factorization of non-Hermitian matrix."
    On newer versions, it accepts the call and performs a true `LDLᵀ` (not `LDLᴴ`).
    This behavior isn't documented by Apple; the exact version threshold is
    unknown.
"""
function factor!(aa_fact::AAFactorization{T},
            kind::SparseFactorization_t = SparseFactorizationTBD) where T<:vTypes
    if aa_fact._factorization.status == SparseYetToBeFactored
        if kind == SparseFactorizationTBD
            # Cholesky for Hermitian (incl. real symmetric); LU for square
            # non-Hermitian (macOS 15.5+); QR otherwise.
            nrow, ncol = size(aa_fact.matrixObj)
            lu_ok = something(_macos_version[], v"0.0.0") >= v"15.5"
            kind = _is_hermitian_attr(aa_fact.matrixObj) ? SparseFactorizationCholesky :
                   (nrow == ncol && lu_ok)               ? SparseFactorizationLU :
                                                           SparseFactorizationQR
        end
        fact = SparseFactor(kind, aa_fact.matrixObj.matrix)
        if fact.status == SparseStatusOk
            aa_fact._factorization = fact
        else
            # The factorization failed (e.g. singular matrix). Free whatever
            # SparseFactor allocated for the failed attempt, but leave the
            # original yet-to-be-factored placeholder in place so (a) the
            # finalizer never calls SparseCleanup on a non-Ok object and
            # (b) a subsequent factor! retries instead of silently reusing a
            # broken factorization.
            SparseCleanup(fact)
            _libsparse_throw(fact.status, "factor")
        end
    end
end

# Translate a libSparse status code into a typed Julia exception. SparseStatusOk
# is a no-op so callers can write `_libsparse_throw(result.status, ...)` after
# every libSparse call without branching.
function _libsparse_throw(status::SparseStatus_t, op::AbstractString)
    status == SparseStatusOk && return
    status == SparseMatrixIsSingular &&
        throw(LinearAlgebra.SingularException(0))
    status == SparseParameterError &&
        throw(ArgumentError("libSparse $(op) failed: parameter error " *
            "(possibly: matrix lacks the properties required by this factorization, " *
            "or factorization is unavailable on this macOS version)"))
    status == SparseStatusFailed &&
        error("libSparse $(op) failed (check that the matrix has the correct " *
              "properties for this factorization)")
    status == SparseInternalError &&
        error("libSparse $(op) failed: internal error")
    error("libSparse $(op) failed: status=$(status)")
end

"""
    solve(f::AAFactorization, b::StridedVecOrMat)

Solve the linear system `Ax = b` using Apple's Sparse Solvers, returning the solution `x`.
The factorization is computed lazily on the first call if not already factored.
Equivalent to `f \\ b`.
"""
function solve(aa_fact::AAFactorization{T}, b::StridedVecOrMat{T}) where T<:vTypes
    size(aa_fact.matrixObj)[1] != size(b, 1) && throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(aa_fact.matrixObj)[1]) and $(size(b, 1))"))
    factor!(aa_fact)
    x = Array{T}(undef, size(aa_fact.matrixObj)[2], size(b)[2:end]...)
    SparseSolve(aa_fact._factorization, b, x)
    return x
end

"""
    solve!(f::AAFactorization, xb::StridedVecOrMat)

Solve the linear system `Ax = b` in-place, overwriting `xb` with the solution.
On input `xb` contains the right-hand side `b`; on output it contains the solution `x`.
Equivalent to `ldiv!(f, xb)`.
"""
function solve!(aa_fact::AAFactorization{T}, xb::StridedVecOrMat{T}) where T<:vTypes
    ((xb isa StridedMatrix) &&
        (size(aa_fact.matrixObj)[1] != size(aa_fact.matrixObj)[2])) &&
        throw(ArgumentError("Can't in-place " *
                "solve: x and b are different sizes and Julia cannot resize a matrix."
            )
        )
    factor!(aa_fact)
    SparseSolve(aa_fact._factorization, xb)
    return xb # written in imitation of KLU.jl, which also returns
end

# `\` would normally come from `AAFactorization <: LinearAlgebra.Factorization`,
# but core keeps the type free of any LinearAlgebra dependency, so define it here.
Base.:\(aa_fact::AAFactorization{T}, b::StridedVecOrMat{T}) where T<:vTypes =
        solve(aa_fact, b)

LinearAlgebra.ldiv!(aa_fact::AAFactorization{T}, xb::StridedVecOrMat{T}) where T<:vTypes =
        solve!(aa_fact, xb)

function LinearAlgebra.ldiv!(x::StridedVecOrMat{T},
                            aa_fact::AAFactorization{T},
                            b::StridedVecOrMat{T}) where T<:vTypes
    size(aa_fact.matrixObj)[1] != size(b, 1) && throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(aa_fact.matrixObj)[1]) and $(size(b, 1))"))
    size(aa_fact.matrixObj)[2] != size(x, 1) && throw(DimensionMismatch(
        "Matrix and output size mismatch: got "
        * "$(size(aa_fact.matrixObj)[2]) and $(size(x, 1))"))
    factor!(aa_fact)
    SparseSolve(aa_fact._factorization, b, x)
    return x
end

# ===========================================================================
# refactor! (numeric refactor reusing a symbolic factorization)
# ===========================================================================

"""
    refactor!(f::AAFactorization, A::AASparseMatrix)
    refactor!(f::AAFactorization, A::SparseMatrixCSC)

Recompute the numeric factorization stored in `f` using the values of `A`,
**reusing the existing symbolic factorization** (the fill-reducing ordering and
sparsity analysis). `A` must have the **same sparsity pattern** as the matrix
`f` was originally factored from; only its numeric values may differ.

This is substantially cheaper than building a fresh [`AAFactorization`](@ref)
whenever a matrix changes values but not structure — the common case in Newton
iterations, implicit time stepping, and parameter sweeps.

`f` must already hold a completed factorization (call [`factor!`](@ref) or
[`solve`](@ref) at least once first); otherwise an `ArgumentError` is thrown.
The factorization object `f` is mutated in place and returned. After
`refactor!`, `f` references `A`'s data, so keep `A` alive as long as `f` is used.

!!! warning
    Apple's libSparse does not validate that the new pattern matches the old
    one. Passing a matrix with a different sparsity pattern is undefined
    behavior. The dimensions and stored nonzero count are checked here as a
    cheap guard, but an identical count with a different pattern will not be
    caught.
"""
function refactor!(aa_fact::AAFactorization{T}, A::AASparseMatrix{T}) where T<:vTypes
    status = aa_fact._factorization.status
    status == SparseStatusOk || throw(ArgumentError(
        "refactor! requires an already-computed factorization; call factor! or " *
        "solve first (current status: $status)"))
    old = aa_fact.matrixObj
    size(old) == size(A) || throw(DimensionMismatch(
        "refactor! requires matching dimensions: factored matrix is $(size(old)), " *
        "new matrix is $(size(A))"))
    length(old._nzval) == length(A._nzval) || throw(ArgumentError(
        "refactor! requires the same number of stored nonzeros (same sparsity " *
        "pattern): factored matrix has $(length(old._nzval)), new matrix has " *
        "$(length(A._nzval))"))

    family = _refactor_family(aa_fact._factorization.symbolicFactorization.type)
    nfopts = Ref(SparseNumericFactorOptions(T))
    factref = Ref(aa_fact._factorization)
    GC.@preserve A begin
        _sparse_refactor!(factref, A.matrix, nfopts, Val(family))
    end
    new_fact = factref[]
    _libsparse_throw(new_fact.status, "refactor")
    aa_fact._factorization = new_fact
    # Keep the new matrix alive: its CSC buffers now back the numeric factors'
    # input, and dropping the old wrapper avoids confusion about which values
    # the factorization reflects.
    aa_fact.matrixObj = A
    return aa_fact
end

# ===========================================================================
# Installation into the parent module + forwarding on load
# ===========================================================================

function __init__()
    # BLAS/LAPACK forwarding. Mirror the macOS-version guard that the core
    # __init__ used before forwarding was moved here.
    ver = get_macos_version()
    # dsptrf has a bug in the initial release of the $NEWLAPACK symbols in 13.3;
    # require macOS 13.4 for ILP64, a correct LAPACK, and threading APIs.
    if ver === nothing || ver < v"13.4"
        @info "AppleAccelerate.jl needs macOS 13.4 or later for BLAS/LAPACK forwarding"
        return
    end
    AppleAccelerate.load_accelerate(; clear = false, load_ilp64 = true)
end

end # @static if Sys.isapple()

end # module
