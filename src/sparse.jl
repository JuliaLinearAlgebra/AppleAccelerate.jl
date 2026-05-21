## sparse.jl ##

using SparseArrays
using Libdl

# To find most recent header:
# "$(xcrun --show-sdk-path)/System/Library/Frameworks/Accelerate.framework/
# Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/Sparse/Solve.h"

# unused. But if I can find a library that supports C-style structs
# of packed bitflags with enum fields, then I'll want them.
#=
@enum SparseTriangle_t::UInt8 begin
    SparseUpperTriangle = 0
    SparseLowerTriangle = 1
end

@enum SparseKind_t::UInt32 begin
    SparseOrdinary = 0
    SparseTriangular = 1
    SparseUnitTriangular = 2
    SparseSymmetric = 3
end=#

@enum SparseFactorization_t::UInt8 begin
    SparseFactorizationCholesky = 0
    SparseFactorizationLDLT = 1
    SparseFactorizationLDLTUnpivoted = 2
    SparseFactorizationLDLTSBK = 3
    SparseFactorizationLDLTTPP = 4
    SparseFactorizationQR = 40
    SparseFactorizationCholeskyAtA = 41
    # LU variants require macOS 15.5+. Calling them on older versions returns
    # SparseParameterError from libSparse, which factor!'s status check surfaces.
    SparseFactorizationLU = 80
    SparseFactorizationLUUnpivoted = 81
    SparseFactorizationLUSPP = 82
    SparseFactorizationLUTPP = 83
    SparseFactorizationTBD = 64 # my own addition.
end

@enum SparseOrder_t::UInt8 begin
    SparseOrderDefault = 0
    SparseOrderUser = 1
    SparseOrderAMD = 2
    SparseOrderMetis = 3
    SparseOrderCOLAMD = 4
    SparseOrderMTMetis = 5 # macOS 26+
end

@enum SparseScaling_t::UInt8 begin
    SpraseScalingDefault = 0
    SparseScalingUser = 1
    SparseScalingEquilibriationInf = 2
    # macOS 26+. Hungarian-and-ordering is only valid in a combined
    # symbolic+numeric SparseFactor call, and only for LU.
    SparseScalingHungarianScalingOnly = 3
    SparseScalingHungarianScalingAndOrdering = 4
end

@enum SparseStatus_t::Int32 begin
    SparseStatusOk = 0
    SparseStatusFailed = -1
    SparseMatrixIsSingular = -2
    SparseInternalError = -3
    SparseParameterError = -4
    SparseYetToBeFactored = -5 # my own addition.
    SparseStatusReleased = -2147483647
end
# Apple renamed SparseStatusFailed to SparseFactorizationFailed in newer SDKs.
const SparseFactorizationFailed = SparseStatusFailed

@enum SparseControl_t::UInt32 begin
    SparseDefaultControl = 0
end

# can't implement SparseAttributes directly. Workaround:
const att_type = Cuint
const ATT_TRANSPOSE = att_type(1)
const ATT_UPPER_TRIANGLE = att_type(0)
const ATT_LOWER_TRIANGLE = att_type(2)
const ATT_ORDINARY = att_type(0)
const ATT_TRIANGULAR = att_type(4)
const ATT_UNIT_TRIANGULAR = att_type(8)
const ATT_SYMMETRIC = att_type(12)
const ATT_ALLOCATED_BY_SPARSE = att_type(1) << 15
const ATT_TRI_LOWER = ATT_TRIANGULAR | ATT_LOWER_TRIANGLE
const ATT_TRI_UPPER = ATT_TRIANGULAR | ATT_UPPER_TRIANGLE
const ATT_KIND_MASK = att_type(12)
const ATT_TRIANGLE_MASK = att_type(4)

const vTypes = Union{Cfloat, Cdouble}

struct SparseMatrixStructure
    rowCount::Cint
    columnCount::Cint
    columnStarts::Ptr{Clong} # Apple uses different types for indices
    rowIndices::Ptr{Cint} # whereas SparseMatrixCSC uses the same
    attributes::att_type
    blockSize::UInt8
end

struct SparseNumericFactorOptions
    control::SparseControl_t
    scalingMethod::SparseScaling_t
    scaling::Ptr{Cvoid}
    pivotTolerance::Float64
    zeroTolerance::Float64
end

# default tolerance parameters found in SolveImplementation.h header.
SparseNumericFactorOptions(T::Type) = SparseNumericFactorOptions(
    SparseDefaultControl,
    SpraseScalingDefault,
    C_NULL,
    T == Cfloat ? 0.1 : 0.01,
    eps(T)*1e-4
)

# Note: to use system defaults, values of malloc/free should be
# nil, which is *not the same* as C_NULL.
struct SparseSymbolicFactorOptions
    control::SparseControl_t
    orderMethod::SparseOrder_t
    order::Ptr{Cvoid}
    ignoreRowsAndColumns::Ptr{Cvoid}
    malloc::Ptr{Cvoid} # arg: Csize_t
    free::Ptr{Cvoid} # arg: Ptr{Cvoid}
    reportError::Ptr{Cvoid} # arg: Cstring, assuming null-terminated.
end

SparseSymbolicFactorOptions() = SparseSymbolicFactorOptions(
    SparseDefaultControl,
    SparseOrderDefault,
    C_NULL,
    C_NULL,
    @cfunction(Libc.malloc, Ptr{Cvoid}, (Csize_t,)),
    @cfunction(Libc.free, Cvoid, (Ptr{Cvoid},)),
    @cfunction(text->error(unsafe_string(text)), Cvoid, (Cstring, ))
)

struct DenseVector{T<:vTypes}
    count::Cint
    data::Ptr{T}
end

struct SparseMatrix{T<:vTypes}
    structure::SparseMatrixStructure
    data::Ptr{T}
end

struct DenseMatrix{T<:vTypes}
    rowCount::Cint
    columnCount::Cint
    columnStride::Cint
    attributes::att_type
    data::Ptr{T}
end

struct SparseOpaqueSymbolicFactorization
    status::SparseStatus_t
    rowCount::Cint
    columnCount::Cint
    attributes::att_type
    blockSize::UInt8
    type::SparseFactorization_t
    factorization::Ptr{Cvoid}
    workspaceSize_Float::Csize_t
    workspaceSize_Double::Csize_t
    factorSize_Float::Csize_t
    factorSize_Double::Csize_t
end

SparseOpaqueSymbolicFactorization() = SparseOpaqueSymbolicFactorization(
    SparseYetToBeFactored,
    0, 0,
    ATT_ORDINARY,
    0,
    SparseFactorizationQR,
    C_NULL,
    0, 0, 0, 0
)

# T isn't used in struct body, but I need it for type dispatch.
# The header has _Double, _Float variants (with identical bodies).
struct SparseOpaqueFactorization{T<:vTypes}
    status::SparseStatus_t
    attributes::att_type
    symbolicFactorization::SparseOpaqueSymbolicFactorization
    userFactorStorage::Bool
    numericFactorization::Ptr{Cvoid}
    solveWorkspaceRequiredStatic::Csize_t
    solveWorkspaceRequiredPerRHS::Csize_t
end

# placeholder object, for type stability reasons
SparseOpaqueFactorization(T::Type) = SparseOpaqueFactorization{T}(
    SparseYetToBeFactored,
    ATT_ORDINARY,
    SparseOpaqueSymbolicFactorization(),
    false,
    C_NULL,
    0, 0
)

# note: I haven't implemented anything involving Subfactor, Preconditioner, or IterativeMethod
const LIBSPARSE = "/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/"*
                    "vecLib.framework/libSparse.dylib"

Base.cconvert(::Type{DenseMatrix{T}}, m::StridedMatrix{T}) where T<:vTypes = m

function Base.unsafe_convert(::Type{DenseMatrix{T}}, m::StridedMatrix{T}) where T<:vTypes
    stride(m, 1) == 1 || throw(ArgumentError("Stride of first dimension must be 1"))
    # hard-coded ATT_ORDINARY
    return DenseMatrix{T}(size(m)..., stride(m, 2), ATT_ORDINARY, pointer(m))
end

Base.cconvert(::Type{DenseVector{T}}, v::StridedVector{T}) where T<:vTypes = v

function Base.unsafe_convert(::Type{DenseVector{T}}, v::StridedVector{T}) where T<:vTypes
    stride(v, 1) == 1 || throw(ArgumentError("Stride of first dimension must be 1"))
    return DenseVector{T}(size(v)[1], pointer(v))
end

# generates jlName([args])::retType = @ccall libSparse.cName([args])::retType
# where [args] = arg1::jlArgTypes[1], etc.
# Omit parameters for parametrized types, and instead pass as param.
# e.g. don't pass SparseMatrix{Cfloat} in jlArgTypes; just SparseMatrix, and param = Cfloat.
macro generateDemangled(jlName, cName, param, retType, jlArgTypes...)
    # all the parameterized types that we'll be using as arguments,
    # where we want param to be inserted as parameter.
    local ALL_PARAM_TYPES = [:(StridedVector), :(StridedMatrix),
                                :(SparseMatrix), :(SparseOpaqueFactorization), :(Ptr)]
    # all the enums we'll be using as arguments.
    local ALL_ENUMS = [:(SparseFactorization_t), :(SparseOrder_t), :(SparseScaling_t),
                                                    :(SparseStatus_t), :(SparseControl_t)]
    # Assemble arg/ret types: put param inside type (if parametrized type),
    # and convert to type for @ccall if different from jl type.
    local jlArgTypesWParams = Vector{Any}()
    local cArgTypesWParams = Vector{Any}()
    for T in jlArgTypes
        if T in ALL_PARAM_TYPES
            push!(jlArgTypesWParams, Expr(:curly, esc(T), esc(param)))
            if T == :(StridedVector)
                push!(cArgTypesWParams, Expr(:curly, DenseVector, esc(param)))
            elseif T == :(StridedMatrix)
                push!(cArgTypesWParams, Expr(:curly, DenseMatrix,  esc(param)))
            else
                push!(cArgTypesWParams, Expr(:curly, esc(T), esc(param)))
            end
        elseif T in ALL_ENUMS
            push!(jlArgTypesWParams, esc(T))
            push!(cArgTypesWParams, Cuint)
        else
            push!(jlArgTypesWParams, esc(T))
            push!(cArgTypesWParams, esc(T))
        end
    end
    # due to my use case, I don't need to worry about return type translation.
    local retTypeWParam = retType in ALL_PARAM_TYPES ? Expr(:curly, esc(retType), esc(param)) : esc(retType)
    # assemble julia call
    local jlArgExprs = map(enumerate(jlArgTypesWParams)) do (i, T)
        Expr(:(::), esc(Symbol("arg$i")), T)
    end
    local jlCall = Expr(:(::), Expr(:call, esc(jlName), jlArgExprs...), retTypeWParam)
    # Build ccall directly instead of using @ccall macro to avoid variable name conflicts
    local cArgNames = [Symbol("arg$i") for i in 1:length(cArgTypesWParams)]
    local LIBSPARSE = "/System/Library/Frameworks/Accelerate.framework/Versions"*
                    "/A/Frameworks/vecLib.framework/libSparse.dylib"
    # Construct ccall((:cName, LIBSPARSE), retType, (argTypes...), args...)
    local cCallExpr = Expr(:call, :ccall,
                          Expr(:tuple, esc(cName), LIBSPARSE),
                          retTypeWParam,
                          Expr(:tuple, cArgTypesWParams...),
                          cArgNames...)
    return Expr(Symbol("="), jlCall, cCallExpr)
end

# sparse * (dense matrix)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply18SparseMatrix_Float17DenseMatrix_FloatS0_,
    Cfloat,
    Cvoid,
    SparseMatrix, StridedMatrix, StridedMatrix)

@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply19SparseMatrix_Double18DenseMatrix_DoubleS0_,
    Cdouble,
    Cvoid,
    SparseMatrix, StridedMatrix, StridedMatrix)

# sparse * (dense vector)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply18SparseMatrix_Float17DenseVector_FloatS0_,
    Cfloat,
    Cvoid,
    SparseMatrix, StridedVector, StridedVector)

@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply19SparseMatrix_Double18DenseVector_DoubleS0_,
    Cdouble,
    Cvoid,
    SparseMatrix, StridedVector, StridedVector)

# scalar * sparse * (dense matrix)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyf18SparseMatrix_Float17DenseMatrix_FloatS0_,
    Cfloat,
    Cvoid,
    Cfloat, SparseMatrix, StridedMatrix, StridedMatrix)

@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyd19SparseMatrix_Double18DenseMatrix_DoubleS0_,
    Cdouble,
    Cvoid,
    Cdouble,  SparseMatrix, StridedMatrix, StridedMatrix)

# scalar * sparse * (dense vector)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyf18SparseMatrix_Float17DenseVector_FloatS0_,
    Cfloat,
    Cvoid,
    Cfloat, SparseMatrix, StridedVector, StridedVector)

@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyd19SparseMatrix_Double18DenseVector_DoubleS0_,
    Cdouble,
    Cvoid,
    Cdouble,  SparseMatrix, StridedVector, StridedVector)

# dense matrix += sparse * (dense matrix)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd18SparseMatrix_Float17DenseMatrix_FloatS0_,
    Cfloat,
    Cvoid,
    SparseMatrix, StridedMatrix, StridedMatrix)

@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd19SparseMatrix_Double18DenseMatrix_DoubleS0_,
    Cdouble,
    Cvoid,
    SparseMatrix, StridedMatrix, StridedMatrix)

# dense vector += sparse * (dense vector)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd18SparseMatrix_Float17DenseVector_FloatS0_,
    Cfloat,
    Cvoid,
    SparseMatrix, StridedVector, StridedVector)

@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd19SparseMatrix_Double18DenseVector_DoubleS0_,
    Cdouble,
    Cvoid,
    SparseMatrix, StridedVector, StridedVector)

# dense matrix += scalar * sparse * (dense matrix)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddf18SparseMatrix_Float17DenseMatrix_FloatS0_,
    Cfloat,
    Cvoid,
    Cfloat, SparseMatrix, StridedMatrix, StridedMatrix)

@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddd19SparseMatrix_Double18DenseMatrix_DoubleS0_,
    Cdouble,
    Cvoid,
    Cdouble, SparseMatrix, StridedMatrix, StridedMatrix)

# dense vector += scalar * sparse * (dense vector)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddf18SparseMatrix_Float17DenseVector_FloatS0_,
    Cfloat,
    Cvoid,
    Cfloat, SparseMatrix, StridedVector, StridedVector)

@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddd19SparseMatrix_Double18DenseVector_DoubleS0_,
    Cdouble,
    Cvoid,
    Cdouble, SparseMatrix, StridedVector, StridedVector)

# transpose of matrix
@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose18SparseMatrix_Float,
    Cfloat,
    SparseMatrix,
    SparseMatrix)

@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose19SparseMatrix_Double,
    Cdouble,
    SparseMatrix,
    SparseMatrix)
# (skipped: 2 subfactor transposes)
# transpose of opaque factorization
@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose31SparseOpaqueFactorization_Float,
    Cfloat,
    SparseOpaqueFactorization,
    SparseOpaqueFactorization)

@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose32SparseOpaqueFactorization_Double,
    Cdouble,
    SparseOpaqueFactorization,
    SparseOpaqueFactorization)

# TODO: these 4 SparseConvertFromCoord functions are unused and untested.
@generateDemangled(SparseConvertFromCoordinate,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKf,
    Cfloat,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr)

@generateDemangled(SparseConvertFromCoordinate,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKd,
    Cdouble,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr)

@generateDemangled(SparseConvertFromCoord,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKdPvS4_,
    Cfloat,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr, Ptr{Cvoid}, Ptr{Cvoid})

@generateDemangled(SparseConvertFromCoord,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKfPvS4_,
    Cdouble,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr, Ptr{Cvoid}, Ptr{Cvoid})

# cleanup SparseMatrix
@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup18SparseMatrix_Float,
    Cfloat,
    Cvoid,
    SparseMatrix)

@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup19SparseMatrix_Double,
    Cdouble,
    Cvoid,
    SparseMatrix)

# cleanup OpaqueFactorization.
@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup31SparseOpaqueFactorization_Float,
    Cfloat,
    Cvoid,
    SparseOpaqueFactorization
)

@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup32SparseOpaqueFactorization_Double,
    Cdouble,
    Cvoid,
    SparseOpaqueFactorization
)

# factor: factorization type, matrix
# additional outer function for error handling with the added TBD factorization type.
# The "NoErrors" part: if the call throws, Julia won't recognize the thrown
# object (it's objective-C, os_log_error) and just reports "illegal instruction"
@generateDemangled(_SparseFactorNoErrors_inner,
    :_Z12SparseFactorh18SparseMatrix_Float,
    Cfloat,
    SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix
)

@generateDemangled(_SparseFactorNoErrors_inner,
    :_Z12SparseFactorh19SparseMatrix_Double,
    Cdouble,
    SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix
)

function SparseFactorNoErrors(arg1::SparseFactorization_t,
                            arg2::SparseMatrix{T}) where T <: vTypes
    arg1 != SparseFactorizationTBD || throw(ArgumentError("Factorization type must be specified"))
    _SparseFactorNoErrors_inner(arg1, arg2)
end

# factor with non-default options
# additional outer function for error handling with the added TBD factorization type.
@generateDemangled(_SparseFactor_inner,
    :_Z12SparseFactorh18SparseMatrix_Float27SparseSymbolicFactorOptions26SparseNumericFactorOptions,
    Cfloat,
    SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix, SparseSymbolicFactorOptions, SparseNumericFactorOptions
)

@generateDemangled(_SparseFactor_inner,
    :_Z12SparseFactorh19SparseMatrix_Double27SparseSymbolicFactorOptions26SparseNumericFactorOptions,
    Cdouble,
    SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix, SparseSymbolicFactorOptions, SparseNumericFactorOptions
)

function SparseFactor(arg1::SparseFactorization_t,
                        arg2::SparseMatrix{T},
                        arg3::SparseSymbolicFactorOptions,
                        arg4::SparseNumericFactorOptions) where T <: vTypes
    arg1 != SparseFactorizationTBD || throw(ArgumentError("Factorization type must be specified"))
    _SparseFactor_inner(arg1, arg2, arg3, arg4)
end

# in-place solve, matrix RHS
@generateDemangled(SparseSolve,
    :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseMatrix_Float,
    Cfloat,
    Cvoid,
    SparseOpaqueFactorization, StridedMatrix
)

@generateDemangled(SparseSolve,
    :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_Double,
    Cdouble,
    Cvoid,
    SparseOpaqueFactorization, StridedMatrix
)

# solve, matrix RHS (modifies 3nd argument instead of returning)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseMatrix_FloatS0_,
    Cfloat,
    Cvoid,
    SparseOpaqueFactorization, StridedMatrix, StridedMatrix
)

@generateDemangled(SparseSolve,
    :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_DoubleS0_,
    Cdouble,
    Cvoid,
    SparseOpaqueFactorization, StridedMatrix, StridedMatrix
)

# in-place solve, vector RHS. The resize call prevents me from doing @generateDemangled.
for T in (Cdouble, Cfloat)
    local sparseSolveVecInPlace = T == Cfloat ? :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseVector_Float :
                                :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseVector_Double
    @eval function SparseSolve(arg1::SparseOpaqueFactorization{$T},
                                arg2::StridedVector{$T})
        @ccall LIBSPARSE.$sparseSolveVecInPlace(arg1::SparseOpaqueFactorization{$T},
                                            arg2::DenseVector{$T})::Cvoid
        resize!(arg2, arg1.symbolicFactorization.columnCount)
    end
end

# solve, vector RHS. modifies 3rd argument instead of returning
@generateDemangled(SparseSolve,
    :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseVector_FloatS0_,
    Cfloat,
    Cvoid,
    SparseOpaqueFactorization, StridedVector,StridedVector
)

@generateDemangled(SparseSolve,
    :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseVector_DoubleS0_,
    Cdouble,
    Cvoid,
    SparseOpaqueFactorization, StridedVector,StridedVector
)

# calls to SparseFactor default to the version with error handling.
SparseFactor(arg1::SparseFactorization_t,
                        arg2::SparseMatrix{T},
                        noErrors::Bool = false) where T <: vTypes =
                noErrors ? SparseFactor(arg1, arg2) :
                SparseFactor(arg1, arg2, SparseSymbolicFactorOptions(),
                                         SparseNumericFactorOptions(T))

# could replace with generateDemangled macrocalls, but not parametrized so not much simpler.
SparseCleanup(arg1::SparseOpaqueSymbolicFactorization) = @ccall (
        LIBSPARSE._Z13SparseCleanup33SparseOpaqueSymbolicFactorization(
                                arg1::SparseOpaqueSymbolicFactorization)::Cvoid
)

function SparseFactor(arg1::SparseFactorization_t,
        arg2::SparseMatrixStructure,
        noErrors::Bool = false)
    noErrors ? SparseFactorNoErrors(arg1, arg2) :
            SparseFactor(arg1, arg2, SparseSymbolicFactorOptions())
end

# keep just in case libSparse's default malloc/free cooperate better
# with the library code (due to being from Objective-C, not from C)
function SparseFactorNoErrors(arg1::SparseFactorization_t, arg2::SparseMatrixStructure)
    @ccall LIBSPARSE._Z12SparseFactorh21SparseMatrixStructure(
        arg1::Cuint, arg2::SparseMatrixStructure
    )::SparseOpaqueSymbolicFactorization
end

function SparseFactor(arg1::SparseFactorization_t,
            arg2::SparseMatrixStructure,
            arg3::SparseSymbolicFactorOptions)
    @ccall LIBSPARSE._Z12SparseFactorh21SparseMatrixStructure27SparseSymbolicFactorOptions(
                arg1::Cuint,
                arg2::SparseMatrixStructure,
                arg3::SparseSymbolicFactorOptions
    )::SparseOpaqueSymbolicFactorization
end

# ============================================================
# AASparseMatrix
# ============================================================

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
LinearAlgebra.istriu(M::AASparseMatrix) = istri(M) && (M.matrix.structure.attributes &
                                        ATT_LOWER_TRIANGLE == ATT_UPPER_TRIANGLE)
LinearAlgebra.istril(M::AASparseMatrix) = istri(M) && (M.matrix.structure.attributes &
                                        ATT_LOWER_TRIANGLE == ATT_LOWER_TRIANGLE)

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

# ============================================================
# AAFactorization
# ============================================================

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
"""
    factor!(f::AAFactorization, [type::SparseFactorization_t])

Explicitly compute the factorization. If `type` is not specified, Cholesky is used
for symmetric matrices and QR for non-symmetric. Called automatically by [`solve`](@ref)
if the factorization has not yet been computed.
"""
function factor!(aa_fact::AAFactorization{T},
            kind::SparseFactorization_t = SparseFactorizationTBD) where T<:vTypes
    if aa_fact._factorization.status == SparseYetToBeFactored
        if kind == SparseFactorizationTBD
            # Cholesky for symmetric; on macOS 15.5+ LU for square non-symmetric
            # (faster than QR for square solves); QR otherwise.
            nrow, ncol = size(aa_fact.matrixObj)
            lu_ok = something(_macos_version[], v"0.0.0") >= v"15.5"
            kind = issymmetric(aa_fact.matrixObj) ? SparseFactorizationCholesky :
                   (nrow == ncol && lu_ok)        ? SparseFactorizationLU :
                                                    SparseFactorizationQR
        end
        aa_fact._factorization = SparseFactor(kind, aa_fact.matrixObj.matrix)
        _libsparse_throw(aa_fact._factorization.status, "factor")
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
    size(aa_fact.matrixObj)[2] != size(b, 1) && throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(aa_fact.matrixObj)[2]) and $(size(b, 1))"))
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

LinearAlgebra.ldiv!(aa_fact::AAFactorization{T}, xb::StridedVecOrMat{T}) where T<:vTypes =
        solve!(aa_fact, xb)

function LinearAlgebra.ldiv!(x::StridedVecOrMat{T},
                            aa_fact::AAFactorization{T},
                            b::StridedVecOrMat{T}) where T<:vTypes
    size(aa_fact.matrixObj)[2] != size(b, 1) && throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(aa_fact.matrixObj)[2]) and $(size(b, 1))"))
    size(aa_fact.matrixObj)[2] != size(x, 1) && throw(DimensionMismatch(
        "Matrix and output size mismatch: got "
        * "$(size(aa_fact.matrixObj)[2]) and $(size(x, 1))"))
    factor!(aa_fact)
    SparseSolve(aa_fact._factorization, b, x)
    return x
end
