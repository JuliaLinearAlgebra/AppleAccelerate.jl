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
    SparseScalingDefault = 0
    SparseScalingUser = 1
    SparseScalingEquilibriationInf = 2
    # macOS 26+. Hungarian-and-ordering is only valid in a combined
    # symbolic+numeric SparseFactor call, and only for LU. Untested.
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

# Julia can't represent C bitfields directly, so we work at the bit level: define
# constants and shift to appropriate offsets. In the header, 2 separate structs:
# SparseAttributes_t and SparseAttributesComplex_t, but similar enough to treat together.
const att_type = Cuint

# Kind takes up 1 additional bit for complex (Hermitian).
const _WIDTH_KIND_REAL    = 2
const _WIDTH_KIND_COMPLEX = 3

# Bit positions.
const _SHIFT_TRANSPOSE           = 0
const _SHIFT_TRIANGLE            = 1
const _SHIFT_KIND                = 2
const _SHIFT_CONJUGATE_TRANSPOSE = _SHIFT_KIND + _WIDTH_KIND_COMPLEX
const _SHIFT_ALLOCATED           = 15

# Bitfields, lsb to msb:
# (1) bool transpose: 1 bit
const ATT_TRANSPOSE           = att_type(1) << _SHIFT_TRANSPOSE

# (2) SparseTriangle_t: 1 bit
const ATT_UPPER_TRIANGLE      = att_type(0) << _SHIFT_TRIANGLE
const ATT_LOWER_TRIANGLE      = att_type(1) << _SHIFT_TRIANGLE

# (3) SparseKind_t: 2 bits for real, 3 bits for complex (to accommodate Hermitian)
#     Ordinary=0, Triangular=1, UnitTriangular=2, Symmetric=3, Hermitian=7 (complex only).
const ATT_ORDINARY        = att_type(0) << _SHIFT_KIND
const ATT_TRIANGULAR      = att_type(1) << _SHIFT_KIND
const ATT_UNIT_TRIANGULAR = att_type(2) << _SHIFT_KIND
const ATT_SYMMETRIC       = att_type(3) << _SHIFT_KIND
const ATT_HERMITIAN       = att_type(7) << _SHIFT_KIND

# (4) bool conjugate_transpose: 1 bit, complex only
const ATT_CONJUGATE_TRANSPOSE = att_type(1) << _SHIFT_CONJUGATE_TRANSPOSE

# _reserved: middle bits "for future expansion." skip.

# (5) bool _allocatedBySparse: highest 1 bit.
const ATT_ALLOCATED_BY_SPARSE = att_type(1) << _SHIFT_ALLOCATED

# Masks for extracting fields
const ATT_TRIANGLE_MASK     = att_type(1) << _SHIFT_TRIANGLE
const ATT_KIND_MASK         = ((att_type(1) << _WIDTH_KIND_REAL)    - 1) << _SHIFT_KIND
const ATT_KIND_MASK_COMPLEX = ((att_type(1) << _WIDTH_KIND_COMPLEX) - 1) << _SHIFT_KIND
attr_mask(::Type{<:Complex}) = ATT_KIND_MASK_COMPLEX
attr_mask(::Type{<:Real})    = ATT_KIND_MASK

# Convenience combinations
const ATT_TRI_LOWER = ATT_TRIANGULAR | ATT_LOWER_TRIANGLE
const ATT_TRI_UPPER = ATT_TRIANGULAR | ATT_UPPER_TRIANGLE

const vTypes = Union{Cfloat, Cdouble, ComplexF32, ComplexF64}

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
# For complex types, use the real-part's eps and precision tier.
SparseNumericFactorOptions(T::Type) = SparseNumericFactorOptions(
    SparseDefaultControl,
    SparseScalingDefault,
    C_NULL,
    real(T) == Float32 ? 0.1 : 0.01,
    eps(real(T))*1e-4
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

# complex variants of SparseMultiply (macOS 15.5+). Note: scaled-multiply mangling
# uses Cf/Cd (not f/d) since the scalar is float/double complex.
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS0_,
    ComplexF32, Cvoid, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS0_,
    ComplexF64, Cvoid, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply26SparseMatrix_Complex_Float25DenseVector_Complex_FloatS0_,
    ComplexF32, Cvoid, SparseMatrix, StridedVector, StridedVector)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiply27SparseMatrix_Complex_Double26DenseVector_Complex_DoubleS0_,
    ComplexF64, Cvoid, SparseMatrix, StridedVector, StridedVector)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyCf26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS1_,
    ComplexF32, Cvoid, ComplexF32, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyCd27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS1_,
    ComplexF64, Cvoid, ComplexF64, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyCf26SparseMatrix_Complex_Float25DenseVector_Complex_FloatS1_,
    ComplexF32, Cvoid, ComplexF32, SparseMatrix, StridedVector, StridedVector)
@generateDemangled(SparseMultiply,
    :_Z14SparseMultiplyCd27SparseMatrix_Complex_Double26DenseVector_Complex_DoubleS1_,
    ComplexF64, Cvoid, ComplexF64, SparseMatrix, StridedVector, StridedVector)

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

# complex variants of SparseMultiplyAdd (macOS 15.5+). Scaled mangling uses Cf/Cd.
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS0_,
    ComplexF32, Cvoid, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS0_,
    ComplexF64, Cvoid, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd26SparseMatrix_Complex_Float25DenseVector_Complex_FloatS0_,
    ComplexF32, Cvoid, SparseMatrix, StridedVector, StridedVector)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAdd27SparseMatrix_Complex_Double26DenseVector_Complex_DoubleS0_,
    ComplexF64, Cvoid, SparseMatrix, StridedVector, StridedVector)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddCf26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS1_,
    ComplexF32, Cvoid, ComplexF32, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddCd27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS1_,
    ComplexF64, Cvoid, ComplexF64, SparseMatrix, StridedMatrix, StridedMatrix)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddCf26SparseMatrix_Complex_Float25DenseVector_Complex_FloatS1_,
    ComplexF32, Cvoid, ComplexF32, SparseMatrix, StridedVector, StridedVector)
@generateDemangled(SparseMultiplyAdd,
    :_Z17SparseMultiplyAddCd27SparseMatrix_Complex_Double26DenseVector_Complex_DoubleS1_,
    ComplexF64, Cvoid, ComplexF64, SparseMatrix, StridedVector, StridedVector)

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

# complex variants of SparseGetTranspose (macOS 15.5+)
@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose26SparseMatrix_Complex_Float,
    ComplexF32, SparseMatrix, SparseMatrix)
@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose27SparseMatrix_Complex_Double,
    ComplexF64, SparseMatrix, SparseMatrix)
@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose39SparseOpaqueFactorization_Complex_Float,
    ComplexF32, SparseOpaqueFactorization, SparseOpaqueFactorization)
@generateDemangled(SparseGetTranspose,
    :_Z18SparseGetTranspose40SparseOpaqueFactorization_Complex_Double,
    ComplexF64, SparseOpaqueFactorization, SparseOpaqueFactorization)

# conjugate-transpose: complex-only (no real analogue, since for real it == transpose).
# Flips the conjugate_transpose attribute bit; underlying CSC data unchanged.
@generateDemangled(SparseGetConjugateTranspose,
    :_Z27SparseGetConjugateTranspose26SparseMatrix_Complex_Float,
    ComplexF32, SparseMatrix, SparseMatrix)
@generateDemangled(SparseGetConjugateTranspose,
    :_Z27SparseGetConjugateTranspose27SparseMatrix_Complex_Double,
    ComplexF64, SparseMatrix, SparseMatrix)
@generateDemangled(SparseGetConjugateTranspose,
    :_Z27SparseGetConjugateTranspose39SparseOpaqueFactorization_Complex_Float,
    ComplexF32, SparseOpaqueFactorization, SparseOpaqueFactorization)
@generateDemangled(SparseGetConjugateTranspose,
    :_Z27SparseGetConjugateTranspose40SparseOpaqueFactorization_Complex_Double,
    ComplexF64, SparseOpaqueFactorization, SparseOpaqueFactorization)

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

# complex variants of SparseCleanup (macOS 15.5+)
@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup26SparseMatrix_Complex_Float,
    ComplexF32, Cvoid, SparseMatrix)
@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup27SparseMatrix_Complex_Double,
    ComplexF64, Cvoid, SparseMatrix)
@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup39SparseOpaqueFactorization_Complex_Float,
    ComplexF32, Cvoid, SparseOpaqueFactorization)
@generateDemangled(SparseCleanup,
    :_Z13SparseCleanup40SparseOpaqueFactorization_Complex_Double,
    ComplexF64, Cvoid, SparseOpaqueFactorization)

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

# complex variants of _SparseFactorNoErrors_inner (macOS 15.5+)
@generateDemangled(_SparseFactorNoErrors_inner,
    :_Z12SparseFactorh26SparseMatrix_Complex_Float,
    ComplexF32, SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix)
@generateDemangled(_SparseFactorNoErrors_inner,
    :_Z12SparseFactorh27SparseMatrix_Complex_Double,
    ComplexF64, SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix)

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

# complex variants of _SparseFactor_inner (macOS 15.5+)
@generateDemangled(_SparseFactor_inner,
    :_Z12SparseFactorh26SparseMatrix_Complex_Float27SparseSymbolicFactorOptions26SparseNumericFactorOptions,
    ComplexF32, SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix, SparseSymbolicFactorOptions, SparseNumericFactorOptions)
@generateDemangled(_SparseFactor_inner,
    :_Z12SparseFactorh27SparseMatrix_Complex_Double27SparseSymbolicFactorOptions26SparseNumericFactorOptions,
    ComplexF64, SparseOpaqueFactorization,
    SparseFactorization_t, SparseMatrix, SparseSymbolicFactorOptions, SparseNumericFactorOptions)

function SparseFactor(arg1::SparseFactorization_t,
                        arg2::SparseMatrix{T},
                        arg3::SparseSymbolicFactorOptions,
                        arg4::SparseNumericFactorOptions) where T <: vTypes
    arg1 != SparseFactorizationTBD || throw(ArgumentError("Factorization type must be specified"))
    _SparseFactor_inner(arg1, arg2, arg3, arg4)
end

# TODO: SparseSolve variants that take a trailing `void*` workspace pointer, for repeated
# factorizations. would require adding a workspace buffer to AAFactorization.

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

# complex variants of in-place matrix solve (macOS 15.5+)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve39SparseOpaqueFactorization_Complex_Float25DenseMatrix_Complex_Float,
    ComplexF32, Cvoid, SparseOpaqueFactorization, StridedMatrix)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve40SparseOpaqueFactorization_Complex_Double26DenseMatrix_Complex_Double,
    ComplexF64, Cvoid, SparseOpaqueFactorization, StridedMatrix)

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

# complex variants of solve-with-output, matrix RHS (macOS 15.5+)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve39SparseOpaqueFactorization_Complex_Float25DenseMatrix_Complex_FloatS0_,
    ComplexF32, Cvoid, SparseOpaqueFactorization, StridedMatrix, StridedMatrix)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve40SparseOpaqueFactorization_Complex_Double26DenseMatrix_Complex_DoubleS0_,
    ComplexF64, Cvoid, SparseOpaqueFactorization, StridedMatrix, StridedMatrix)

# in-place solve, vector RHS. The resize call prevents me from doing @generateDemangled.
const _SOLVE_VEC_INPLACE_SYM = Dict(
    Cfloat     => :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseVector_Float,
    Cdouble    => :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseVector_Double,
    ComplexF32 => :_Z11SparseSolve39SparseOpaqueFactorization_Complex_Float25DenseVector_Complex_Float,
    ComplexF64 => :_Z11SparseSolve40SparseOpaqueFactorization_Complex_Double26DenseVector_Complex_Double,
)
for (T, sym) in _SOLVE_VEC_INPLACE_SYM
    @eval function SparseSolve(arg1::SparseOpaqueFactorization{$T},
                                arg2::StridedVector{$T})
        nrows = Int(arg1.symbolicFactorization.rowCount)
        ncols = Int(arg1.symbolicFactorization.columnCount)
        # libSparse reads the RHS and writes the solution into the same buffer,
        # which must therefore hold at least max(rows, cols) elements. Validate
        # before the ccall so an undersized buffer can't be written out of bounds.
        length(arg2) >= max(nrows, ncols) ||
            throw(DimensionMismatch(
                "in-place vector solve needs a buffer of at least " *
                "max(rows, cols) = $(max(nrows, ncols)) elements; got $(length(arg2))"))
        @ccall LIBSPARSE.$sym(arg1::SparseOpaqueFactorization{$T},
                              arg2::DenseVector{$T})::Cvoid
        if length(arg2) != ncols
            # Non-square system: trim the buffer down to the solution length.
            # Only a real Vector can be resized; a view cannot.
            arg2 isa Vector ||
                throw(ArgumentError(
                    "in-place solve of a non-square system needs a resizable " *
                    "Vector RHS (got a $(typeof(arg2))); use the allocating " *
                    "`solve`/`ldiv!(x, f, b)` instead"))
            resize!(arg2, ncols)
        end
        return arg2
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

# complex variants of solve-with-output, vector RHS (macOS 15.5+)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve39SparseOpaqueFactorization_Complex_Float25DenseVector_Complex_FloatS0_,
    ComplexF32, Cvoid, SparseOpaqueFactorization, StridedVector, StridedVector)
@generateDemangled(SparseSolve,
    :_Z11SparseSolve40SparseOpaqueFactorization_Complex_Double26DenseVector_Complex_DoubleS0_,
    ComplexF64, Cvoid, SparseOpaqueFactorization, StridedVector, StridedVector)

# calls to SparseFactor default to the version with error handling.
SparseFactor(arg1::SparseFactorization_t,
                        arg2::SparseMatrix{T},
                        noErrors::Bool = false) where T <: vTypes =
                noErrors ? SparseFactorNoErrors(arg1, arg2) :
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
# TODO: libSparse also exposes symbolic-only SparseFactor variants taking
# `SparseMatrixStructureComplex` (e.g. _Z12SparseFactorh28SparseMatrixStructureComplex)
# for complex matrices. Not wrapped here — the symbolic-only path is barely
# used today, and adding it would require introducing a complex-attributes
# variant of SparseMatrixStructure.
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
# Refactorization (numeric refactor reusing a symbolic factorization)
# ============================================================
#
# `SparseRefactor` in Apple's header is an inline C++ helper that dispatches on
# the symbolic factorization's `type` and calls one of the exported
# `_SparseRefactor{Symmetric,QR,LU}_{Float,Double,...}` symbols, passing a
# scratch workspace of `symbolicFactorization.workspaceSize`. We replicate that
# dispatch in Julia. The numeric factorization is updated *in place*: given a
# new SparseMatrix with the *same sparsity pattern* as the one originally
# factored, it recomputes the numeric factors without redoing the (expensive)
# symbolic analysis/ordering. This is the workhorse for repeated solves where
# only the matrix values change (e.g. Newton iterations, time stepping).
#
# These low-level symbols take the matrix by value, a pointer to the opaque
# factorization (mutated in place), a pointer to numeric options, and a
# workspace pointer.
const _REFACTOR_SYM = Dict(
    # (T, kind) => (symbol)
    (Cfloat,     :Symmetric) => :_SparseRefactorSymmetric_Float,
    (Cdouble,    :Symmetric) => :_SparseRefactorSymmetric_Double,
    (ComplexF32, :Symmetric) => :_SparseRefactorSymmetric_Complex_Float,
    (ComplexF64, :Symmetric) => :_SparseRefactorSymmetric_Complex_Double,
    (Cfloat,     :QR)        => :_SparseRefactorQR_Float,
    (Cdouble,    :QR)        => :_SparseRefactorQR_Double,
    (ComplexF32, :QR)        => :_SparseRefactorQR_Complex_Float,
    (ComplexF64, :QR)        => :_SparseRefactorQR_Complex_Double,
    (Cfloat,     :LU)        => :_SparseRefactorLU_Float,
    (Cdouble,    :LU)        => :_SparseRefactorLU_Double,
    (ComplexF32, :LU)        => :_SparseRefactorLU_Complex_Float,
    (ComplexF64, :LU)        => :_SparseRefactorLU_Complex_Double,
)

# Map a SparseFactorization_t (read from the symbolic factorization) to the
# refactor family. Mirrors the switch in SparseRefactor (SolveImplementationTyped.h).
function _refactor_family(t::SparseFactorization_t)
    if t in (SparseFactorizationQR, SparseFactorizationCholeskyAtA)
        return :QR
    elseif t in (SparseFactorizationLU, SparseFactorizationLUUnpivoted,
                 SparseFactorizationLUSPP, SparseFactorizationLUTPP)
        return :LU
    else
        # Cholesky and the LDLT family all go through the "symmetric" refactor.
        return :Symmetric
    end
end

# Per-type workspace size of the symbolic factorization.
_refactor_workspace_size(s::SparseOpaqueSymbolicFactorization, ::Type{<:Union{Cfloat,ComplexF32}}) =
    Int(s.workspaceSize_Float)
_refactor_workspace_size(s::SparseOpaqueSymbolicFactorization, ::Type{<:Union{Cdouble,ComplexF64}}) =
    Int(s.workspaceSize_Double)

# Low-level refactor: recompute the numeric factorization of `fact` in place
# using the values of `matrix`. `matrix` must share the sparsity pattern of the
# matrix `fact` was originally built from. Returns the (mutated) factorization.
for ((T, kind), sym) in _REFACTOR_SYM
    @eval function _sparse_refactor!(fact::Base.RefValue{SparseOpaqueFactorization{$T}},
                                     matrix::SparseMatrix{$T},
                                     nfopts::Base.RefValue{SparseNumericFactorOptions},
                                     ::Val{$(QuoteNode(kind))})
        symb = fact[].symbolicFactorization
        wsize = _refactor_workspace_size(symb, $T)
        workspace = Libc.malloc(wsize)
        workspace != C_NULL || throw(OutOfMemoryError())
        try
            # The plain-C `_SparseRefactor*` symbols take the matrix BY POINTER
            # (`SparseMatrix_Double *`), unlike the C++ `SparseFactor` entry point
            # which takes it by value. Passing the struct by value misaligns the
            # x86_64 SysV argument registers (it lands on the stack), which crashed
            # libSparse with SIGILL on Intel macOS while happening to work on arm64
            # (large structs are passed indirectly there). Pass via `Ref` so ccall
            # hands over a GC-rooted pointer matching the generated ABI.
            @ccall LIBSPARSE.$sym(matrix::Ref{SparseMatrix{$T}},
                                  fact::Ptr{SparseOpaqueFactorization{$T}},
                                  nfopts::Ptr{SparseNumericFactorOptions},
                                  workspace::Ptr{Cvoid})::Cvoid
        finally
            Libc.free(workspace)
        end
        return fact[]
    end
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
    if attributes == ATT_ORDINARY
        if ishermitian(sparseM)
            kind = T <: Complex ? ATT_HERMITIAN : ATT_SYMMETRIC
            return AASparseMatrix(tril(sparseM), kind | ATT_LOWER_TRIANGLE)
        elseif istril(sparseM) || istriu(sparseM)
            attributes = istril(sparseM) ? ATT_TRI_LOWER : ATT_TRI_UPPER
        end
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

# Apple stores the underlying CSC dimensions and uses the transpose attribute
# bit to "implicitly transpose" the matrix on subsequent ops. Julia callers
# expect `size(transpose(A))` to reverse the dimensions, so we honor the bit
# here — matching the contract of `LinearAlgebra.Adjoint` / `Transpose`.
function Base.size(M::AASparseMatrix)
    nrow = Int(M.matrix.structure.rowCount)
    ncol = Int(M.matrix.structure.columnCount)
    return (M.matrix.structure.attributes & ATT_TRANSPOSE) != 0 ?
        (ncol, nrow) : (nrow, ncol)
end

LinearAlgebra.issymmetric(M::AASparseMatrix{T}) where T = (M.matrix.structure.attributes &
                                                attr_mask(T)) == ATT_SYMMETRIC
LinearAlgebra.ishermitian(M::AASparseMatrix{<:Complex}) = (M.matrix.structure.attributes &
                                                ATT_KIND_MASK_COMPLEX) == ATT_HERMITIAN
LinearAlgebra.ishermitian(M::AASparseMatrix{<:Real}) = issymmetric(M)
istri(M::AASparseMatrix{T}) where T = (M.matrix.structure.attributes
                                    & attr_mask(T)) == ATT_TRIANGULAR
LinearAlgebra.istriu(M::AASparseMatrix) = istri(M) && (M.matrix.structure.attributes &
                                        ATT_TRIANGLE_MASK == ATT_UPPER_TRIANGLE)
LinearAlgebra.istril(M::AASparseMatrix) = istri(M) && (M.matrix.structure.attributes &
                                        ATT_TRIANGLE_MASK == ATT_LOWER_TRIANGLE)

# Both bits encode adjoint; ct alone is the "no meaning" state that
# libSparse treats as identity (reachable via adjoint(adjoint(M))).
_maybe_conjugate(attr::att_type, ::Type{<:Complex}) =
    (attr & (ATT_TRANSPOSE | ATT_CONJUGATE_TRANSPOSE)) ==
    (ATT_TRANSPOSE | ATT_CONJUGATE_TRANSPOSE)
_maybe_conjugate(::att_type, ::Type{<:Real}) = false

function Base.getindex(M::AASparseMatrix{T}, i::Int, j::Int) where T<:vTypes
    ((size(M)[1] >= i >= 1) && (size(M)[2] >= j >= 1)) || throw(BoundsError(M, (i, j)))
    # (i,j) is the *logical* index; map back to the raw CSC layout.
    attrs = M.matrix.structure.attributes
    transposed = (attrs & ATT_TRANSPOSE) != 0
    raw_i, raw_j = transposed ? (j, i) : (i, j)
    (startCol, endCol) = (M._colptr[raw_j], M._colptr[raw_j+1]-1) .+ 1
    rowsInCol = @view M._rowval[startCol:endCol]
    ind = searchsortedfirst(rowsInCol, raw_i-1)
    if ind <= length(rowsInCol) && rowsInCol[ind] == raw_i-1
        val = M._nzval[startCol+ind-1]
        return _maybe_conjugate(attrs, T) ? conj(val) : val
    end
    return zero(eltype(M))
end

# libSparse traps when transposing an adjoint (T,T): would land on conj(A).
function Base.transpose(M::AASparseMatrix)
    attrs = M.matrix.structure.attributes
    !((attrs & ATT_TRANSPOSE) != 0 && (attrs & ATT_CONJUGATE_TRANSPOSE) != 0) ||
        throw(ArgumentError("cannot transpose an already-adjointed AASparseMatrix " *
            "(libSparse cannot represent the resulting `conj(A)` view; " *
            "apply `conj` to the underlying data instead)"))
    AASparseMatrix(SparseGetTranspose(M.matrix), M._colptr, M._rowval, M._nzval)
end

# Real: same as transpose. Complex: flip ct (and toggle t); CSC data shared.
Base.adjoint(M::AASparseMatrix{<:Real}) = transpose(M)
function Base.adjoint(M::AASparseMatrix{T}) where T<:Complex
    attrs = M.matrix.structure.attributes
    # libSparse traps on adjoint of a transpose-only matrix.
    !((attrs & ATT_TRANSPOSE) != 0 && (attrs & ATT_CONJUGATE_TRANSPOSE) == 0) ||
        throw(ArgumentError("cannot adjoint an already-transposed AASparseMatrix " *
            "(libSparse cannot represent the resulting `conj(A)` view; " *
            "apply `conj` to the underlying data instead)"))
    AASparseMatrix(SparseGetConjugateTranspose(M.matrix),
                   M._colptr, M._rowval, M._nzval)
end

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
    (size(x, 1) == size(A)[2] && size(y, 1) == size(A)[1] &&
        size(x)[2:end] == size(y)[2:end]) || throw(DimensionMismatch(
        "Dimension mismatch in y += A*x: A is $(size(A)), "
        * "x is $(size(x)), y is $(size(y))"))
    SparseMultiplyAdd(A.matrix, x, y)
end

"""
Computes y += alpha*A*x in place. Note that this modifies its LAST argument.
"""
function muladd!(alpha::T, A::AASparseMatrix{T},
                x::StridedVecOrMat{T}, y::StridedVecOrMat{T}) where T<:vTypes
    (size(x, 1) == size(A)[2] && size(y, 1) == size(A)[1] &&
        size(x)[2:end] == size(y)[2:end]) || throw(DimensionMismatch(
        "Dimension mismatch in y += alpha*A*x: A is $(size(A)), "
        * "x is $(size(x)), y is $(size(y))"))
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
            kind = ishermitian(aa_fact.matrixObj) ? SparseFactorizationCholesky :
                   (nrow == ncol && lu_ok)        ? SparseFactorizationLU :
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

refactor!(aa_fact::AAFactorization{T}, A::SparseMatrixCSC{T, Int64}) where T<:vTypes =
    refactor!(aa_fact, AASparseMatrix(A))

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
