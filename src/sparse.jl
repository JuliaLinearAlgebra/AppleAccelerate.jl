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

# SparseConvertFromCoordinate: non-workspace variant (Apple allocates internally) and
# workspace variant `SparseConvertFromCoord` (caller-owned storage). The workspace
# variant backs the COO->CSC `AASparseMatrix(I, J, V, m, n)` constructor below; the
# non-workspace pair remains available but currently has no idiomatic caller.
# (Float32/Float64 dispatch on the two workspace-taking variants below was corrected
#  so Cfloat binds the PKf symbol and Cdouble the PKd symbol.)
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
    Cdouble,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr, Ptr{Cvoid}, Ptr{Cvoid})

@generateDemangled(SparseConvertFromCoord,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKfPvS4_,
    Cfloat,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr, Ptr{Cvoid}, Ptr{Cvoid})

# complex COO->CSC workspace variants (macOS 15.5+). The library mangles the
# complex value pointer as C `float _Complex`/`double _Complex` (`PKCf`/`PKCd`),
# which adds one substitution vs the real symbols, so the trailing void* reuses
# `S5_` (not `S4_`). The attribute struct is still the plain `SparseAttributes_t`
# (18); its 3-bit Hermitian kind field fits the same 32-bit word.
@generateDemangled(SparseConvertFromCoord,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKCfPvS5_,
    ComplexF32,
    SparseMatrix,
    Cint, Cint, Clong, Cuchar, att_type, Ptr{Cint}, Ptr{Cint}, Ptr, Ptr{Cvoid}, Ptr{Cvoid})

@generateDemangled(SparseConvertFromCoord,
    :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKCdPvS5_,
    ComplexF64,
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

"""
    AASparseMatrix(I, J, V, m, n, attributes = ATT_ORDINARY)

Build an `m`×`n` `AASparseMatrix` from coordinate (COO) triplets: `V[k]` is the
value at row `I[k]`, column `J[k]`, with 1-based indices as in
`SparseArrays.sparse(I, J, V, m, n)`. Triplets may be given in any order (they are
sorted internally) and duplicate coordinates are summed. Uses Accelerate's
[`SparseConvertFromCoordinate`](https://developer.apple.com/documentation/accelerate/sparseconvertfromcoordinate).

`V` must be `Float32`, `Float64`, `ComplexF32`, or `ComplexF64` (complex requires
macOS 15.5+). Pass `attributes` (e.g. `ATT_SYMMETRIC | ATT_LOWER_TRIANGLE`, or
`ATT_HERMITIAN | ATT_LOWER_TRIANGLE` for complex) to build a matrix with
symmetry/triangle structure; Accelerate coerces the input to conform.
"""
function AASparseMatrix(I::AbstractVector{<:Integer}, J::AbstractVector{<:Integer},
                        V::AbstractVector{T}, m::Integer, n::Integer,
                        attributes::att_type = ATT_ORDINARY) where T<:vTypes
    length(I) == length(J) == length(V) ||
        throw(DimensionMismatch("I, J, and V must have equal lengths " *
                                "($(length(I)), $(length(J)), $(length(V)))"))
    nnz_in = length(V)
    # Accelerate's workspace SparseConvertFromCoordinate expects the coordinates in
    # column-major order (its rowCount-sized workspace only supports a single
    # column-ordered pass); sort by (column, row) so any input order works, matching
    # `SparseArrays.sparse`. Duplicate coordinates end up adjacent and Accelerate sums them.
    perm = sortperm(collect(zip(J, I)))
    # Accelerate takes 0-based Cint coordinate indices.
    rows0 = Vector{Cint}(undef, nnz_in)
    cols0 = Vector{Cint}(undef, nnz_in)
    vals  = Vector{T}(undef, nnz_in)
    @inbounds for (dst, k) in enumerate(perm)
        1 <= I[k] <= m || throw(ArgumentError("row index I[$k] = $(I[k]) out of range 1:$m"))
        1 <= J[k] <= n || throw(ArgumentError("column index J[$k] = $(J[k]) out of range 1:$n"))
        rows0[dst] = Cint(I[k] - 1)
        cols0[dst] = Cint(J[k] - 1)
        vals[dst]  = V[k]
    end
    blockSize = 1

    # Caller-owned storage/workspace, sized per Accelerate's Solve.h contract
    # (blockCount = number of input entries, an upper bound on the deduplicated result).
    storage_bytes = 48 + (Int(n) + 1) * sizeof(Clong) + nnz_in * sizeof(Cint) +
                    nnz_in * blockSize * blockSize * sizeof(T)
    storage   = Vector{UInt8}(undef, storage_bytes)
    workspace = Vector{UInt8}(undef, Int(m) * sizeof(Cint))

    colptr = Vector{Clong}(undef, Int(n) + 1)
    local rowval::Vector{Cint}
    local nzval::Vector{T}
    GC.@preserve rows0 cols0 vals storage workspace begin
        sm = SparseConvertFromCoord(Cint(m), Cint(n), Clong(nnz_in), Cuchar(blockSize),
                                    attributes, pointer(rows0), pointer(cols0), pointer(vals),
                                    Ptr{Cvoid}(pointer(storage)), Ptr{Cvoid}(pointer(workspace)))
        # The returned CSC arrays point into `storage`; copy them into owned Vectors
        # so the AASparseMatrix follows the standard lifetime model and `storage`
        # (and its scratch) can be collected. Indices come back 0-based, which is
        # exactly what the advanced AASparseMatrix constructor expects.
        copyto!(colptr, unsafe_wrap(Array, sm.structure.columnStarts, Int(n) + 1))
        nnz_out = Int(colptr[end])
        rowval = Vector{Cint}(undef, nnz_out)
        nzval  = Vector{T}(undef, nnz_out)
        copyto!(rowval, unsafe_wrap(Array, sm.structure.rowIndices, nnz_out))
        copyto!(nzval,  unsafe_wrap(Array, sm.data, nnz_out))
    end
    return AASparseMatrix(Int(m), Int(n), colptr, rowval, nzval, attributes)
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

# Look up a (raw_i, raw_j) entry directly in the stored CSC data, returning
# `nothing` if there is no stored value at that position.
function _csc_lookup(M::AASparseMatrix{T}, raw_i::Int, raw_j::Int) where T<:vTypes
    (startCol, endCol) = (M._colptr[raw_j], M._colptr[raw_j+1]-1) .+ 1
    rowsInCol = @view M._rowval[startCol:endCol]
    ind = searchsortedfirst(rowsInCol, raw_i-1)
    if ind <= length(rowsInCol) && rowsInCol[ind] == raw_i-1
        return M._nzval[startCol+ind-1]
    end
    return nothing
end

function Base.getindex(M::AASparseMatrix{T}, i::Int, j::Int) where T<:vTypes
    ((size(M)[1] >= i >= 1) && (size(M)[2] >= j >= 1)) || throw(BoundsError(M, (i, j)))
    # (i,j) is the *logical* index; map back to the raw CSC layout.
    attrs = M.matrix.structure.attributes
    transposed = (attrs & ATT_TRANSPOSE) != 0
    raw_i, raw_j = transposed ? (j, i) : (i, j)
    val = _csc_lookup(M, raw_i, raw_j)
    if val !== nothing
        return _maybe_conjugate(attrs, T) ? conj(val) : val
    end
    # Symmetric/Hermitian matrices store only one triangle, so an element in the
    # empty triangle must be read from its mirror (j,i), conjugating for
    # Hermitian-complex.
    if issymmetric(M) || ishermitian(M)
        mirrored = _csc_lookup(M, raw_j, raw_i)
        if mirrored !== nothing
            herm = (attrs & ATT_KIND_MASK_COMPLEX) == ATT_HERMITIAN
            return herm ? conj(mirrored) : mirrored
        end
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
    # M.matrix carries raw pointers into M's CSC buffers; keep M alive for the call.
    t = GC.@preserve M SparseGetTranspose(M.matrix)
    AASparseMatrix(t, M._colptr, M._rowval, M._nzval)
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
    ct = GC.@preserve M SparseGetConjugateTranspose(M.matrix)
    AASparseMatrix(ct, M._colptr, M._rowval, M._nzval)
end

function Base.:(*)(A::AASparseMatrix{T}, x::StridedVecOrMat{T}) where T<:vTypes
    size(x)[1] == size(A)[2] || throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(A)[2]) and $(size(x, 1))"))
    y = Array{T}(undef, size(A)[1], size(x)[2:end]...)
    # A.matrix is a C value struct holding raw pointers into A's CSC buffers,
    # which nothing else roots across the ccall; keep A alive for the call.
    GC.@preserve A SparseMultiply(A.matrix, x, y)
    return y
end

function Base.:(*)(alpha::T, A::AASparseMatrix{T},
                            x::StridedVecOrMat{T}) where T<:vTypes
    size(x)[1] == size(A)[2] || throw(DimensionMismatch(
        "Matrix and right-hand side size mismatch: got "
        * "$(size(A)[2]) and $(size(x, 1))"))
    y = Array{T}(undef, size(A)[1], size(x)[2:end]...)
    GC.@preserve A SparseMultiply(alpha, A.matrix, x, y)
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
    GC.@preserve A SparseMultiplyAdd(A.matrix, x, y)
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
    GC.@preserve A SparseMultiplyAdd(alpha, A.matrix, x, y)
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

# ============================================================
# Low-level numeric refactorization kernels
# ============================================================
# These call the plain-C `_SparseRefactor*` symbols directly. They recompute the
# numeric factorization in place, reusing the symbolic factorization (the
# fill-reducing ordering and sparsity analysis) — the common case in Newton
# iterations, implicit time stepping, and parameter sweeps where a matrix's
# values change but its sparsity pattern does not.

const _REFACTOR_SYM = Dict(
    # (T, kind) => (symbol)
    (Cfloat,     :Symmetric) => :_SparseRefactorSymmetric_Float,
    (Cdouble,    :Symmetric) => :_SparseRefactorSymmetric_Double,
    (ComplexF32, :Symmetric) => :_SparseRefactorSymmetric_Complex_Float,
    (ComplexF64, :Symmetric) => :_SparseRefactorSymmetric_Complex_Double,
    # Complex-Hermitian Cholesky/LDLT refactor: libSparse's *Symmetric_Complex*
    # kernels reject Hermitian factorizations ("only applies to SparseSymmetric
    # matrices"), so a separate Hermitian family is required. (Real Hermitian ==
    # symmetric, so real types stay in the :Symmetric family above.)
    (ComplexF32, :Hermitian) => :_SparseRefactorHermitian_Complex_Float,
    (ComplexF64, :Hermitian) => :_SparseRefactorHermitian_Complex_Double,
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
        # Transient scratch for the refactor call (libSparse does not retain it).
        # A GC-managed buffer rooted across the ccall via `GC.@preserve`; Julia
        # array data is 16-byte aligned, matching the solver's expectation.
        workspace = Vector{UInt8}(undef, wsize)
        GC.@preserve workspace begin
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
                                  pointer(workspace)::Ptr{Cvoid})::Cvoid
        end
        return fact[]
    end
end

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
        # matrixObj.matrix holds raw pointers into the CSC buffers; keep the
        # owning wrapper alive across the factorization call.
        fact = GC.@preserve aa_fact SparseFactor(kind, aa_fact.matrixObj.matrix)
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
    _libsparse_throw(aa_fact._factorization.status, "solve")
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
    _libsparse_throw(aa_fact._factorization.status, "solve")
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
    _libsparse_throw(aa_fact._factorization.status, "solve")
    return x
end

# ============================================================
# refactor! (numeric refactor reusing a symbolic factorization)
# ============================================================

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

For LU/Cholesky/LDLᵀ factorizations you can equivalently use the
`LinearAlgebra`-style spellings `lu!(f, A)`, `cholesky!(f, A)`, and `ldlt!(f, A)`
(matching `SparseArrays`' symbolic-reuse API); they require `f` to already hold a
factorization of the matching kind and delegate here. QR reuse has no stdlib
spelling, so use `refactor!` directly for it.

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

    # The symbolic factorization's `attributes` field is zeroed, so the
    # Hermitian-vs-symmetric distinction (needed only for complex Cholesky/LDLT
    # refactor, which has dedicated Hermitian kernels) must come from the NEW
    # matrix `A`. Real Hermitian == symmetric, so real types stay :Symmetric.
    family = (T <: Complex && ishermitian(A)) ? :Hermitian :
             _refactor_family(aa_fact._factorization.symbolicFactorization.type)
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

# LinearAlgebra-style spellings of `refactor!`, matching how SparseArrays exposes
# symbolic-factorization reuse (`lu!(F, A)` for UMFPACK, `cholesky!(F, A)` /
# `ldlt!(F, A)` for CHOLMOD). They dispatch on our own `AAFactorization`, so they
# are not type piracy and do not touch Julia's `SparseMatrixCSC` factorizations.
# Each requires `F` to already hold a factorization of the matching kind, then
# delegates to `refactor!`. (QR reuse has no stdlib spelling — use `refactor!`.)
_is_lu_kind(t::SparseFactorization_t) = _refactor_family(t) == :LU
_is_cholesky_kind(t::SparseFactorization_t) = t == SparseFactorizationCholesky
_is_ldlt_kind(t::SparseFactorization_t) =
    t in (SparseFactorizationLDLT, SparseFactorizationLDLTUnpivoted,
          SparseFactorizationLDLTSBK, SparseFactorizationLDLTTPP)

function _assert_refactor_kind(aa_fact::AAFactorization, pred, name::AbstractString)
    f = aa_fact._factorization
    # If F isn't factored yet, let refactor! raise the clearer "call factor! first".
    f.status == SparseStatusOk || return
    pred(f.symbolicFactorization.type) || throw(ArgumentError(
        "$name requires F to hold the matching factorization kind; it holds " *
        "$(f.symbolicFactorization.type). Build F with that factorization first."))
    return
end

for (fn, pred) in ((:lu!, :_is_lu_kind),
                   (:cholesky!, :_is_cholesky_kind),
                   (:ldlt!, :_is_ldlt_kind))
    for AT in (:(AASparseMatrix{T}), :(SparseMatrixCSC{T, Int64}))
        @eval LinearAlgebra.$fn(aa_fact::AAFactorization{T}, A::$AT) where T<:vTypes =
            (_assert_refactor_kind(aa_fact, $pred, $(string(fn))); refactor!(aa_fact, A))
    end
end

# ============================================================
# AASparseMatrix -> SparseMatrixCSC round-trip
# ============================================================

"""
    SparseMatrixCSC(A::AASparseMatrix)

Materialize an Accelerate `AASparseMatrix` back into a standard Julia
`SparseMatrixCSC`. The transpose/adjoint/symmetric/Hermitian/triangular
attribute bits are honored, so the result equals the logical matrix `A`
represents (`size`, `getindex`).
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

    sym  = issymmetric(A)
    herm = T <: Complex && ishermitian(A)
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

# ============================================================
# Factorization entry points on AASparseMatrix
# ============================================================
# These dispatch on our own `AASparseMatrix`, NOT on `SparseMatrixCSC`, so they
# do not pirate or override Julia's UMFPACK/CHOLMOD sparse solvers. The explicit
# path for a Julia sparse matrix is e.g. `lu(AASparseMatrix(A))`.

"""
    lu(A::AASparseMatrix) -> AAFactorization

LU factorization of `A` via Apple's Sparse solvers (requires macOS 15.5+ for the
libSparse LU path). Dispatches on `AASparseMatrix`, so it does not override
Julia's `lu(::SparseMatrixCSC)`.
"""
LinearAlgebra.lu(A::AASparseMatrix) = (f = AAFactorization(A); factor!(f, SparseFactorizationLU); f)

"""
    cholesky(A::AASparseMatrix) -> AAFactorization

Cholesky factorization of `A` via Apple's Sparse solvers. `A` must be symmetric
(real) or Hermitian (complex) and positive definite. Dispatches on
`AASparseMatrix`, so it does not override Julia's `cholesky(::SparseMatrixCSC)`.
"""
LinearAlgebra.cholesky(A::AASparseMatrix) = (f = AAFactorization(A); factor!(f, SparseFactorizationCholesky); f)

"""
    qr(A::AASparseMatrix) -> AAFactorization

QR factorization of `A` via Apple's Sparse solvers. Works for rectangular
matrices and solves least-squares systems through `\\`. Dispatches on
`AASparseMatrix`, so it does not override Julia's `qr(::SparseMatrixCSC)`.
"""
LinearAlgebra.qr(A::AASparseMatrix) = (f = AAFactorization(A); factor!(f, SparseFactorizationQR); f)

"""
    ldlt(A::AASparseMatrix) -> AAFactorization

LDLᵀ factorization of `A` via Apple's Sparse solvers, for symmetric/Hermitian
(indefinite) matrices. Dispatches on `AASparseMatrix`, so it does not override
Julia's `ldlt(::SparseMatrixCSC)`.
"""
LinearAlgebra.ldlt(A::AASparseMatrix) = (f = AAFactorization(A); factor!(f, SparseFactorizationLDLT); f)

# ============================================================
# Iterative-solver option structs
# ============================================================
# Field-by-field mirrors of libSparse's option structs (verified against the
# raw layer: sizes 40/48/72, offsets exact). A zeroed field selects the library
# default (e.g. maxIterations 0 → 100, rtol 0 → sqrt(eps)).

"""Options for the conjugate-gradient iterative solver. Zeroed fields select
libSparse defaults (`maxIterations` 0 → 100, `rtol` 0 → `sqrt(eps)`)."""
struct SparseCGOptions
    reportError::Ptr{Cvoid}
    maxIterations::Cint
    atol::Cdouble
    rtol::Cdouble
    reportStatus::Ptr{Cvoid}
end
SparseCGOptions(; maxIterations::Integer = 0, atol::Real = 0.0, rtol::Real = 0.0) =
    SparseCGOptions(C_NULL, Cint(maxIterations), Float64(atol), Float64(rtol), C_NULL)

"""Options for the GMRES iterative solver. `variant` is `0`=DQGMRES (default),
`1`=GMRES, `2`=FGMRES; `nvec` is the number of orthogonalization vectors."""
struct SparseGMRESOptions
    reportError::Ptr{Cvoid}
    variant::UInt8
    nvec::Cint
    maxIterations::Cint
    atol::Cdouble
    rtol::Cdouble
    reportStatus::Ptr{Cvoid}
end
SparseGMRESOptions(; variant::Integer = 0, nvec::Integer = 0, maxIterations::Integer = 0,
                   atol::Real = 0.0, rtol::Real = 0.0) =
    SparseGMRESOptions(C_NULL, UInt8(variant), Cint(nvec), Cint(maxIterations),
                       Float64(atol), Float64(rtol), C_NULL)

"""Options for the LSMR least-squares iterative solver. `lambda` is the
(Tikhonov) damping factor; `nvec` the number of local-reorthogonalization
vectors; `convergenceTest` `0`=default, `1`=Fong-Saunders."""
struct SparseLSMROptions
    reportError::Ptr{Cvoid}
    lambda::Cdouble
    nvec::Cint
    convergenceTest::Cint
    atol::Cdouble
    rtol::Cdouble
    btol::Cdouble
    conditionLimit::Cdouble
    maxIterations::Cint
    reportStatus::Ptr{Cvoid}
end
SparseLSMROptions(; lambda::Real = 0.0, nvec::Integer = 0, convergenceTest::Integer = 0,
                  atol::Real = 0.0, rtol::Real = 0.0, btol::Real = 0.0,
                  conditionLimit::Real = 0.0, maxIterations::Integer = 0) =
    SparseLSMROptions(C_NULL, Float64(lambda), Cint(nvec), Cint(convergenceTest),
                      Float64(atol), Float64(rtol), Float64(btol),
                      Float64(conditionLimit), Cint(maxIterations), C_NULL)

# The 264-byte `SparseIterativeMethod`: a method tag (CG=0, GMRES=1, LSMR=2) at
# offset 0 and a 256-byte options union at offset 8. We pack the chosen option
# struct's bytes into the union payload.
struct SparseIterativeMethod
    method::Cint
    _reserved::Cint
    options::NTuple{256,UInt8}
end

function _pack_iter_options(opt)::NTuple{256,UInt8}
    buf = zeros(UInt8, 256)
    r = Ref(opt)
    GC.@preserve r buf begin
        unsafe_copyto!(pointer(buf),
                       Ptr{UInt8}(Base.unsafe_convert(Ptr{typeof(opt)}, r)),
                       sizeof(opt))
        return unsafe_load(Ptr{NTuple{256,UInt8}}(pointer(buf)))
    end
end

_iter_method(::Val{:cg},    o::SparseCGOptions)    = SparseIterativeMethod(Cint(0), Cint(0), _pack_iter_options(o))
_iter_method(::Val{:gmres}, o::SparseGMRESOptions) = SparseIterativeMethod(Cint(1), Cint(0), _pack_iter_options(o))
_iter_method(::Val{:lsmr},  o::SparseLSMROptions)  = SparseIterativeMethod(Cint(2), Cint(0), _pack_iter_options(o))

@enum SparseIterativeStatus_t::Int32 begin
    SparseIterativeConverged = 0
    SparseIterativeMaxIterations = 1
    SparseIterativeParameterError = -1
    SparseIterativeIllConditioned = -2
    SparseIterativeInternalError = -99
end

# ============================================================
# Preconditioners (opaque handle)
# ============================================================

const SparsePreconditioner_t = Cint
const SparsePreconditionerNone       = SparsePreconditioner_t(0)
const SparsePreconditionerUser       = SparsePreconditioner_t(1)
const SparsePreconditionerDiagonal   = SparsePreconditioner_t(2)
const SparsePreconditionerDiagScaling = SparsePreconditioner_t(3)

# Layout mirrors libSparse's SparseOpaquePreconditioner_* (type, mem, apply).
# All fields are pointers/ints, so the layout is identical for real and complex T
# (T is a phantom used only for dispatch).
struct SparseOpaquePreconditioner{T<:vTypes}
    type::SparsePreconditioner_t
    mem::Ptr{Cvoid}
    apply::Ptr{Cvoid}
end

"""Opaque preconditioner handle for the iterative solvers. Build with
[`AAPreconditioner`](@ref). Wraps a libSparse `SparseOpaquePreconditioner`; the
backing memory is released by a finalizer."""
mutable struct AAPreconditioner{T<:vTypes}
    _p::SparseOpaquePreconditioner{T}
    _matrix::AASparseMatrix{T}   # keep the source matrix rooted for the handle's life
end

const _PRECOND_CREATE = Dict(
    Cfloat     => :_SparseCreatePreconditioner_Float,
    Cdouble    => :_SparseCreatePreconditioner_Double,
    ComplexF32 => :_SparseCreatePreconditioner_Complex_Float,
    ComplexF64 => :_SparseCreatePreconditioner_Complex_Double,
)
const _PRECOND_RELEASE = Dict(
    Cfloat     => :_SparseReleaseOpaquePreconditioner_Float,
    Cdouble    => :_SparseReleaseOpaquePreconditioner_Double,
    ComplexF32 => :_SparseReleaseOpaquePreconditioner_Complex_Float,
    ComplexF64 => :_SparseReleaseOpaquePreconditioner_Complex_Double,
)

function _precond_symbol(p::Symbol)
    p === :diagonal    && return SparsePreconditionerDiagonal
    p === :diagscaling && return SparsePreconditionerDiagScaling
    throw(ArgumentError("unknown preconditioner $(repr(p)); use :diagonal or :diagscaling"))
end

for (T, create) in _PRECOND_CREATE
    release = _PRECOND_RELEASE[T]
    @eval function _create_preconditioner(type::SparsePreconditioner_t, A::AASparseMatrix{$T})
        # `_SparseCreatePreconditioner` reads A's structure/data (raw pointers into
        # A's CSC buffers), so keep A rooted across the call.
        p = GC.@preserve A begin
            @ccall LIBSPARSE.$create(type::SparsePreconditioner_t,
                A.matrix::Ref{SparseMatrix{$T}})::SparseOpaquePreconditioner{$T}
        end
        obj = AAPreconditioner{$T}(p, A)
        finalizer(obj) do o
            if o._p.type != SparsePreconditionerNone
                # Release the internally-allocated backing store. A stack `Ref`
                # copy is fine: the C call frees `->mem`, which is shared by the
                # copy. The finalizer runs once, so there is no double free.
                pref = Ref(o._p)
                GC.@preserve pref begin
                    @ccall LIBSPARSE.$release(
                        Base.unsafe_convert(Ptr{SparseOpaquePreconditioner{$T}}, pref)::Ptr{SparseOpaquePreconditioner{$T}})::Cvoid
                end
            end
        end
        return obj
    end
end

"""
    AAPreconditioner(A::AASparseMatrix; kind = :diagonal)

Construct a preconditioner for `A` to accelerate the iterative solvers
([`solve`](@ref) with a `method` keyword). `kind` is `:diagonal` (Jacobi) or
`:diagscaling`. `Float32`/`Float64`, and (macOS 15.5+) `ComplexF32`/`ComplexF64`.
"""
function AAPreconditioner(A::AASparseMatrix{T}; kind::Symbol = :diagonal) where {T<:vTypes}
    return _create_preconditioner(_precond_symbol(kind), A)
end

# ============================================================
# Iterative solve (CG / GMRES / LSMR)
# ============================================================
# The public iterative `SparseSolve(method, A, B, X[, precond])` entry points are
# exported (mangled) symbols; they build the operator-apply block internally, so
# we call them directly. A is passed by value (raw pointers into the CSC buffers
# → GC.@preserve the wrapper); B/X convert to DenseMatrix (rooting their arrays).

const _ITER_SOLVE_SYMS = Dict(
    # T => (base, enumPrecond, opaquePrecond)
    Cfloat => (
        :_Z11SparseSolve21SparseIterativeMethod18SparseMatrix_Float17DenseMatrix_FloatS1_,
        :_Z11SparseSolve21SparseIterativeMethod18SparseMatrix_Float17DenseMatrix_FloatS1_i,
        :_Z11SparseSolve21SparseIterativeMethod18SparseMatrix_Float17DenseMatrix_FloatS1_32SparseOpaquePreconditioner_Float),
    Cdouble => (
        :_Z11SparseSolve21SparseIterativeMethod19SparseMatrix_Double18DenseMatrix_DoubleS1_,
        :_Z11SparseSolve21SparseIterativeMethod19SparseMatrix_Double18DenseMatrix_DoubleS1_i,
        :_Z11SparseSolve21SparseIterativeMethod19SparseMatrix_Double18DenseMatrix_DoubleS1_33SparseOpaquePreconditioner_Double),
    # complex variants (macOS 15.5+). CG requires a Hermitian(-positive-definite)
    # matrix; GMRES/LSMR are general.
    ComplexF32 => (
        :_Z11SparseSolve21SparseIterativeMethod26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS1_,
        :_Z11SparseSolve21SparseIterativeMethod26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS1_i,
        :_Z11SparseSolve21SparseIterativeMethod26SparseMatrix_Complex_Float25DenseMatrix_Complex_FloatS1_40SparseOpaquePreconditioner_Complex_Float),
    ComplexF64 => (
        :_Z11SparseSolve21SparseIterativeMethod27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS1_,
        :_Z11SparseSolve21SparseIterativeMethod27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS1_i,
        :_Z11SparseSolve21SparseIterativeMethod27SparseMatrix_Complex_Double26DenseMatrix_Complex_DoubleS1_41SparseOpaquePreconditioner_Complex_Double),
)
for (T, (base, enumsym, opaquesym)) in _ITER_SOLVE_SYMS
    @eval _iter_solve!(m::SparseIterativeMethod, A::SparseMatrix{$T},
            B::StridedMatrix{$T}, X::StridedMatrix{$T}) =
        @ccall LIBSPARSE.$base(m::SparseIterativeMethod, A::SparseMatrix{$T},
            B::DenseMatrix{$T}, X::DenseMatrix{$T})::SparseIterativeStatus_t
    @eval _iter_solve!(m::SparseIterativeMethod, A::SparseMatrix{$T},
            B::StridedMatrix{$T}, X::StridedMatrix{$T}, pre::SparsePreconditioner_t) =
        @ccall LIBSPARSE.$enumsym(m::SparseIterativeMethod, A::SparseMatrix{$T},
            B::DenseMatrix{$T}, X::DenseMatrix{$T}, pre::SparsePreconditioner_t)::SparseIterativeStatus_t
    @eval _iter_solve!(m::SparseIterativeMethod, A::SparseMatrix{$T},
            B::StridedMatrix{$T}, X::StridedMatrix{$T}, pre::SparseOpaquePreconditioner{$T}) =
        @ccall LIBSPARSE.$opaquesym(m::SparseIterativeMethod, A::SparseMatrix{$T},
            B::DenseMatrix{$T}, X::DenseMatrix{$T}, pre::SparseOpaquePreconditioner{$T})::SparseIterativeStatus_t
end

function _check_iter_status(status::SparseIterativeStatus_t)
    status == SparseIterativeConverged && return
    status == SparseIterativeMaxIterations &&
        (@warn "iterative solver reached its iteration limit without converging"; return)
    status == SparseIterativeParameterError &&
        throw(ArgumentError("iterative solve failed: parameter error " *
            "(check the matrix/RHS dimensions and that CG's matrix is symmetric)"))
    status == SparseIterativeIllConditioned &&
        error("iterative solve failed: matrix is ill-conditioned")
    error("iterative solve failed: status = $status")
end

"""
    solve(A::AASparseMatrix, b; method, preconditioner = :none,
          atol = 0, rtol = 0, maxiter = 0, nvec = 0, variant = :dqgmres, lambda = 0)

Solve `A x = b` with an iterative Krylov method instead of a direct
factorization, returning `x`. `method` (required) is one of:

  - `:cg`    — conjugate gradient; `A` must be symmetric positive-definite.
  - `:gmres` — GMRES; for square non-symmetric / indefinite systems. `variant`
    is `:dqgmres` (default), `:gmres`, or `:fgmres`; `nvec` sets the number of
    orthogonalization vectors.
  - `:lsmr`  — LSMR least-squares; for rectangular or singular systems.
    `lambda` applies Tikhonov damping.

`atol`/`rtol` are the absolute/relative convergence tolerances (0 selects the
library default), `maxiter` the iteration cap (0 → 100). `preconditioner` may be
`:none`, `:diagonal`, `:diagscaling`, or an [`AAPreconditioner`](@ref). `b` may
be a vector or a matrix (multiple right-hand sides). `Float32`/`Float64`, and
(macOS 15.5+) `ComplexF32`/`ComplexF64` — for complex, `:cg` requires a Hermitian
positive-definite matrix.
"""
function solve(A::AASparseMatrix{T}, b::StridedVecOrMat{T};
               method::Symbol,
               preconditioner::Union{Symbol,AAPreconditioner} = :none,
               atol::Real = 0, rtol::Real = 0, maxiter::Integer = 0,
               nvec::Integer = 0, variant::Symbol = :dqgmres,
               lambda::Real = 0) where {T<:vTypes}
    m, n = size(A)
    size(b, 1) == m || throw(DimensionMismatch(
        "right-hand side has $(size(b, 1)) rows; A has $m rows"))
    if method === :cg
        meth = _iter_method(Val(:cg), SparseCGOptions(; maxIterations = maxiter, atol, rtol))
    elseif method === :gmres
        variantcode = variant === :dqgmres ? 0 : variant === :gmres ? 1 :
                      variant === :fgmres ? 2 :
                      throw(ArgumentError("unknown GMRES variant $(repr(variant))"))
        meth = _iter_method(Val(:gmres),
            SparseGMRESOptions(; variant = variantcode, nvec, maxIterations = maxiter, atol, rtol))
    elseif method === :lsmr
        meth = _iter_method(Val(:lsmr),
            SparseLSMROptions(; lambda, nvec, maxIterations = maxiter, atol, rtol))
    else
        throw(ArgumentError("unknown iterative method $(repr(method)); use :cg, :gmres, or :lsmr"))
    end
    nrhs = b isa AbstractVector ? 1 : size(b, 2)
    B = reshape(b, m, nrhs)
    X = zeros(T, n, nrhs)          # zero initial guess (also the output buffer)
    status = GC.@preserve A begin
        if preconditioner isa AAPreconditioner
            preconditioner isa AAPreconditioner{T} ||
                throw(ArgumentError("preconditioner eltype does not match A"))
            GC.@preserve preconditioner _iter_solve!(meth, A.matrix, B, X, preconditioner._p)
        elseif preconditioner === :none
            _iter_solve!(meth, A.matrix, B, X)
        else
            _iter_solve!(meth, A.matrix, B, X, _precond_symbol(preconditioner))
        end
    end
    _check_iter_status(status)
    return b isa AbstractVector ? vec(X) : X
end

solve(A::SparseMatrixCSC{T,Int64}, b::StridedVecOrMat{T}; kw...) where {T<:vTypes} =
    solve(AASparseMatrix(A), b; kw...)

# ============================================================
# Preallocated / thread-safe factorization solve workspace
# ============================================================
# The direct-solve C++ entry points have overloads taking a caller-owned
# `void* workspace`, letting repeated/concurrent solves avoid the per-call
# internal allocation. Size the buffer with `solve_workspace_size`.

const _SOLVE_WS_SYMS = Dict(
    # T => (outOfPlace(B,X,ws), inPlace(XB,ws))
    Cfloat  => (:_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseMatrix_FloatS0_Pv,
                :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseMatrix_FloatPv),
    Cdouble => (:_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_DoubleS0_Pv,
                :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_DoublePv),
    # complex variants (macOS 15.5+)
    ComplexF32 => (:_Z11SparseSolve39SparseOpaqueFactorization_Complex_Float25DenseMatrix_Complex_FloatS0_Pv,
                   :_Z11SparseSolve39SparseOpaqueFactorization_Complex_Float25DenseMatrix_Complex_FloatPv),
    ComplexF64 => (:_Z11SparseSolve40SparseOpaqueFactorization_Complex_Double26DenseMatrix_Complex_DoubleS0_Pv,
                   :_Z11SparseSolve40SparseOpaqueFactorization_Complex_Double26DenseMatrix_Complex_DoublePv),
)
for (T, (outsym, insym)) in _SOLVE_WS_SYMS
    @eval _solve_ws!(f::SparseOpaqueFactorization{$T}, B::StridedMatrix{$T},
            X::StridedMatrix{$T}, ws::Vector{UInt8}) =
        GC.@preserve ws (@ccall LIBSPARSE.$outsym(f::SparseOpaqueFactorization{$T},
            B::DenseMatrix{$T}, X::DenseMatrix{$T}, pointer(ws)::Ptr{Cvoid})::Cvoid)
    @eval _solve_ws!(f::SparseOpaqueFactorization{$T}, XB::StridedMatrix{$T},
            ws::Vector{UInt8}) =
        GC.@preserve ws (@ccall LIBSPARSE.$insym(f::SparseOpaqueFactorization{$T},
            XB::DenseMatrix{$T}, pointer(ws)::Ptr{Cvoid})::Cvoid)
end

"""
    solve_workspace_size(f::AAFactorization, nrhs = 1) -> Int

Number of bytes of scratch a workspace buffer must hold for a
[`solve!`](@ref)-with-workspace call on `f` with `nrhs` right-hand sides. The
factorization must already be computed (call [`factor!`](@ref) first). Use to
size the `workspace` argument of the preallocated-workspace `solve!`.
"""
function solve_workspace_size(f::AAFactorization, nrhs::Integer = 1)
    fac = f._factorization
    fac.status == SparseStatusOk || throw(ArgumentError(
        "solve_workspace_size requires a computed factorization; call factor! first"))
    return Int(fac.solveWorkspaceRequiredStatic) + Int(nrhs) * Int(fac.solveWorkspaceRequiredPerRHS)
end

"""
    solve!(f::AAFactorization, b, x, workspace::Vector{UInt8})

Solve `A x = b` writing the result into `x`, using the caller-supplied
`workspace` scratch buffer instead of allocating internally — for repeated or
concurrent solves that reuse a factorization without per-call allocation. Size
`workspace` with [`solve_workspace_size`](@ref)`(f, size(b, 2))`. Each
concurrent thread must use its **own** workspace and its own `x`. Returns `x`.

    solve!(f::AAFactorization, xb, workspace::Vector{UInt8})

In-place variant: `xb` holds the right-hand side on input and the solution on
output (square systems only).
"""
function solve!(aa_fact::AAFactorization{T}, b::StridedVecOrMat{T},
                x::StridedVecOrMat{T}, workspace::Vector{UInt8}) where {T<:vTypes}
    m, n = size(aa_fact.matrixObj)
    size(b, 1) == m || throw(DimensionMismatch("RHS has $(size(b,1)) rows; A has $m rows"))
    size(x, 1) == n || throw(DimensionMismatch("solution has $(size(x,1)) rows; A has $n cols"))
    nrhs = b isa AbstractVector ? 1 : size(b, 2)
    factor!(aa_fact)
    need = solve_workspace_size(aa_fact, nrhs)
    length(workspace) >= need || throw(ArgumentError(
        "workspace too small: need $need bytes, got $(length(workspace)); " *
        "size it with solve_workspace_size(f, nrhs)"))
    B = reshape(b, m, nrhs); X = reshape(x, n, nrhs)
    _solve_ws!(aa_fact._factorization, B, X, workspace)
    _libsparse_throw(aa_fact._factorization.status, "solve")
    return x
end

function solve!(aa_fact::AAFactorization{T}, xb::StridedVecOrMat{T},
                workspace::Vector{UInt8}) where {T<:vTypes}
    m, n = size(aa_fact.matrixObj)
    m == n || throw(ArgumentError("in-place workspace solve requires a square system"))
    size(xb, 1) == n || throw(DimensionMismatch("RHS has $(size(xb,1)) rows; A is $m×$n"))
    nrhs = xb isa AbstractVector ? 1 : size(xb, 2)
    factor!(aa_fact)
    need = solve_workspace_size(aa_fact, nrhs)
    length(workspace) >= need || throw(ArgumentError(
        "workspace too small: need $need bytes, got $(length(workspace))"))
    XB = reshape(xb, n, nrhs)
    _solve_ws!(aa_fact._factorization, XB, workspace)
    _libsparse_throw(aa_fact._factorization.status, "solve")
    return xb
end

# ============================================================
# Subfactor extraction & application (Q, R, L, D, P, ...)
# ============================================================

const SparseSubfactor_t = UInt8
const SparseSubfactorInvalid = SparseSubfactor_t(0)
const SparseSubfactorP    = SparseSubfactor_t(1)
const SparseSubfactorS    = SparseSubfactor_t(2)
const SparseSubfactorL    = SparseSubfactor_t(3)
const SparseSubfactorD    = SparseSubfactor_t(4)
const SparseSubfactorPLPS = SparseSubfactor_t(5)
const SparseSubfactorQ    = SparseSubfactor_t(6)
const SparseSubfactorR    = SparseSubfactor_t(7)
const SparseSubfactorRP   = SparseSubfactor_t(8)

# 128-byte layout: attributes@0, contents@4, factor@8, wsStatic@112, wsPerRHS@120
# (verified against the raw layer). Layout is identical for real and complex T
# (the embedded SparseOpaqueFactorization{T} has a phantom T, no T-typed field).
struct SparseOpaqueSubfactor{T<:vTypes}
    attributes::att_type
    contents::SparseSubfactor_t
    factor::SparseOpaqueFactorization{T}
    workspaceRequiredStatic::Csize_t
    workspaceRequiredPerRHS::Csize_t
end

"""Handle to a sub-factor (`Q`, `R`, `L`, `D`, `P`, …) of an
[`AAFactorization`](@ref), created with [`subfactor`](@ref). Holds a *borrowed*
reference into its parent factorization (kept alive by this object), so it needs
no separate release. Apply it with `*` (multiply) or `\\` (solve)."""
mutable struct AASubfactor{T<:vTypes}
    _sub::SparseOpaqueSubfactor{T}
    _parent::AAFactorization{T}
end

const _SUBF_WS_SYMS = Dict(
    Cfloat     => :_SparseGetWorkspaceRequired_Float,
    Cdouble    => :_SparseGetWorkspaceRequired_Double,
    ComplexF32 => :_SparseGetWorkspaceRequired_Complex_Float,
    ComplexF64 => :_SparseGetWorkspaceRequired_Complex_Double,
)
for (T, sym) in _SUBF_WS_SYMS
    @eval function _subfactor_ws(contents::SparseSubfactor_t, f::SparseOpaqueFactorization{$T})
        ws = Ref{Csize_t}(0); wp = Ref{Csize_t}(0)
        @ccall LIBSPARSE.$sym(contents::SparseSubfactor_t, f::SparseOpaqueFactorization{$T},
            ws::Ptr{Csize_t}, wp::Ptr{Csize_t})::Cvoid
        return (ws[], wp[])
    end
end

# Sub-factor view attributes, mirroring libSparse's SparseCreateSubfactor:
# L is lower-triangular, R/RP upper-triangular, everything else ordinary.
function _subfactor_attributes(c::SparseSubfactor_t)
    c == SparseSubfactorL && return ATT_TRI_LOWER
    (c == SparseSubfactorR || c == SparseSubfactorRP) && return ATT_TRI_UPPER
    return ATT_ORDINARY
end

"""
    subfactor(f::AAFactorization, which) -> AASubfactor

Extract an individual factor of the factorization `f` for direct application.
`which` is one of the `SparseSubfactor*` constants:
`SparseSubfactorQ`/`SparseSubfactorR` (from QR),
`SparseSubfactorL`/`SparseSubfactorD`/`SparseSubfactorP` (from Cholesky/LDLᵀ).
`f` is factored if necessary. Apply the result with `sub * x` (multiply by the
factor) or `sub \\ b` (solve against it). `Float32`/`Float64`, and (macOS 15.5+)
`ComplexF32`/`ComplexF64`.
"""
function subfactor(f::AAFactorization{T}, which::SparseSubfactor_t) where {T<:vTypes}
    factor!(f)
    fac = f._factorization
    fac.status == SparseStatusOk ||
        throw(ArgumentError("subfactor requires a completed factorization"))
    ws, wp = GC.@preserve f _subfactor_ws(which, fac)
    sub = SparseOpaqueSubfactor{T}(_subfactor_attributes(which), which, fac, ws, wp)
    return AASubfactor{T}(sub, f)
end

# The subfactor multiply/solve C++ entry points allocate their own scratch
# (sized from the workspace fields we filled in), so we call them directly.
const _SUBF_APPLY_SYMS = Dict(
    Cfloat  => (:_Z14SparseMultiply27SparseOpaqueSubfactor_Float17DenseMatrix_FloatS0_,
                :_Z11SparseSolve27SparseOpaqueSubfactor_Float17DenseMatrix_FloatS0_),
    Cdouble => (:_Z14SparseMultiply28SparseOpaqueSubfactor_Double18DenseMatrix_DoubleS0_,
                :_Z11SparseSolve28SparseOpaqueSubfactor_Double18DenseMatrix_DoubleS0_),
    # complex variants (macOS 15.5+)
    ComplexF32 => (:_Z14SparseMultiply35SparseOpaqueSubfactor_Complex_Float25DenseMatrix_Complex_FloatS0_,
                   :_Z11SparseSolve35SparseOpaqueSubfactor_Complex_Float25DenseMatrix_Complex_FloatS0_),
    ComplexF64 => (:_Z14SparseMultiply36SparseOpaqueSubfactor_Complex_Double26DenseMatrix_Complex_DoubleS0_,
                   :_Z11SparseSolve36SparseOpaqueSubfactor_Complex_Double26DenseMatrix_Complex_DoubleS0_),
)
for (T, (mulsym, solvesym)) in _SUBF_APPLY_SYMS
    @eval _subfactor_mul!(s::SparseOpaqueSubfactor{$T}, X::StridedMatrix{$T}, Y::StridedMatrix{$T}) =
        @ccall LIBSPARSE.$mulsym(s::SparseOpaqueSubfactor{$T},
            X::DenseMatrix{$T}, Y::DenseMatrix{$T})::Cvoid
    @eval _subfactor_solve!(s::SparseOpaqueSubfactor{$T}, B::StridedMatrix{$T}, X::StridedMatrix{$T}) =
        @ccall LIBSPARSE.$solvesym(s::SparseOpaqueSubfactor{$T},
            B::DenseMatrix{$T}, X::DenseMatrix{$T})::Cvoid
end

# Dimensions of the operand a sub-factor acts on, mirroring libSparse's
# `_SparseSubFactorGetDimn`: the parent is stored with rows ≥ cols; every
# sub-factor is n×n EXCEPT the Q of a QR (m×n); a transpose bit swaps them. As a
# map it takes an n-vector to an m-vector (`m×n` acting on the right).
function _subfactor_dimn(s::SparseOpaqueSubfactor)
    symb = s.factor.symbolicFactorization
    bs = max(Int(symb.blockSize), 1)
    m = Int(symb.rowCount) * bs
    n = Int(symb.columnCount) * bs
    m < n && ((m, n) = (n, m))    # parent always factored with m ≥ n
    isQ = symb.type == SparseFactorizationQR && s.contents == SparseSubfactorQ
    isQ || (m = n)                # all sub-factors are n×n except Q of QR
    (s.attributes & ATT_TRANSPOSE) != 0 && ((m, n) = (n, m))
    return (m, n)
end

"""
    sub * x

Multiply the vector/matrix `x` by the extracted sub-factor `sub` (e.g. form
`Q*x`). If `sub` acts as an `m×n` operator, `x` must have `n` rows and the
result has `m` rows.
"""
function Base.:(*)(sub::AASubfactor{T}, x::StridedVecOrMat{T}) where {T<:vTypes}
    m, n = _subfactor_dimn(sub._sub)
    size(x, 1) == n || throw(DimensionMismatch(
        "sub-factor multiply expects $(n) rows; got $(size(x, 1))"))
    nrhs = x isa AbstractVector ? 1 : size(x, 2)
    X = reshape(x, n, nrhs)
    Y = zeros(T, m, nrhs)
    GC.@preserve sub _subfactor_mul!(sub._sub, X, Y)
    return x isa AbstractVector ? vec(Y) : Y
end

"""
    sub \\ b

Solve a system against the extracted sub-factor `sub` (e.g. triangular solve
`R \\ b`, or apply `Qᵀ` via `Q \\ b`). If `sub` acts as an `m×n` operator, `b`
must have `m` rows and the result has `n` rows.
"""
function Base.:(\)(sub::AASubfactor{T}, b::StridedVecOrMat{T}) where {T<:vTypes}
    m, n = _subfactor_dimn(sub._sub)
    size(b, 1) == m || throw(DimensionMismatch(
        "sub-factor solve expects $(m) rows; got $(size(b, 1))"))
    nrhs = b isa AbstractVector ? 1 : size(b, 2)
    B = reshape(b, m, nrhs)
    X = zeros(T, n, nrhs)
    GC.@preserve sub _subfactor_solve!(sub._sub, B, X)
    return b isa AbstractVector ? vec(X) : X
end

# ============================================================
# Partial / low-rank LU refactorization update
# ============================================================
# Distinct from the full numeric `refactor!`: recompute only the L/U values that
# a from-scratch LU would change, given a small set of modified (row, col)
# entries. Requires a *pivotless* LU factorization (LUUnpivoted/LUSPP/LUTPP).

const _UPDATE_LU_SYMS = Dict(
    Cfloat     => :_SparseUpdatePartialRefactorLU_Float,
    Cdouble    => :_SparseUpdatePartialRefactorLU_Double,
    ComplexF32 => :_SparseUpdatePartialRefactorLU_Complex_Float,
    ComplexF64 => :_SparseUpdatePartialRefactorLU_Complex_Double,
)
for (T, sym) in _UPDATE_LU_SYMS
    @eval function _update_partial_lu!(fref::Base.RefValue{SparseOpaqueFactorization{$T}},
            updatedIndices::Vector{Cint}, newMatrix::SparseMatrix{$T})
        updateCount = Cint(length(updatedIndices) ÷ 2)
        GC.@preserve updatedIndices begin
            @ccall LIBSPARSE.$sym(
                Base.unsafe_convert(Ptr{SparseOpaqueFactorization{$T}}, fref)::Ptr{SparseOpaqueFactorization{$T}},
                updateCount::Cint, pointer(updatedIndices)::Ptr{Cint},
                newMatrix::SparseMatrix{$T})::Cvoid
        end
        return fref[]
    end
end

"""
    update_partial_lu!(f::AAFactorization, updated, A_new)

Apply a partial LU refactorization to `f` in place: recompute only the factor
values that a from-scratch LU of `A_new` would change, given that just the
entries listed in `updated` differ from the originally-factored matrix. `updated`
is a vector of 1-based `(row, col)` tuples of the modified positions; `A_new` is
a full copy of the matrix with those entries at their new values and the **same
sparsity pattern** as the original.

`f` must hold a **pivotless** LU factorization (`SparseFactorizationLUUnpivoted`,
`SparseFactorizationLUSPP`, or `SparseFactorizationLUTPP`) — build it with
`factor!(f, SparseFactorizationLUUnpivoted)`. Requires macOS 15.5+.
`Float32`/`Float64`/`ComplexF32`/`ComplexF64`. Returns `f`.

Distinct from [`refactor!`](@ref), which recomputes the entire numeric
factorization.
"""
function update_partial_lu!(f::AAFactorization{T},
        updated::AbstractVector{<:Tuple{Integer,Integer}},
        A_new::AASparseMatrix{T}) where {T<:vTypes}
    fac = f._factorization
    fac.status == SparseStatusOk || throw(ArgumentError(
        "update_partial_lu! requires a completed factorization; call factor! first"))
    t = fac.symbolicFactorization.type
    t in (SparseFactorizationLUUnpivoted, SparseFactorizationLUSPP, SparseFactorizationLUTPP) ||
        throw(ArgumentError("update_partial_lu! requires a pivotless LU factorization " *
            "(LUUnpivoted/LUSPP/LUTPP); f holds $t"))
    size(f.matrixObj) == size(A_new) || throw(DimensionMismatch(
        "update matrix size $(size(A_new)) ≠ factored size $(size(f.matrixObj))"))
    # Flatten to the 0-based (row, col) pair list libSparse expects.
    idx = Vector{Cint}(undef, 2 * length(updated))
    @inbounds for (k, (i, j)) in enumerate(updated)
        idx[2k-1] = Cint(i - 1)
        idx[2k]   = Cint(j - 1)
    end
    fref = Ref(fac)
    GC.@preserve A_new begin
        _update_partial_lu!(fref, idx, A_new.matrix)
    end
    newfac = fref[]
    _libsparse_throw(newfac.status, "update_partial_lu")
    f._factorization = newfac
    f.matrixObj = A_new
    return f
end

# ============================================================
# Factorization introspection (read options back)
# ============================================================

const _NUMOPTS_SYMS = Dict(
    Cfloat     => :_SparseGetOptionsFromNumericFactor_Float,
    Cdouble    => :_SparseGetOptionsFromNumericFactor_Double,
    ComplexF32 => :_SparseGetOptionsFromNumericFactor_Complex_Float,
    ComplexF64 => :_SparseGetOptionsFromNumericFactor_Complex_Double,
)
for (T, sym) in _NUMOPTS_SYMS
    @eval function _numeric_options(f::SparseOpaqueFactorization{$T})
        fref = Ref(f)
        GC.@preserve fref begin
            @ccall LIBSPARSE.$sym(
                Base.unsafe_convert(Ptr{SparseOpaqueFactorization{$T}}, fref)::Ptr{SparseOpaqueFactorization{$T}})::SparseNumericFactorOptions
        end
    end
end

"""
    numeric_options(f::AAFactorization) -> SparseNumericFactorOptions

Read back the numeric-factorization options (scaling method, pivot/zero
tolerances) that libSparse recorded for the completed factorization `f`. `f`
must already be factored. `Float32`/`Float64`/`ComplexF32`/`ComplexF64`.
"""
function numeric_options(f::AAFactorization{T}) where {T<:vTypes}
    f._factorization.status == SparseStatusOk || throw(ArgumentError(
        "numeric_options requires a completed factorization; call factor! first"))
    return _numeric_options(f._factorization)
end

"""
    symbolic_options(f::AAFactorization) -> SparseSymbolicFactorOptions

Read back the symbolic-factorization options (ordering method, etc.) recorded
for the completed factorization `f`.
"""
function symbolic_options(f::AAFactorization{T}) where {T<:vTypes}
    f._factorization.status == SparseStatusOk || throw(ArgumentError(
        "symbolic_options requires a completed factorization; call factor! first"))
    symb = Ref(f._factorization.symbolicFactorization)
    GC.@preserve symb begin
        return @ccall LIBSPARSE._SparseGetOptionsFromSymbolicFactor(
            Base.unsafe_convert(Ptr{SparseOpaqueSymbolicFactorization}, symb)::Ptr{SparseOpaqueSymbolicFactorization})::SparseSymbolicFactorOptions
    end
end
