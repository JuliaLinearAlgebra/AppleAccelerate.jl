using SparseArrays
using Libdl

# header found at: /Library/Developer/CommandLineTools/SDKs/MacOSX15.sdk/System/
# Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/
# Versions/A/Headers/Sparse/Solve.h (replace MacOSX15 with correct version number)

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
    SparseFactorizationTBD = 64 # my own addition.
end

@enum SparseOrder_t::UInt8 begin
    SparseOrderDefault = 0
    SparseOrderUser = 1
    SparseOrderAMD = 2
    SparseOrderMetis = 3
    SparseOrderCOLAMD = 4
end

@enum SparseScaling_t::UInt8 begin
    SpraseScalingDefault = 0
    SparseScalingUser = 1
    SparseScalingEquilibriationInf = 2
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
SparseFactorNoErrors(arg1::SparseFactorization_t, arg2::SparseMatrixStructure) = @ccall (
    LIBSPARSE._Z12SparseFactorh21SparseMatrixStructure(
        arg1::Cuint, arg2::SparseMatrixStructure
    )::SparseOpaqueSymbolicFactorization
)

SparseFactor(arg1::SparseFactorization_t,
            arg2::SparseMatrixStructure,
            arg3::SparseSymbolicFactorOptions) = @ccall(
     LIBSPARSE._Z12SparseFactorh21SparseMatrixStructure27SparseSymbolicFactorOptions(
                arg1::Cuint,
                arg2::SparseMatrixStructure,
                arg3::SparseSymbolicFactorOptions
    )::SparseOpaqueSymbolicFactorization
)