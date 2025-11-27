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

for T in (Cfloat, Cdouble)
    # write a macro to make this less copy-paste heavy? Most of the function
    # bodies are straight @ccalls, replacing Strided with Dense, but a few aren't.

    # also may want to add guardrails in case the mangled function names
    # change at some later date.

    # see https://discourse.julialang.org/t/apple-accelerate-sparse-solvers/110175
    local dmMultMangled = T == Cfloat ? :_Z14SparseMultiply18SparseMatrix_Float17DenseMatrix_FloatS0_ :
                        :_Z14SparseMultiply19SparseMatrix_Double18DenseMatrix_DoubleS0_
    @eval begin
        function SparseMultiply(arg1::SparseMatrix{$T}, arg2::StridedMatrix{$T},
                                                            arg3::StridedMatrix{$T})
            @ccall LIBSPARSE.$dmMultMangled(arg1::SparseMatrix{$T}, arg2::DenseMatrix{$T},
                                                arg3::DenseMatrix{$T})::Cvoid
        end
    end

    local dvMultMangled = T == Cfloat ? :_Z14SparseMultiply18SparseMatrix_Float17DenseVector_FloatS0_ :
                        :_Z14SparseMultiply19SparseMatrix_Double18DenseVector_DoubleS0_
    @eval begin
        function SparseMultiply(arg1::SparseMatrix{$T}, arg2::StridedVector{$T},
                                                        arg3::StridedVector{$T})
            @ccall LIBSPARSE.$dvMultMangled(arg1::SparseMatrix{$T}, arg2::DenseVector{$T},
                                            arg3::DenseVector{$T})::Cvoid
        end
    end

    local sdmMultMangled = T == Cfloat ? :_Z14SparseMultiplyf18SparseMatrix_Float17DenseMatrix_FloatS0_ :
                                :_Z14SparseMultiplyd19SparseMatrix_Double18DenseMatrix_DoubleS0_
    @eval begin
        function SparseMultiply(arg1::$T, arg2::SparseMatrix{$T}, arg3::StridedMatrix{$T},
                                                    arg4::StridedMatrix{$T})
            @ccall LIBSPARSE.$sdmMultMangled(arg1::$T, arg2::SparseMatrix{$T},
                                                arg3::DenseMatrix{$T}, arg4::DenseMatrix{$T})::Cvoid
        end
    end

    local sdvMultMangled = T == Cfloat ? :_Z14SparseMultiplyf18SparseMatrix_Float17DenseVector_FloatS0_ :
                                        :_Z14SparseMultiplyd19SparseMatrix_Double18DenseVector_DoubleS0_
    @eval SparseMultiply(arg1::$T,
                        arg2::SparseMatrix{$T},
                        arg3::StridedVector{$T},
                        arg4::StridedVector{$T}) = @ccall (
            LIBSPARSE.$sdvMultMangled(arg1::$T,
                                    arg2::SparseMatrix{$T},
                                    arg3::DenseVector{$T},
                                    arg4::DenseVector{$T})::Cvoid
        )

    local dmMultAddMangled = T == Cfloat ? :_Z17SparseMultiplyAdd18SparseMatrix_Float17DenseMatrix_FloatS0_ :
                                            :_Z17SparseMultiplyAdd19SparseMatrix_Double18DenseMatrix_DoubleS0_
    @eval SparseMultiplyAdd(arg1::SparseMatrix{$T},
                            arg2::StridedMatrix{$T},
                            arg3::StridedMatrix{$T}) = @ccall (
            LIBSPARSE.$dmMultAddMangled(arg1::SparseMatrix{$T},
                                        arg2::DenseMatrix{$T},
                                        arg3::DenseMatrix{$T})::Cvoid
        )

    local dvMultAddMangled = T == Cfloat ? :_Z17SparseMultiplyAdd18SparseMatrix_Float17DenseVector_FloatS0_ :
                                            :_Z17SparseMultiplyAdd19SparseMatrix_Double18DenseVector_DoubleS0_
    @eval SparseMultiplyAdd(arg1::SparseMatrix{$T},
                            arg2::StridedVector{$T},
                            arg3::StridedVector{$T}) = @ccall(
            LIBSPARSE.$dvMultAddMangled(arg1::SparseMatrix{$T},
                                        arg2::DenseVector{$T},
                                        arg3::DenseVector{$T})::Cvoid
        )

    local sdmMultAddMangled = T == Cfloat ? :_Z17SparseMultiplyAddf18SparseMatrix_Float17DenseMatrix_FloatS0_ :
                                        :_Z17SparseMultiplyAddd19SparseMatrix_Double18DenseMatrix_DoubleS0_
    @eval SparseMultiplyAdd(arg0::$T,
                            arg1::SparseMatrix{$T},
                            arg2::StridedMatrix{$T},
                            arg3::StridedMatrix{$T}) = @ccall (
            LIBSPARSE.$sdmMultAddMangled(arg0::$T,
                                        arg1::SparseMatrix{$T},
                                        arg2::DenseMatrix{$T},
                                        arg3::DenseMatrix{$T})::Cvoid
        )

    local sdvMultAddMangled = T == Cfloat ? :_Z17SparseMultiplyAddf18SparseMatrix_Float17DenseVector_FloatS0_ :
                                        :_Z17SparseMultiplyAddd19SparseMatrix_Double18DenseVector_DoubleS0_
    @eval SparseMultiplyAdd(arg0::$T,
                            arg1::SparseMatrix{$T},
                            arg2::StridedVector{$T},
                            arg3::StridedVector{$T}) = @ccall (
            LIBSPARSE.$sdvMultAddMangled(arg0::$T,
                                        arg1::SparseMatrix{$T},
                                        arg2::DenseVector{$T},
                                        arg3::DenseVector{$T})::Cvoid
        )

    local mTransposeMangled = T == Cfloat ? :_Z18SparseGetTranspose18SparseMatrix_Float :
                                            :_Z18SparseGetTranspose19SparseMatrix_Double
    @eval SparseGetTranspose(arg1::SparseMatrix{$T}) = @ccall(
        LIBSPARSE.$mTransposeMangled(arg1::SparseMatrix{$T})::SparseMatrix{$T}
    )
    # skipped: 2 subfactor transposes.
    local ofTransposeMangled = T == Cfloat ? :_Z18SparseGetTranspose31SparseOpaqueFactorization_Float :
                                        :_Z18SparseGetTranspose32SparseOpaqueFactorization_Double
    @eval SparseGetTranspose(arg1::SparseOpaqueFactorization{$T}) = @ccall (
            LIBSPARSE.$ofTransposeMangled(arg1::SparseOpaqueFactorization{$T})::SparseOpaqueFactorization{$T}
    )

    # TODO: these SparseConvertFromCoord functions are unused and untested.
    local convertCoordMangled = T == Cfloat ? :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKf :
                                            :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKd
    @eval begin
        function SparseConvertFromCoord(arg1::Cint, arg2::Cint, arg3::Clong, arg4::Cuchar, arg5::$att_type,
                                        arg6::Ptr{Cint}, arg7::Ptr{Cint}, arg8::Ptr{$T})
            @ccall LIBSPARSE.$convertCoordMangled(arg1::Cint, arg2::Cint, arg3::Clong, arg4::Cuchar, arg5::$att_type,
                                                    arg6::Ptr{Cint}, arg7::Ptr{Cint}, arg8::Ptr{$T})::SparseMatrix{$T}
        end
    end

    local uConvertCoordMangled = T == Cfloat ? :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKdPvS4_ :
                                            :_Z27SparseConvertFromCoordinateiilh18SparseAttributes_tPKiS1_PKfPvS4_
    @eval begin
        function SparseConvertFromCoord(arg1::Cint, arg2::Cint, arg3::Clong, arg4::Cuchar, arg5::$att_type,
                    arg6::Ptr{Cint}, arg7::Ptr{Cint}, arg8::Ptr{$T}, arg9::Ptr{Cvoid}, arg10::Ptr{Cvoid})
            @ccall LIBSPARSE.$uConvertCoordMangled(arg1::Cint, arg2::Cint, arg3::Clong, arg4::Cuchar, arg5::$att_type,
                    arg6::Ptr{Cint}, arg7::Ptr{Cint}, arg8::Ptr{$T}, arg9::Ptr{Cvoid}, arg10::Ptr{Cvoid})::SparseMatrix{$T}
        end
    end

    local mCleanup = T == Cfloat ? :_Z13SparseCleanup18SparseMatrix_Float :
                                     :_Z13SparseCleanup19SparseMatrix_Double
    @eval SparseCleanup(arg1::SparseMatrix{$T}) = @ccall (
            LIBSPARSE.$mCleanup(arg1::SparseMatrix{$T})::Cvoid
    )
    
    local ofCleanup = T == Cfloat ? :_Z13SparseCleanup31SparseOpaqueFactorization_Float :
                                    :_Z13SparseCleanup32SparseOpaqueFactorization_Double
    @eval SparseCleanup(arg1::SparseOpaqueFactorization{$T}) = @ccall (
            LIBSPARSE.$ofCleanup(arg1::SparseOpaqueFactorization{$T})::Cvoid
    )

    @eval SparseFactor(arg1::SparseFactorization_t,
                        arg2::SparseMatrix{$T},
                        noErrors::Bool = false) = 
                noErrors ? SparseFactor(arg1, arg2) :
                SparseFactor(arg1, arg2, SparseSymbolicFactorOptions(), SparseNumericFactorOptions($T)) 

    local sparseFactorMatrix = T == Cfloat ? :_Z12SparseFactorh18SparseMatrix_Float :
                                                :_Z12SparseFactorh19SparseMatrix_Double
    @eval function SparseFactorNoErrors(arg1::SparseFactorization_t,
                                arg2::SparseMatrix{$T})::SparseOpaqueFactorization{$T}
        arg1 != SparseFactorizationTBD || throw(ArgumentError("Factorization type must be specified"))
        @ccall(LIBSPARSE.$sparseFactorMatrix(arg1::Cuint,
                                        arg2::SparseMatrix{$T})::SparseOpaqueFactorization{$T})
    end

    local sparseFactorMatrixOpts = T == Cfloat ? :_Z12SparseFactorh18SparseMatrix_Float27SparseSymbolicFactorOptions26SparseNumericFactorOptions :
                                           :_Z12SparseFactorh19SparseMatrix_Double27SparseSymbolicFactorOptions26SparseNumericFactorOptions     
    @eval function SparseFactor(arg1::SparseFactorization_t,
                    arg2::SparseMatrix{$T},
                    arg3::SparseSymbolicFactorOptions,
                    arg4::SparseNumericFactorOptions)
        arg1 != SparseFactorizationTBD || throw(ArgumentError("Factorization type must be specified"))
        @ccall(LIBSPARSE.$sparseFactorMatrixOpts(
                    arg1::Cuint,
                    arg2::SparseMatrix{$T},
                    arg3::SparseSymbolicFactorOptions,
                    arg4::SparseNumericFactorOptions
            )::SparseOpaqueFactorization{$T}
        )
    end

    local sparseSolveInplace = T == Cfloat ? :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseMatrix_Float :
                                            :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_Double
    @eval function SparseSolve(arg1::SparseOpaqueFactorization{$T},
                            arg2::StridedMatrix{$T})
        @ccall LIBSPARSE.$sparseSolveInplace(arg1::SparseOpaqueFactorization{$T},
                                                arg2::DenseMatrix{$T})::Cvoid
    end

    local sparseSolve = T == Cfloat ? :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseMatrix_FloatS0_ :
                                :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_DoubleS0_
    @eval SparseSolve(arg1::SparseOpaqueFactorization{$T}, arg2::StridedMatrix{$T},
                        arg3::StridedMatrix{$T}) = @ccall (
        LIBSPARSE.$sparseSolve(arg1::SparseOpaqueFactorization{$T}, arg2::DenseMatrix{$T}, 
                                    arg3::DenseMatrix{$T})::Cvoid
    )

    local sparseSolveVecInPlace = T == Cfloat ? :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseVector_Float :
                                :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseVector_Double
    @eval function SparseSolve(arg1::SparseOpaqueFactorization{$T},
                                arg2::StridedVector{$T})
        @ccall LIBSPARSE.$sparseSolveVecInPlace(arg1::SparseOpaqueFactorization{$T},
                                            arg2::DenseVector{$T})::Cvoid
        resize!(arg2, arg1.symbolicFactorization.columnCount)
    end

    local sparseSolveVec = T == Cfloat ? :_Z11SparseSolve31SparseOpaqueFactorization_Float17DenseVector_FloatS0_ :
                                    :_Z11SparseSolve32SparseOpaqueFactorization_Double18DenseVector_DoubleS0_
    @eval SparseSolve(arg1::SparseOpaqueFactorization{$T},
                    arg2::StridedVector{$T},
                    arg3::StridedVector{$T}) = @ccall (
        LIBSPARSE.$sparseSolveVec(arg1::SparseOpaqueFactorization{$T},
                                arg2::DenseVector{$T},
                                arg3::DenseVector{$T})::Cvoid
    )
end

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