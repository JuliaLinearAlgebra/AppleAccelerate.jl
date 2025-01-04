using BFloat16s

include("libBNNS.jl")

bnnsdatatype_modifier(::Type{T}) where {T <: Union{AbstractFloat, Bool}} = BNNSDataTypeFloatBit
bnnsdatatype_modifier(::Type{T}) where {T <: Signed} = BNNSDataTypeIntBit
bnnsdatatype_modifier(::Type{T}) where {T <: Unsigned} = BNNSDataTypeUIntBit
bnnsdatatype_modifier(::Type{Bool}) = BNNSDataTypeMiscellaneousBit
bnnsdatatype_modifier(::Type{BFloat16}) = 0x18000

Base.convert(::Type{BNNSDataType}, T) = BNNSDataType(bnnsdatatype_modifier(T) | UInt32(sizeof(T)*8))

function BNNSNDArrayDescriptor(arr::AbstractArray{T, N}) where {T,N}
    N > 8 && throw(ArgumentError("BNNSNDArrays do not support more than 8 dimensions."))


    layout = BNNSDataLayout(UInt32(N) * UInt32(BNNSDataLayoutVector) | 0x8000)
    # layout = datalayout[N]
    sz = ntuple(Val(8)) do i
        Csize_t(get(size(arr), i, 0))
    end
    stride = ntuple(_ -> Csize_t(0), Val(8))
    return GC.@preserve arr BNNSNDArrayDescriptor(BNNSNDArrayFlagBackpropSet,
                                 layout,
                                 sz,
                                 stride,
                                 Ptr{Nothing}(pointer(arr)),
                                 T,
                                 0,
                                 T,
                                 1,
                                 0)
end

# Definitions for the Random extension
function bnns_rng end
function default_rng end
function rand end
function randn end
function rand! end
function randn! end
function seed! end
