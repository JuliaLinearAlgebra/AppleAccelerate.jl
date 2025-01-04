module RandomExt

@static if Sys.isapple()

using BFloat16s
using AppleAccelerate: BNNS
using .BNNS: BNNSFilterParameters,
             BNNSRandomGeneratorMethodAES_CTR,
             BNNSCreateRandomGenerator,
             BNNSCreateRandomGeneratorWithSeed,
             BNNSRandomGeneratorStateSize,
             BNNSRandomGeneratorSetState,
             BNNSRandomGeneratorGetState,
             BNNSNDArrayDescriptor,
             BNNSRandomFillNormalFloat,
             BNNSRandomFillUniformFloat,
             BNNSRandomFillUniformInt
using Random: Random, AbstractRNG

"""
    RNG()

A random number generator using AppleAccelerate's BNNS functionality.
"""
mutable struct RNG <: AbstractRNG
    ptr::Ptr{Nothing}
    function RNG(filter_parameters::Union{Nothing, BNNSFilterParameters}=nothing)
        params = isnothing(filter_parameters) ? Ptr{BNNSFilterParameters}(0) : [filter_parameters]
        res = new(BNNSCreateRandomGenerator(BNNSRandomGeneratorMethodAES_CTR, params))
        # finalizer(res) do
        #     BNNSDestroyRandomGenerator(res.ptr)
        # end
        return res
    end
    function RNG(seed::Integer, filter_parameters::Union{Nothing, BNNSFilterParameters}=nothing)
        seed = seed%UInt64
        params = isnothing(filter_parameters) ? Ptr{BNNSFilterParameters}(0) : [filter_parameters]
        res = new(BNNSCreateRandomGeneratorWithSeed(BNNSRandomGeneratorMethodAES_CTR, seed, params))
        # finalizer(res) do
        #     BNNSDestroyRandomGenerator(res.ptr)
        # end
        return res
    end
end

BNNS.bnns_rng() = RNG()
BNNS.bnns_rng(seed::Integer) = RNG(seed)

@static if isdefined(Base, :Memory) #VERSION >= v"1.11"
    function _get_rng_state(rng::RNG)
        stateSize = BNNSRandomGeneratorStateSize(rng.ptr)
        state = Memory{UInt8}(undef, Int64(stateSize))
        BNNSRandomGeneratorGetState(rng.ptr, stateSize, state)
        return state
    end
else
    function _get_rng_state(rng::RNG)
        stateSize = BNNSRandomGeneratorStateSize(rng.ptr)
        state = Vector{UInt8}(undef, Int64(stateSize))
        BNNSRandomGeneratorGetState(rng.ptr, stateSize, state)
        return state
    end
end

function Base.copy!(dest::RNG, src::RNG)
    state = _get_rng_state(src)
    BNNSRandomGeneratorSetState(dest.ptr, length(state), state)
    return dest
end

function Base.copy(rng::RNG)
    newrng = RNG()
    return copy!(newrng, rng)
end

Base.:(==)(rng1::RNG, rng2::RNG) = _get_rng_state(rng1) == _get_rng_state(rng2)

function Random.seed!(rng::RNG, seed::Integer)
    return copy!(rng, RNG(seed))
end

function Random.seed!(rng::RNG)
    return copy!(rng, RNG())
end

const GLOBAL_RNG = Ref{RNG}()
function BNNS.default_rng()
    if !isassigned(GLOBAL_RNG)
        GLOBAL_RNG[] = BNNS.bnns_rng()
    end
    return GLOBAL_RNG[]
end

const BNNSInt = Union{UInt8,Int8,UInt16,Int16,UInt32,Int32,UInt64,Int64}
const BNNSFloat = Union{Float16, Float32, BFloat16}

const BNNSUniform = Union{<:BNNSInt,<:BNNSFloat}
const BNNSNormal = BNNSFloat

function Random.rand!(rng::RNG, A::DenseArray{T}) where {T<:BNNSInt}
    isempty(A) && return A
    desc = Ref(BNNSNDArrayDescriptor(A))
    res = BNNSRandomFillUniformInt(rng.ptr, desc, typemin(signed(T)), typemax(signed(T)))
    @assert res == 0
    return A
end
function Random.rand!(rng::RNG, A::DenseArray{T}) where {T<:BNNSFloat}
    isempty(A) && return A
    desc = Ref(BNNSNDArrayDescriptor(A))
    res = BNNSRandomFillUniformFloat(rng.ptr, desc, T(0), T(1))
    @assert res == 0
    return A
end
function Random.randn!(rng::RNG, A::DenseArray{T}) where {T<:BNNSFloat}
    isempty(A) && return A
    desc = Ref(BNNSNDArrayDescriptor(A))
    res = BNNSRandomFillNormalFloat(rng.ptr, desc, Float32(0), Float32(1))
    @assert res == 0
    return A
end

# Out of place
Random.rand(rng::RNG, ::Type{T}, dims::Dims) where T <: BNNSUniform =
    Random.rand!(rng, Array{T,length(dims)}(undef, dims...))
Random.randn(rng::RNG, ::Type{T}, dims::Dims) where T <: BNNSNormal =
    Random.randn!(rng, Array{T,length(dims)}(undef, dims...))

# support all dimension specifications
Random.rand(rng::RNG, ::Type{T}, dim1::Integer, dims::Integer...) where T <: BNNSUniform =
    Random.rand!(rng, Array{T,length(dims) + 1}(undef, dim1, dims...))
Random.randn(rng::RNG, ::Type{T}, dim1::Integer, dims::Integer...) where T <: BNNSNormal =
    Random.randn!(rng, Array{T,length(dims) + 1}(undef, dim1, dims...))

# untyped out-of-place
Random.rand(rng::RNG, dim1::Integer, dims::Integer...) =
    Random.rand!(rng, Array{Float32,length(dims) + 1}(undef, dim1, dims...))
Random.randn(rng::RNG, dim1::Integer, dims::Integer...) =
    Random.randn!(rng, Array{Float32,length(dims) + 1}(undef, dim1, dims...))

# scalars
Random.rand(rng::RNG, T::Union{Type{Float16}, Type{Float32}, Type{BFloat16},
Type{Int8}, Type{UInt8},
Type{Int16}, Type{UInt16},
Type{Int32}, Type{UInt32},
Type{Int64}, Type{UInt64}}=Float32) = Random.rand(rng, T, 1)[1]

# This is the only way I could fix method ambiguity
Random.randn(rng::RNG, T::Type{BFloat16}) = Random.randn(rng, T, 1)[1]
Random.randn(rng::RNG, T::Type{Float16}) = Random.randn(rng, T, 1)[1]
Random.randn(rng::RNG, T::Type{Float32}) = Random.randn(rng, T, 1)[1]
Random.randn(rng::RNG) = Random.randn(rng, Float32)


# GPUArrays out-of-place
function BNNS.rand(::Type{T}, dims::Dims) where T <: BNNSUniform
    return Random.rand!(BNNS.default_rng(), Array{T,length(dims)}(undef, dims...))
end
function BNNS.randn(::Type{T}, dims::Dims) where T <: BNNSNormal
    return Random.randn!(BNNS.default_rng(), Array{T,length(dims)}(undef, dims...))
end

# support all dimension specifications
function BNNS.rand(::Type{T}, dim1::Integer, dims::Integer...) where T <: BNNSUniform
    return Random.rand!(BNNS.default_rng(), Array{T,length(dims) + 1}(undef, dim1, dims...))
end
function BNNS.randn(::Type{T}, dim1::Integer, dims::Integer...) where T <: BNNSNormal
    return Random.randn!(BNNS.default_rng(), Array{T,length(dims) + 1}(undef, dim1, dims...))
end

# untyped out-of-place
BNNS.rand(dim1::Integer, dims::Integer...) =
    Random.rand!(BNNS.default_rng(), Array{Float32,length(dims) + 1}(undef, dim1, dims...))
BNNS.randn(dim1::Integer, dims::Integer...) =
    Random.randn!(BNNS.default_rng(), Array{Float32,length(dims) + 1}(undef, dim1, dims...))

# scalars
BNNS.rand(T::Type=Float32) = BNNS.rand(T, 1)[1]
BNNS.randn(T::Type=Float32) = BNNS.randn(T, 1)[1]

# seeding
function BNNS.seed!(seed=Base.rand(UInt64))
    Random.seed!(BNNS.default_rng(), seed)
end

end
end # module
