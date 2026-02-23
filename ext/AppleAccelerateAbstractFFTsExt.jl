module AppleAccelerateAbstractFFTsExt

@static if Sys.isapple()

using AbstractFFTs
using AbstractFFTs: Plan, ScaledPlan
using LinearAlgebra
using AppleAccelerate
using AppleAccelerate: FFTSetup

# Forward FFT plan wrapping vDSP
mutable struct VDSPFFTPlan{T<:Union{Float32,Float64},N,G} <: Plan{Complex{T}}
    setup::FFTSetup{T}
    sz::NTuple{N,Int}
    region::G
    pinv::Plan{Complex{T}}
    function VDSPFFTPlan{T}(setup::FFTSetup{T}, sz::NTuple{N,Int}, region::G) where {T,N,G}
        new{T,N,G}(setup, sz, region)
    end
end

# Backward (unnormalized inverse) FFT plan wrapping vDSP
mutable struct VDSPBFFTPlan{T<:Union{Float32,Float64},N,G} <: Plan{Complex{T}}
    setup::FFTSetup{T}
    sz::NTuple{N,Int}
    region::G
    pinv::Plan{Complex{T}}
    function VDSPBFFTPlan{T}(setup::FFTSetup{T}, sz::NTuple{N,Int}, region::G) where {T,N,G}
        new{T,N,G}(setup, sz, region)
    end
end

Base.size(p::VDSPFFTPlan) = p.sz
Base.size(p::VDSPBFFTPlan) = p.sz

AbstractFFTs.AdjointStyle(::VDSPFFTPlan) = AbstractFFTs.FFTAdjointStyle()
AbstractFFTs.AdjointStyle(::VDSPBFFTPlan) = AbstractFFTs.FFTAdjointStyle()

function _check_vdsp_fft(x::AbstractVector)
    n = length(x)
    n > 0 || throw(ArgumentError("input array must be non-empty"))
    ispow2(n) || throw(ArgumentError("vDSP FFT requires power-of-2 length, got $n"))
    return nothing
end

function _check_vdsp_fft(x::AbstractMatrix)
    nrows, ncols = size(x)
    nrows > 0 && ncols > 0 || throw(ArgumentError("input matrix must be non-empty"))
    ispow2(nrows) || throw(ArgumentError("vDSP 2D FFT requires power-of-2 row count, got $nrows"))
    ispow2(ncols) || throw(ArgumentError("vDSP 2D FFT requires power-of-2 column count, got $ncols"))
    return nothing
end

# 1D plans
function AbstractFFTs.plan_fft(x::Vector{Complex{T}}, region=1:1; kwargs...) where {T<:Union{Float32,Float64}}
    _check_vdsp_fft(x)
    setup = FFTSetup{T}(length(x))
    VDSPFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_bfft(x::Vector{Complex{T}}, region=1:1; kwargs...) where {T<:Union{Float32,Float64}}
    _check_vdsp_fft(x)
    setup = FFTSetup{T}(length(x))
    VDSPBFFTPlan{T}(setup, size(x), region)
end

# 2D plans
function AbstractFFTs.plan_fft(x::Matrix{Complex{T}}, region=1:2; kwargs...) where {T<:Union{Float32,Float64}}
    _check_vdsp_fft(x)
    nrows, ncols = size(x)
    setup = FFTSetup{T}(max(nrows, ncols))
    VDSPFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_bfft(x::Matrix{Complex{T}}, region=1:2; kwargs...) where {T<:Union{Float32,Float64}}
    _check_vdsp_fft(x)
    nrows, ncols = size(x)
    setup = FFTSetup{T}(max(nrows, ncols))
    VDSPBFFTPlan{T}(setup, size(x), region)
end

# 1D execution
function Base.:*(p::VDSPFFTPlan{T,1}, x::AbstractVector) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.fft(convert(Vector{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPBFFTPlan{T,1}, x::AbstractVector) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.bfft(convert(Vector{Complex{T}}, x), p.setup)
end

# 2D execution
function Base.:*(p::VDSPFFTPlan{T,2}, x::AbstractMatrix) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.fft(convert(Matrix{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPBFFTPlan{T,2}, x::AbstractMatrix) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.bfft(convert(Matrix{Complex{T}}, x), p.setup)
end

# 1D mul!
function LinearAlgebra.mul!(y::AbstractVector{Complex{T}}, p::VDSPFFTPlan{T,1}, x::AbstractVector{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.fft(convert(Vector{Complex{T}}, x), p.setup))
end

function LinearAlgebra.mul!(y::AbstractVector{Complex{T}}, p::VDSPBFFTPlan{T,1}, x::AbstractVector{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.bfft(convert(Vector{Complex{T}}, x), p.setup))
end

# 2D mul!
function LinearAlgebra.mul!(y::AbstractMatrix{Complex{T}}, p::VDSPFFTPlan{T,2}, x::AbstractMatrix{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.fft(convert(Matrix{Complex{T}}, x), p.setup))
end

function LinearAlgebra.mul!(y::AbstractMatrix{Complex{T}}, p::VDSPBFFTPlan{T,2}, x::AbstractMatrix{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.bfft(convert(Matrix{Complex{T}}, x), p.setup))
end

function AbstractFFTs.plan_inv(p::VDSPFFTPlan{T}) where {T}
    bplan = VDSPBFFTPlan{T}(p.setup, p.sz, p.region)
    N = AbstractFFTs.normalization(Complex{T}, p.sz, p.region)
    bplan.pinv = ScaledPlan(p, N)
    ScaledPlan(bplan, N)
end

function AbstractFFTs.plan_inv(p::VDSPBFFTPlan{T}) where {T}
    fplan = VDSPFFTPlan{T}(p.setup, p.sz, p.region)
    N = AbstractFFTs.normalization(Complex{T}, p.sz, p.region)
    fplan.pinv = ScaledPlan(p, N)
    ScaledPlan(fplan, N)
end

end # @static if Sys.isapple()

end # module
