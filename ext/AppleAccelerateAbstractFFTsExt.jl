module AppleAccelerateAbstractFFTsExt

@static if Sys.isapple()

using AbstractFFTs
using AbstractFFTs: Plan, ScaledPlan
using LinearAlgebra
using AppleAccelerate
using AppleAccelerate: FFTSetup

# --- Validity checks ---
# Used in plan creation to decide whether to use vDSP or fall through.
_can_vdsp(x::AbstractVector) = length(x) > 0 && ispow2(length(x))
function _can_vdsp(x::AbstractMatrix)
    nrows, ncols = size(x)
    nrows > 0 && ncols > 0 && ispow2(nrows) && ispow2(ncols)
end
_can_vdsp(x::AbstractArray) = false  # 3D+ not supported

# Used in plan execution to guard against misuse of an existing plan.
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

# Error with a clear message when vDSP can't handle the input.
function _vdsp_unsupported(x)
    if ndims(x) > 2
        throw(ArgumentError("vDSP only supports 1D and 2D FFTs (got $(ndims(x))D array). Use FFTW.jl for higher-dimensional transforms."))
    elseif length(x) == 0
        throw(ArgumentError("vDSP FFT requires non-empty arrays."))
    else
        throw(ArgumentError("vDSP FFT requires power-of-2 dimensions (got size $(size(x))). Use FFTW.jl for arbitrary sizes."))
    end
end

# === Out-of-place plan types ===

mutable struct VDSPFFTPlan{T<:Union{Float32,Float64},N,G} <: Plan{Complex{T}}
    setup::FFTSetup{T}
    sz::NTuple{N,Int}
    region::G
    pinv::Plan{Complex{T}}
    function VDSPFFTPlan{T}(setup::FFTSetup{T}, sz::NTuple{N,Int}, region::G) where {T,N,G}
        new{T,N,G}(setup, sz, region)
    end
end

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

# --- Out-of-place plan creation (1D, 2D, N-d) ---

function AbstractFFTs.plan_fft(x::Vector{Complex{T}}, region=1:1; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    setup = FFTSetup{T}(length(x))
    VDSPFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_bfft(x::Vector{Complex{T}}, region=1:1; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    setup = FFTSetup{T}(length(x))
    VDSPBFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_fft(x::Matrix{Complex{T}}, region=1:2; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    nrows, ncols = size(x)
    setup = FFTSetup{T}(max(nrows, ncols))
    VDSPFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_bfft(x::Matrix{Complex{T}}, region=1:2; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    nrows, ncols = size(x)
    setup = FFTSetup{T}(max(nrows, ncols))
    VDSPBFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_fft(x::Array{Complex{T},N}, region=1:N; kwargs...) where {T<:Union{Float32,Float64},N}
    _vdsp_unsupported(x)
end

function AbstractFFTs.plan_bfft(x::Array{Complex{T},N}, region=1:N; kwargs...) where {T<:Union{Float32,Float64},N}
    _vdsp_unsupported(x)
end

# --- Out-of-place execution ---

function Base.:*(p::VDSPFFTPlan{T,1}, x::AbstractVector) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.fft(convert(Vector{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPBFFTPlan{T,1}, x::AbstractVector) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.bfft(convert(Vector{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPFFTPlan{T,2}, x::AbstractMatrix) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.fft(convert(Matrix{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPBFFTPlan{T,2}, x::AbstractMatrix) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.bfft(convert(Matrix{Complex{T}}, x), p.setup)
end

# --- Out-of-place mul! ---

function LinearAlgebra.mul!(y::AbstractVector{Complex{T}}, p::VDSPFFTPlan{T,1}, x::AbstractVector{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.fft(convert(Vector{Complex{T}}, x), p.setup))
end

function LinearAlgebra.mul!(y::AbstractVector{Complex{T}}, p::VDSPBFFTPlan{T,1}, x::AbstractVector{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.bfft(convert(Vector{Complex{T}}, x), p.setup))
end

function LinearAlgebra.mul!(y::AbstractMatrix{Complex{T}}, p::VDSPFFTPlan{T,2}, x::AbstractMatrix{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.fft(convert(Matrix{Complex{T}}, x), p.setup))
end

function LinearAlgebra.mul!(y::AbstractMatrix{Complex{T}}, p::VDSPBFFTPlan{T,2}, x::AbstractMatrix{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, AppleAccelerate.bfft(convert(Matrix{Complex{T}}, x), p.setup))
end

# --- Out-of-place plan_inv ---

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

# === In-place plan types ===

mutable struct VDSPInplaceFFTPlan{T<:Union{Float32,Float64},N,G} <: Plan{Complex{T}}
    setup::FFTSetup{T}
    sz::NTuple{N,Int}
    region::G
    pinv::Plan{Complex{T}}
    function VDSPInplaceFFTPlan{T}(setup::FFTSetup{T}, sz::NTuple{N,Int}, region::G) where {T,N,G}
        new{T,N,G}(setup, sz, region)
    end
end

mutable struct VDSPInplaceBFFTPlan{T<:Union{Float32,Float64},N,G} <: Plan{Complex{T}}
    setup::FFTSetup{T}
    sz::NTuple{N,Int}
    region::G
    pinv::Plan{Complex{T}}
    function VDSPInplaceBFFTPlan{T}(setup::FFTSetup{T}, sz::NTuple{N,Int}, region::G) where {T,N,G}
        new{T,N,G}(setup, sz, region)
    end
end

Base.size(p::VDSPInplaceFFTPlan) = p.sz
Base.size(p::VDSPInplaceBFFTPlan) = p.sz

AbstractFFTs.AdjointStyle(::VDSPInplaceFFTPlan) = AbstractFFTs.FFTAdjointStyle()
AbstractFFTs.AdjointStyle(::VDSPInplaceBFFTPlan) = AbstractFFTs.FFTAdjointStyle()

# --- In-place plan creation (1D, 2D, N-d) ---

function AbstractFFTs.plan_fft!(x::Vector{Complex{T}}, region=1:1; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    setup = FFTSetup{T}(length(x))
    VDSPInplaceFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_bfft!(x::Vector{Complex{T}}, region=1:1; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    setup = FFTSetup{T}(length(x))
    VDSPInplaceBFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_fft!(x::Matrix{Complex{T}}, region=1:2; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    nrows, ncols = size(x)
    setup = FFTSetup{T}(max(nrows, ncols))
    VDSPInplaceFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_bfft!(x::Matrix{Complex{T}}, region=1:2; kwargs...) where {T<:Union{Float32,Float64}}
    _can_vdsp(x) || _vdsp_unsupported(x)
    nrows, ncols = size(x)
    setup = FFTSetup{T}(max(nrows, ncols))
    VDSPInplaceBFFTPlan{T}(setup, size(x), region)
end

function AbstractFFTs.plan_fft!(x::Array{Complex{T},N}, region=1:N; kwargs...) where {T<:Union{Float32,Float64},N}
    _vdsp_unsupported(x)
end

function AbstractFFTs.plan_bfft!(x::Array{Complex{T},N}, region=1:N; kwargs...) where {T<:Union{Float32,Float64},N}
    _vdsp_unsupported(x)
end

# --- In-place execution ---

function Base.:*(p::VDSPInplaceFFTPlan{T,1}, x::AbstractVector) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.fft!(convert(Vector{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPInplaceBFFTPlan{T,1}, x::AbstractVector) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.bfft!(convert(Vector{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPInplaceFFTPlan{T,2}, x::AbstractMatrix) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.fft!(convert(Matrix{Complex{T}}, x), p.setup)
end

function Base.:*(p::VDSPInplaceBFFTPlan{T,2}, x::AbstractMatrix) where {T}
    _check_vdsp_fft(x)
    AppleAccelerate.bfft!(convert(Matrix{Complex{T}}, x), p.setup)
end

# --- In-place mul! ---

function LinearAlgebra.mul!(y::AbstractVector{Complex{T}}, p::VDSPInplaceFFTPlan{T,1}, x::AbstractVector{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, x)
    AppleAccelerate.fft!(convert(Vector{Complex{T}}, y), p.setup)
end

function LinearAlgebra.mul!(y::AbstractVector{Complex{T}}, p::VDSPInplaceBFFTPlan{T,1}, x::AbstractVector{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, x)
    AppleAccelerate.bfft!(convert(Vector{Complex{T}}, y), p.setup)
end

function LinearAlgebra.mul!(y::AbstractMatrix{Complex{T}}, p::VDSPInplaceFFTPlan{T,2}, x::AbstractMatrix{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, x)
    AppleAccelerate.fft!(convert(Matrix{Complex{T}}, y), p.setup)
end

function LinearAlgebra.mul!(y::AbstractMatrix{Complex{T}}, p::VDSPInplaceBFFTPlan{T,2}, x::AbstractMatrix{Complex{T}}) where {T}
    _check_vdsp_fft(x)
    copyto!(y, x)
    AppleAccelerate.bfft!(convert(Matrix{Complex{T}}, y), p.setup)
end

# --- In-place plan_inv ---

function AbstractFFTs.plan_inv(p::VDSPInplaceFFTPlan{T}) where {T}
    bplan = VDSPInplaceBFFTPlan{T}(p.setup, p.sz, p.region)
    N = AbstractFFTs.normalization(Complex{T}, p.sz, p.region)
    bplan.pinv = ScaledPlan(p, N)
    ScaledPlan(bplan, N)
end

function AbstractFFTs.plan_inv(p::VDSPInplaceBFFTPlan{T}) where {T}
    fplan = VDSPInplaceFFTPlan{T}(p.setup, p.sz, p.region)
    N = AbstractFFTs.normalization(Complex{T}, p.sz, p.region)
    fplan.pinv = ScaledPlan(p, N)
    ScaledPlan(fplan, N)
end

end # @static if Sys.isapple()

end # module
