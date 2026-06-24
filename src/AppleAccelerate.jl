module AppleAccelerate
using Libdl, LinearAlgebra

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
const libacc_info_plist = "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Resources/Info.plist"

# Cached macOS version, populated once in __init__()
const _macos_version = Ref{Union{Nothing,VersionNumber}}(nothing)

# VecLib Threading API: /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/vecLib.framework/Headers/thread_api.h
@enum Threading::Cuint begin
    BLAS_THREADING_MULTI_THREADED
    BLAS_THREADING_SINGLE_THREADED
end

function forward_accelerate(interface::Symbol;
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

"""
    load_accelerate(; clear = false, verbose = false, load_ilp64 = true)

Load Accelerate, replacing the current LBT forwarding tables if `clear` is `true`. `clear`
is `false` by default to allow for OpenBLAS to act as a fallback for operations missing
from Accelerate, such as `gemmt`. Attempts to load the ILP64 symbols if `load_ilp64` is
`true`, and errors out if unable.
"""
function load_accelerate(; clear::Bool = false,
                           verbose::Bool = false,
                           load_ilp64::Bool = true)
    libacc_hdl = dlopen_e(libacc)
    if libacc_hdl == C_NULL
        return
    end

    # Check to see if we can load ILP64 symbols
    if load_ilp64 && dlsym_e(libacc_hdl, "dgemm\$NEWLAPACK\$ILP64") == C_NULL
        error("Unable to load ILP64 interface from '$(libacc)'; you are running macOS $(get_macos_version()), you need v13.4+")
    end

    # First, load :lp64 symbols, optionally clearing the current LBT forwarding tables
    forward_accelerate(:lp64; new_lapack=true, clear, verbose)
    if load_ilp64
        forward_accelerate(:ilp64; new_lapack=true, verbose)
    end
end

"""
    _read_macos_version()

Read the macOS version from the system plist and return it as a `VersionNumber`.
Early macOS Tahoe betas reported version 16.x before Apple finalized the version
numbering as 26.x; this function normalizes 16.x → 26.x for correct comparisons.
Returns `nothing` on non-Apple platforms or if the version cannot be determined.
"""
function _read_macos_version()
    @static if !Sys.isapple()
        return nothing
    end

    plist_lines = split(String(read("/System/Library/CoreServices/SystemVersion.plist")), "\n")
    vers_idx = findfirst(l -> occursin("ProductVersion", l), plist_lines)
    if vers_idx === nothing
        return nothing
    end

    m = match(r">([\d\.]+)<", plist_lines[vers_idx+1])
    if m === nothing
        return nothing
    end

    cap = only(m.captures)
    cap === nothing && return nothing
    ver = VersionNumber(cap)
    # Early macOS Tahoe developer betas reported 16.x before Apple settled on 26.x
    if ver.major == 16
        return VersionNumber(26, ver.minor, ver.patch)
    end
    return ver
end

"""
    get_macos_version()

Return the current macOS version as a `VersionNumber`.
Returns `nothing` on non-Apple platforms or if the version cannot be determined.
"""
function get_macos_version()
    ver = _macos_version[]
    if ver === nothing
        # Fallback: read directly if called before __init__ (e.g., during precompilation)
        return _read_macos_version()
    end
    return ver
end

"""
    set_num_threads(n::Integer) -> BlasInt

Set the number of threads used by Accelerate BLAS. If `n == 1`, use single-threaded mode;
if `n > 1`, use multi-threaded mode. Returns the resulting thread count (from
[`get_num_threads`](@ref)). On macOS < 15 where the threading API is unavailable,
warns and returns `1`.
"""
function set_num_threads(n::Integer)
    ver = get_macos_version()
    if ver === nothing || ver < v"15"
        @warn "The Accelerate threading API requires macOS 15 or later; ignoring" maxlog=1
        return LinearAlgebra.BlasInt(1)
    end

    retval::Cint = -1
    if n == 1
        retval = ccall((:BLASSetThreading, libacc), Cint, (Cint,), BLAS_THREADING_SINGLE_THREADED)
    elseif n > 1
        retval = ccall((:BLASSetThreading, libacc), Cint, (Cint,), BLAS_THREADING_MULTI_THREADED)
    end
    @assert retval == 0 "AppleAccelerate: Call to BLASSetThreading failed"
    return get_num_threads()
end

"""
    get_num_threads() -> BlasInt

Return the number of threads used by Accelerate BLAS. Returns `1` for single-threaded mode,
or the actual thread count for multi-threaded mode. On macOS < 15 where the threading API
is unavailable, warns and returns `1`.
"""
function get_num_threads()::LinearAlgebra.BlasInt
    ver = get_macos_version()
    if ver === nothing || ver < v"15"
        @warn "The Accelerate threading API requires macOS 15 or later" maxlog=1
        return LinearAlgebra.BlasInt(1)
    end

    retval::Threading = ccall((:BLASGetThreading, libacc), Threading, ())
    if retval == BLAS_THREADING_SINGLE_THREADED
        return LinearAlgebra.BlasInt(1)
    elseif retval == BLAS_THREADING_MULTI_THREADED
        return ccall((:APPLE_NTHREADS, libacc), LinearAlgebra.BlasInt, ())
    else
        error("AppleAccelerate: BLASGetThreading returned unexpected value: $(retval)")
    end
end

function __init__()
    @static !Sys.isapple() && return
    _macos_version[] = _read_macos_version()
    ver = _macos_version[]
    # dsptrf has a bug in the initial release of the $NEWLAPACK symbols in 13.3.
    # Thus use macOS 13.4 for ILP64, a correct LAPACK, and threading APIs
    if ver === nothing || ver < v"13.4"
        @info "AppleAccelerate.jl needs macOS 13.4 or later for BLAS/LAPACK forwarding"
        return
    end
    load_accelerate(; clear = false, load_ilp64=true)
    # libSparse lives at a hard-coded path (LIBSPARSE in sparse.jl). Probe it once
    # here so a future macOS layout change surfaces a clear diagnostic instead of
    # an opaque dlopen error at the first sparse ccall.
    if dlopen_e(LIBSPARSE) == C_NULL
        @warn "AppleAccelerate.jl: libSparse could not be loaded from '$(LIBSPARSE)'; \
               sparse linear algebra (AASparseMatrix, factorizations, solvers) will be unavailable"
    end
end

@static if Sys.isapple()
    include("array.jl")
    include("complexarray.jl")
    include("dsp.jl")
    include("sparse.jl")
end

end # module
