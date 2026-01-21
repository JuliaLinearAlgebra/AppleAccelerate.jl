module AppleAccelerate
using LinearAlgebra, Libdl

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
const libacc_info_plist = "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Resources/Info.plist"

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
    # Silently exit on non-Accelerate-capable platforms
    @static if !Sys.isapple()
        return
    end
    libacc_hdl = dlopen_e(libacc)
    if libacc_hdl == C_NULL
        return
    end

    # Check to see if we can load ILP64 symbols
    if load_ilp64 && dlsym_e(libacc_hdl, "dgemm\$NEWLAPACK\$ILP64") == C_NULL
        @error "Unable to load ILP64 interface from '$(libacc)'; You are running macOS version $(get_macos_version()), you need v13.4+"
    end

    # First, load :lp64 symbols, optionally clearing the current LBT forwarding tables
    forward_accelerate(:lp64; new_lapack=true, clear, verbose)
    if load_ilp64
        forward_accelerate(:ilp64; new_lapack=true, verbose)
    end
end

function get_macos_version(normalize=true)
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

    ver = VersionNumber(only(m.captures))
    if normalize && ver.major == 16
        return VersionNumber(26, ver.minor, ver.patch)
    end
    return ver
end

function set_num_threads(n::LinearAlgebra.BlasInt)
    @static if Sys.isapple()
        if get_macos_version() < v"15"
            @warn "The threading API is only available in macOS 15 and later"
            return -1
        end
    else
        return -1
    end

    retval::Cint = -1
    if n == 1
        retval = ccall((:BLASSetThreading, libacc), Cint, (Cint,), BLAS_THREADING_SINGLE_THREADED)
    elseif n > 1
        retval = ccall((:BLASSetThreading, libacc), Cint, (Cint,), BLAS_THREADING_MULTI_THREADED)
    end
    @assert retval == 0 "AppleAccelerate: Call to BlasSetThreading failed"
    return nothing
end

function get_num_threads()::LinearAlgebra.BlasInt
    @static if Sys.isapple()
        if get_macos_version() < v"15"
            @warn "The threading API is only available in macOS 15 and later"
            return -1
        end
    else
        return -1
    end

    retval::Threading = ccall((:BLASGetThreading, libacc), Threading, ())
    if retval == BLAS_THREADING_SINGLE_THREADED
        return LinearAlgebra.BlasInt(1)
    elseif retval == BLAS_THREADING_MULTI_THREADED
        return ccall((:APPLE_NTHREADS, libacc), LinearAlgebra.BlasInt, ())
    else
        @error "AppleAccelerate: Call to BlasGetThreading failed"
    end
end

function __init__()
    Sys.isapple() || return
    ver = get_macos_version()
    # Default to loading the ILP64 interface on macOS 13.3+
    # dsptrf has a bug in the initial release of the $NEWLAPACK symbols in 13.3.
    # Thus use macOS 13.4 for ILP64, a correct LAPACK, and threading APIs
    if ver < v"13.4"
        @info "AppleAccelerate.jl needs macOS 13.4 or later"
        return
    end
    load_accelerate(; clear = false, load_ilp64=true)
end

if Sys.isapple()
    include("Util.jl")
    include("Array.jl")
    include("DSP.jl")
    include("libSparse/wrappers.jl")
    include("libSparse/AASparseMatrices.jl")
    include("libSparse/AAFactorizations.jl")
end

end # module
