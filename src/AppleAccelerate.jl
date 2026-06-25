module AppleAccelerate
using Libdl

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
const libacc_info_plist = "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Resources/Info.plist"

# Cached macOS version, populated once in __init__()
const _macos_version = Ref{Union{Nothing,VersionNumber}}(nothing)

# VecLib Threading API: /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/vecLib.framework/Headers/thread_api.h
@enum Threading::Cuint begin
    BLAS_THREADING_MULTI_THREADED
    BLAS_THREADING_SINGLE_THREADED
end

# BLAS/LAPACK forwarding to Accelerate is provided by the LinearAlgebra package
# extension (`ext/AppleAccelerateLinearAlgebraExt.jl`), since it depends on
# `LinearAlgebra.BLAS.lbt_forward`. These are declared here as generic-function
# stubs so the names resolve from core (e.g. for `__init__`/docs), with methods
# added by the extension once `LinearAlgebra` is loaded.
"""
    forward_accelerate(interface::Symbol; new_lapack, clear, verbose)

Forward BLAS/LAPACK symbols to Accelerate via libblastrampoline.
Defined by the LinearAlgebra package extension; requires `using LinearAlgebra`.
"""
function forward_accelerate end

"""
    load_accelerate(; clear = false, verbose = false, load_ilp64 = true)

Load Accelerate, replacing the current LBT forwarding tables if `clear` is `true`. `clear`
is `false` by default to allow for OpenBLAS to act as a fallback for operations missing
from Accelerate, such as `gemmt`. Attempts to load the ILP64 symbols if `load_ilp64` is
`true`, and errors out if unable.

Defined by the LinearAlgebra package extension; requires `using LinearAlgebra`.
With only `using AppleAccelerate`, BLAS/LAPACK forwarding does NOT happen.
"""
function load_accelerate end

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
    set_num_threads(n::Integer) -> Int

Set the number of threads used by Accelerate BLAS. If `n == 1`, use single-threaded mode;
if `n > 1`, use multi-threaded mode. Returns the resulting thread count (from
[`get_num_threads`](@ref)). On macOS < 15 where the threading API is unavailable,
warns and returns `1`.
"""
function set_num_threads(n::Integer)
    n < 1 && throw(ArgumentError("number of threads must be ≥ 1"))
    ver = get_macos_version()
    if ver === nothing || ver < v"15"
        @warn "The Accelerate threading API requires macOS 15 or later; ignoring" maxlog=1
        return Int(1)
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
    get_num_threads() -> Int

Return the number of threads used by Accelerate BLAS. Returns `1` for single-threaded mode,
or the actual thread count for multi-threaded mode. On macOS < 15 where the threading API
is unavailable, warns and returns `1`.
"""
function get_num_threads()::Int
    ver = get_macos_version()
    if ver === nothing || ver < v"15"
        @warn "The Accelerate threading API requires macOS 15 or later" maxlog=1
        return Int(1)
    end

    retval::Threading = ccall((:BLASGetThreading, libacc), Threading, ())
    if retval == BLAS_THREADING_SINGLE_THREADED
        return Int(1)
    elseif retval == BLAS_THREADING_MULTI_THREADED
        # APPLE_NTHREADS returns a C `int` (BlasInt would also work, but core
        # must not depend on LinearAlgebra); read it as Cint and widen to Int.
        return Int(ccall((:APPLE_NTHREADS, libacc), Cint, ()))
    else
        error("AppleAccelerate: BLASGetThreading returned unexpected value: $(retval)")
    end
end

function __init__()
    @static !Sys.isapple() && return
    _macos_version[] = _read_macos_version()
    ver = _macos_version[]
    # dsptrf has a bug in the initial release of the $NEWLAPACK symbols in 13.3.
    # Thus use macOS 13.4 for ILP64, a correct LAPACK, and threading APIs.
    #
    # NOTE: BLAS/LAPACK forwarding is NO LONGER performed here. `using
    # AppleAccelerate` alone does not touch LBT. Forwarding happens in the
    # LinearAlgebra package extension's `__init__` (triggered by
    # `using LinearAlgebra`), which re-applies the same macOS 13.4+ guard.
    if ver === nothing || ver < v"13.4"
        @info "AppleAccelerate.jl needs macOS 13.4 or later for BLAS/LAPACK forwarding " *
              "(load LinearAlgebra to enable it)"
    end
    # libSparse lives at a hard-coded path (LIBSPARSE in sparse.jl). Probe it once
    # here so a future macOS layout change surfaces a clear diagnostic instead of
    # an opaque dlopen error at the first sparse ccall.
    if dlopen_e(LIBSPARSE) == C_NULL
        @warn "AppleAccelerate.jl: libSparse could not be loaded from '$(LIBSPARSE)'; \
               sparse linear algebra (AASparseMatrix, factorizations, solvers) will be unavailable"
    end
end

@static if Sys.isapple()
    # Raw, auto-generated ABI layer (Clang.jl). Do not edit by hand — regenerate
    # with `julia --project=gen gen/generate.jl`. The idiomatic wrappers below
    # call into this submodule instead of embedding ccall strings.
    include("lib/LibAccelerate.jl")

    include("array.jl")
    include("vmath.jl")
    include("complexarray.jl")
    include("dsp.jl")
    include("sparse.jl")
    include("quadrature.jl")
    include("bnns.jl")
end

end # module
