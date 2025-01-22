module AppleAccelerate

using LinearAlgebra, Libdl
#using LAPACK_jll, LAPACK32_jll # Needed if use_external_lapack == true

# For now, only use BLAS from Accelerate (that is to say, vecLib)
global const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
global const libacc_info_plist = "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Resources/Info.plist"

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
    load_accelerate(;clear = true, verbose = false)

Load Accelerate, replacing the current LBT forwarding tables if `clear` is `true`.
Attempts to load the ILP64 symbols if `load_ilp64` is `true`, and errors out if unable.
"""
function load_accelerate(;clear::Bool = true, verbose::Bool = false, load_ilp64::Bool = true, use_external_lapack::Bool = true)
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
        error("Unable to load ILP64 interface from '$(libacc)'; You are running macOS version $(get_macos_version()), you need v13.4+")
    end

    # First, load :lp64 symbols, optionally clearing the current LBT forwarding tables
    forward_accelerate(:lp64; new_lapack=true, clear, verbose)
    if load_ilp64
        forward_accelerate(:ilp64; new_lapack=true, verbose)
    end

    # Next, load an external LAPACK, if requested
    if use_external_lapack
        if load_ilp64
            BLAS.lbt_forward(LAPACK_jll.liblapack_path; suffix_hint="64_", verbose)
        end
        BLAS.lbt_forward(LAPACK32_jll.liblapack32_path; verbose)
    end
end

function get_macos_version()
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

    return VersionNumber(only(m.captures))
end

function __init__()
    Sys.isapple() || return
    ver = get_macos_version()
    # Default to loading the ILP64 interface on macOS 13.3+
    # dsptrf has a bug in the initial release of the $NEWLAPACK symbols in 13.3.
    # Thus use macOS 13.4 for ILP64 and a correct LAPACK
    ver < v"13.4" && return
    load_accelerate(; load_ilp64=true, use_external_lapack=false)
end

if Sys.isapple()
    include("Util.jl")
    include("Array.jl")
    include("DSP.jl")
    include("libSparse/wrappers.jl")
    include("libSparse/AASparseMatrices.jl")
    include("libSparse/AAFactorizations.jl")
    export AASparseMatrix, muladd!
    export AAFactorization, solve, solve!, factor!
end

end # module
