module LibAccelerate

# Accelerate ships as a macOS *system framework*. Unlike the usual Clang.jl + JLL
# flow, there is no artifact to load — every Mac already has it. We ccall it by
# absolute path, mirroring the `libacc` const in the top-level AppleAccelerate module.
const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"


@enum quadrature_status::Int32 begin
    QUADRATURE_SUCCESS = 0
    QUADRATURE_ERROR = -1
    QUADRATURE_INVALID_ARG_ERROR = -2
    QUADRATURE_ALLOC_ERROR = -3
    QUADRATURE_INTERNAL_ERROR = -99
    QUADRATURE_INTEGRATE_MAX_EVAL_ERROR = -101
    QUADRATURE_INTEGRATE_BAD_BEHAVIOUR_ERROR = -102
end

# typedef void ( * quadrature_function_array ) ( void * _Null_unspecified __arg , size_t __n , const double * __x , double * __y )
const quadrature_function_array = Ptr{Cvoid}

@enum quadrature_integrator::UInt32 begin
    QUADRATURE_INTEGRATE_QNG = 0
    QUADRATURE_INTEGRATE_QAG = 1
    QUADRATURE_INTEGRATE_QAGS = 2
end

struct quadrature_integrate_function
    fun::quadrature_function_array
    fun_arg::Ptr{Cvoid}
end

struct quadrature_integrate_options
    integrator::quadrature_integrator
    abs_tolerance::Cdouble
    rel_tolerance::Cdouble
    qag_points_per_interval::Csize_t
    max_intervals::Csize_t
end

function quadrature_integrate(__f, __a, __b, options, status, abs_error, workspace_size, workspace)
    @ccall libacc.quadrature_integrate(__f::Ptr{quadrature_integrate_function}, __a::Cdouble, __b::Cdouble, options::Ptr{quadrature_integrate_options}, status::Ptr{quadrature_status}, abs_error::Ptr{Cdouble}, workspace_size::Csize_t, workspace::Ptr{Cvoid})::Cdouble
end

const QUADRATURE_INTEGRATE_QAG_WORKSPACE_PER_INTERVAL = 32

const QUADRATURE_INTEGRATE_QAGS_WORKSPACE_PER_INTERVAL = 152

end # module
