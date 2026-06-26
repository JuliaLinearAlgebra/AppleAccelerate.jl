module LibAccelerate

# Accelerate ships as a macOS *system framework*. Unlike the usual Clang.jl + JLL
# flow, there is no artifact to load — every Mac already has it. We ccall it by
# absolute path, mirroring the `libacc` const in the top-level AppleAccelerate module.
const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

# BNNSGraph opaque handles. In bnns_graph.h these are anonymous structs that all share the
# same `{ void *data; size_t size }` layout. The shim (gen/shims/bnns_graph_shim.h) lets
# Clang.jl's resolver get past their trailing availability attributes, but Clang.jl still
# does not emit a definition for them (identical anonymous bodies), so the generated
# functions reference these names without defining them. We provide the definitions here;
# the field layout matches the C declarations exactly.
struct bnns_graph_t
    data::Ptr{Cvoid}
    size::Csize_t
end

struct bnns_graph_context_t
    data::Ptr{Cvoid}
    size::Csize_t
end

struct bnns_graph_compile_options_t
    data::Ptr{Cvoid}
    size::Csize_t
end

struct bnns_graph_shape_t
    data::Ptr{Cvoid}
    size::Csize_t
end

struct bnns_user_message_data_t
    data::Ptr{Cvoid}
    size::Csize_t
end

# bnns_graph_argument_t holds an anonymous union of three pointer types
# (BNNSTensor* / BNNSNDArrayDescriptor* / void*) followed by a size. Clang.jl drops structs
# with anonymous unions, so define it here. The union collapses to one pointer-sized field.
struct bnns_graph_argument_t
    data::Ptr{Cvoid}            # union { tensor; descriptor; data_ptr }
    data_ptr_size::Csize_t
end


const vDSP_Length = Culong

const vDSP_Stride = Clong

struct DSPComplex
    real::Cfloat
    imag::Cfloat
end

struct DSPDoubleComplex
    real::Cdouble
    imag::Cdouble
end

struct DSPSplitComplex
    realp::Ptr{Cfloat}
    imagp::Ptr{Cfloat}
end

struct DSPDoubleSplitComplex
    realp::Ptr{Cdouble}
    imagp::Ptr{Cdouble}
end

const FFTDirection = Cint

const FFTRadix = Cint

@enum var"##Ctag#277"::Int32 begin
    kFFTDirection_Forward = 1
    kFFTDirection_Inverse = -1
end

@enum var"##Ctag#278"::UInt32 begin
    kFFTRadix2 = 0
    kFFTRadix3 = 1
    kFFTRadix5 = 2
end

@enum var"##Ctag#279"::UInt32 begin
    vDSP_HALF_WINDOW = 1
    vDSP_HANN_DENORM = 0
    vDSP_HANN_NORM = 2
end

struct vDSP_uint24
    bytes::NTuple{3, UInt8}
end

struct vDSP_int24
    bytes::NTuple{3, UInt8}
end

mutable struct OpaqueFFTSetup end

const FFTSetup = Ptr{OpaqueFFTSetup}

mutable struct OpaqueFFTSetupD end

const FFTSetupD = Ptr{OpaqueFFTSetupD}

mutable struct vDSP_biquad_SetupStruct end

const vDSP_biquad_Setup = Ptr{vDSP_biquad_SetupStruct}

mutable struct vDSP_biquad_SetupStructD end

const vDSP_biquad_SetupD = Ptr{vDSP_biquad_SetupStructD}

mutable struct vDSP_biquadm_SetupStruct end

const vDSP_biquadm_Setup = Ptr{vDSP_biquadm_SetupStruct}

mutable struct vDSP_biquadm_SetupStructD end

const vDSP_biquadm_SetupD = Ptr{vDSP_biquadm_SetupStructD}

function vDSP_create_fftsetup(__Log2n, __Radix)
    @ccall libacc.vDSP_create_fftsetup(__Log2n::vDSP_Length, __Radix::FFTRadix)::FFTSetup
end

function vDSP_destroy_fftsetup(__setup)
    @ccall libacc.vDSP_destroy_fftsetup(__setup::FFTSetup)::Cvoid
end

function vDSP_create_fftsetupD(__Log2n, __Radix)
    @ccall libacc.vDSP_create_fftsetupD(__Log2n::vDSP_Length, __Radix::FFTRadix)::FFTSetupD
end

function vDSP_destroy_fftsetupD(__setup)
    @ccall libacc.vDSP_destroy_fftsetupD(__setup::FFTSetupD)::Cvoid
end

function vDSP_biquad_CreateSetup(__Coefficients, __M)
    @ccall libacc.vDSP_biquad_CreateSetup(__Coefficients::Ptr{Cdouble}, __M::vDSP_Length)::vDSP_biquad_Setup
end

function vDSP_biquad_CreateSetupD(__Coefficients, __M)
    @ccall libacc.vDSP_biquad_CreateSetupD(__Coefficients::Ptr{Cdouble}, __M::vDSP_Length)::vDSP_biquad_SetupD
end

function vDSP_biquad_SetCoefficientsDouble(__setup, __coeffs, __start_sec, __nsec)
    @ccall libacc.vDSP_biquad_SetCoefficientsDouble(__setup::vDSP_biquad_Setup, __coeffs::Ptr{Cdouble}, __start_sec::vDSP_Length, __nsec::vDSP_Length)::Cvoid
end

function vDSP_biquad_SetCoefficientsSingle(__setup, __coeffs, __start_sec, __nsec)
    @ccall libacc.vDSP_biquad_SetCoefficientsSingle(__setup::vDSP_biquad_Setup, __coeffs::Ptr{Cfloat}, __start_sec::vDSP_Length, __nsec::vDSP_Length)::Cvoid
end

function vDSP_biquad_DestroySetup(__setup)
    @ccall libacc.vDSP_biquad_DestroySetup(__setup::vDSP_biquad_Setup)::Cvoid
end

function vDSP_biquad_DestroySetupD(__setup)
    @ccall libacc.vDSP_biquad_DestroySetupD(__setup::vDSP_biquad_SetupD)::Cvoid
end

function vDSP_biquadm_CreateSetup(__coeffs, __M, __N)
    @ccall libacc.vDSP_biquadm_CreateSetup(__coeffs::Ptr{Cdouble}, __M::vDSP_Length, __N::vDSP_Length)::vDSP_biquadm_Setup
end

function vDSP_biquadm_CreateSetupD(__coeffs, __M, __N)
    @ccall libacc.vDSP_biquadm_CreateSetupD(__coeffs::Ptr{Cdouble}, __M::vDSP_Length, __N::vDSP_Length)::vDSP_biquadm_SetupD
end

function vDSP_biquadm_DestroySetup(__setup)
    @ccall libacc.vDSP_biquadm_DestroySetup(__setup::vDSP_biquadm_Setup)::Cvoid
end

function vDSP_biquadm_DestroySetupD(__setup)
    @ccall libacc.vDSP_biquadm_DestroySetupD(__setup::vDSP_biquadm_SetupD)::Cvoid
end

function vDSP_biquadm_CopyState(__dest, __src)
    @ccall libacc.vDSP_biquadm_CopyState(__dest::vDSP_biquadm_Setup, __src::Ptr{vDSP_biquadm_SetupStruct})::Cvoid
end

function vDSP_biquadm_CopyStateD(__dest, __src)
    @ccall libacc.vDSP_biquadm_CopyStateD(__dest::vDSP_biquadm_SetupD, __src::Ptr{vDSP_biquadm_SetupStructD})::Cvoid
end

function vDSP_biquadm_ResetState(__setup)
    @ccall libacc.vDSP_biquadm_ResetState(__setup::vDSP_biquadm_Setup)::Cvoid
end

function vDSP_biquadm_ResetStateD(__setup)
    @ccall libacc.vDSP_biquadm_ResetStateD(__setup::vDSP_biquadm_SetupD)::Cvoid
end

function vDSP_biquadm_SetCoefficientsDouble(__setup, __coeffs, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetCoefficientsDouble(__setup::vDSP_biquadm_Setup, __coeffs::Ptr{Cdouble}, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetCoefficientsDoubleD(__setup, __coeffs, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetCoefficientsDoubleD(__setup::vDSP_biquadm_SetupD, __coeffs::Ptr{Cdouble}, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetTargetsDouble(__setup, __targets, __interp_rate, __interp_threshold, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetTargetsDouble(__setup::vDSP_biquadm_Setup, __targets::Ptr{Cdouble}, __interp_rate::Cfloat, __interp_threshold::Cfloat, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetTargetsDoubleD(__setup, __targets, __interp_rate, __interp_threshold, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetTargetsDoubleD(__setup::vDSP_biquadm_SetupD, __targets::Ptr{Cdouble}, __interp_rate::Cdouble, __interp_threshold::Cdouble, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetCoefficientsSingle(__setup, __coeffs, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetCoefficientsSingle(__setup::vDSP_biquadm_Setup, __coeffs::Ptr{Cfloat}, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetCoefficientsSingleD(__setup, __coeffs, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetCoefficientsSingleD(__setup::vDSP_biquadm_SetupD, __coeffs::Ptr{Cfloat}, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetTargetsSingle(__setup, __targets, __interp_rate, __interp_threshold, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetTargetsSingle(__setup::vDSP_biquadm_Setup, __targets::Ptr{Cfloat}, __interp_rate::Cfloat, __interp_threshold::Cfloat, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetTargetsSingleD(__setup, __targets, __interp_rate, __interp_threshold, __start_sec, __start_chn, __nsec, __nchn)
    @ccall libacc.vDSP_biquadm_SetTargetsSingleD(__setup::vDSP_biquadm_SetupD, __targets::Ptr{Cfloat}, __interp_rate::Cdouble, __interp_threshold::Cdouble, __start_sec::vDSP_Length, __start_chn::vDSP_Length, __nsec::vDSP_Length, __nchn::vDSP_Length)::Cvoid
end

function vDSP_biquadm_SetActiveFilters(__setup, __filter_states)
    @ccall libacc.vDSP_biquadm_SetActiveFilters(__setup::vDSP_biquadm_Setup, __filter_states::Ptr{Bool})::Cvoid
end

function vDSP_biquadm_SetActiveFiltersD(__setup, __filter_states)
    @ccall libacc.vDSP_biquadm_SetActiveFiltersD(__setup::vDSP_biquadm_SetupD, __filter_states::Ptr{Bool})::Cvoid
end

function vDSP_ctoz(__C, __IC, __Z, __IZ, __N)
    @ccall libacc.vDSP_ctoz(__C::Ptr{DSPComplex}, __IC::vDSP_Stride, __Z::Ptr{DSPSplitComplex}, __IZ::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_ctozD(__C, __IC, __Z, __IZ, __N)
    @ccall libacc.vDSP_ctozD(__C::Ptr{DSPDoubleComplex}, __IC::vDSP_Stride, __Z::Ptr{DSPDoubleSplitComplex}, __IZ::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_ztoc(__Z, __IZ, __C, __IC, __N)
    @ccall libacc.vDSP_ztoc(__Z::Ptr{DSPSplitComplex}, __IZ::vDSP_Stride, __C::Ptr{DSPComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_ztocD(__Z, __IZ, __C, __IC, __N)
    @ccall libacc.vDSP_ztocD(__Z::Ptr{DSPDoubleSplitComplex}, __IZ::vDSP_Stride, __C::Ptr{DSPDoubleComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_fft_zip(__Setup, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zip(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zipD(__Setup, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zipD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zipt(__Setup, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zipt(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_ziptD(__Setup, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_ziptD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zop(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zopt(__Setup, __A, __IA, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zopt(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zopD(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zopD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zoptD(__Setup, __A, __IA, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zoptD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zrip(__Setup, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zrip(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zripD(__Setup, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zripD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zript(__Setup, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zript(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zriptD(__Setup, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zriptD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zrop(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zrop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zropD(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zropD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zropt(__Setup, __A, __IA, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zropt(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft_zroptD(__Setup, __A, __IA, __C, __IC, __Buffer, __Log2N, __Direction)
    @ccall libacc.vDSP_fft_zroptD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zip(__Setup, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zip(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zipD(__Setup, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zipD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zipt(__Setup, __C, __IC1, __IC0, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zipt(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC1::vDSP_Stride, __IC0::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_ziptD(__Setup, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_ziptD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zop(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zopD(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zopD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zopt(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zopt(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zoptD(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zoptD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zrip(__Setup, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zrip(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zripD(__Setup, __C, __IC0, __IC1, __Log2N0, __Log2N1, __flag)
    @ccall libacc.vDSP_fft2d_zripD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __flag::FFTDirection)::Cvoid
end

function vDSP_fft2d_zript(__Setup, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zript(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zriptD(__Setup, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __flag)
    @ccall libacc.vDSP_fft2d_zriptD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __flag::FFTDirection)::Cvoid
end

function vDSP_fft2d_zrop(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zrop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zropt(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zropt(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zropD(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zropD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft2d_zroptD(__Setup, __A, __IA0, __IA1, __C, __IC0, __IC1, __Buffer, __Log2N0, __Log2N1, __Direction)
    @ccall libacc.vDSP_fft2d_zroptD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA0::vDSP_Stride, __IA1::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC0::vDSP_Stride, __IC1::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N0::vDSP_Length, __Log2N1::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zip(__Setup, __C, __IC, __IM, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zip(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zipD(__Setup, __C, __IC, __IM, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zipD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zipt(__Setup, __C, __IC, __IM, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zipt(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_ziptD(__Setup, __C, __IC, __IM, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_ziptD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zop(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zopD(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zopD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zopt(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zopt(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zoptD(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zoptD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zrip(__Setup, __C, __IC, __IM, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zrip(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zripD(__Setup, __C, __IC, __IM, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zripD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zript(__Setup, __C, __IC, __IM, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zript(__Setup::FFTSetup, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zriptD(__Setup, __C, __IC, __IM, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zriptD(__Setup::FFTSetupD, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IM::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zrop(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zrop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zropt(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zropt(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Buffer::Ptr{DSPSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zropD(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zropD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fftm_zroptD(__Setup, __A, __IA, __IMA, __C, __IC, __IMC, __Buffer, __Log2N, __M, __Direction)
    @ccall libacc.vDSP_fftm_zroptD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __IMA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __IMC::vDSP_Stride, __Buffer::Ptr{DSPDoubleSplitComplex}, __Log2N::vDSP_Length, __M::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft3_zop(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft3_zop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft3_zopD(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft3_zopD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft5_zop(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft5_zop(__Setup::FFTSetup, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_fft5_zopD(__Setup, __A, __IA, __C, __IC, __Log2N, __Direction)
    @ccall libacc.vDSP_fft5_zopD(__Setup::FFTSetupD, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __Log2N::vDSP_Length, __Direction::FFTDirection)::Cvoid
end

function vDSP_biquad(__Setup, __Delay, __X, __IX, __Y, __IY, __N)
    @ccall libacc.vDSP_biquad(__Setup::Ptr{vDSP_biquad_SetupStruct}, __Delay::Ptr{Cfloat}, __X::Ptr{Cfloat}, __IX::vDSP_Stride, __Y::Ptr{Cfloat}, __IY::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_biquadD(__Setup, __Delay, __X, __IX, __Y, __IY, __N)
    @ccall libacc.vDSP_biquadD(__Setup::Ptr{vDSP_biquad_SetupStructD}, __Delay::Ptr{Cdouble}, __X::Ptr{Cdouble}, __IX::vDSP_Stride, __Y::Ptr{Cdouble}, __IY::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_biquadm(__Setup, __X, __IX, __Y, __IY, __N)
    @ccall libacc.vDSP_biquadm(__Setup::vDSP_biquadm_Setup, __X::Ptr{Ptr{Cfloat}}, __IX::vDSP_Stride, __Y::Ptr{Ptr{Cfloat}}, __IY::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_biquadmD(__Setup, __X, __IX, __Y, __IY, __N)
    @ccall libacc.vDSP_biquadmD(__Setup::vDSP_biquadm_SetupD, __X::Ptr{Ptr{Cdouble}}, __IX::vDSP_Stride, __Y::Ptr{Ptr{Cdouble}}, __IY::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_conv(__A, __IA, __F, __IF, __C, __IC, __N, __P)
    @ccall libacc.vDSP_conv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __F::Ptr{Cfloat}, __IF::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_convD(__A, __IA, __F, __IF, __C, __IC, __N, __P)
    @ccall libacc.vDSP_convD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __F::Ptr{Cdouble}, __IF::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zconv(__A, __IA, __F, __IF, __C, __IC, __N, __P)
    @ccall libacc.vDSP_zconv(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __F::Ptr{DSPSplitComplex}, __IF::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zconvD(__A, __IA, __F, __IF, __C, __IC, __N, __P)
    @ccall libacc.vDSP_zconvD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __F::Ptr{DSPDoubleSplitComplex}, __IF::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_f3x3(__A, __NR, __NC, __F, __C)
    @ccall libacc.vDSP_f3x3(__A::Ptr{Cfloat}, __NR::vDSP_Length, __NC::vDSP_Length, __F::Ptr{Cfloat}, __C::Ptr{Cfloat})::Cvoid
end

function vDSP_f3x3D(__A, __NR, __NC, __F, __C)
    @ccall libacc.vDSP_f3x3D(__A::Ptr{Cdouble}, __NR::vDSP_Length, __NC::vDSP_Length, __F::Ptr{Cdouble}, __C::Ptr{Cdouble})::Cvoid
end

function vDSP_f5x5(__A, __NR, __NC, __F, __C)
    @ccall libacc.vDSP_f5x5(__A::Ptr{Cfloat}, __NR::vDSP_Length, __NC::vDSP_Length, __F::Ptr{Cfloat}, __C::Ptr{Cfloat})::Cvoid
end

function vDSP_f5x5D(__A, __NR, __NC, __F, __C)
    @ccall libacc.vDSP_f5x5D(__A::Ptr{Cdouble}, __NR::vDSP_Length, __NC::vDSP_Length, __F::Ptr{Cdouble}, __C::Ptr{Cdouble})::Cvoid
end

function vDSP_imgfir(__A, __NR, __NC, __F, __C, __P, __Q)
    @ccall libacc.vDSP_imgfir(__A::Ptr{Cfloat}, __NR::vDSP_Length, __NC::vDSP_Length, __F::Ptr{Cfloat}, __C::Ptr{Cfloat}, __P::vDSP_Length, __Q::vDSP_Length)::Cvoid
end

function vDSP_imgfirD(__A, __NR, __NC, __F, __C, __P, __Q)
    @ccall libacc.vDSP_imgfirD(__A::Ptr{Cdouble}, __NR::vDSP_Length, __NC::vDSP_Length, __F::Ptr{Cdouble}, __C::Ptr{Cdouble}, __P::vDSP_Length, __Q::vDSP_Length)::Cvoid
end

function vDSP_mtrans(__A, __IA, __C, __IC, __M, __N)
    @ccall libacc.vDSP_mtrans(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length)::Cvoid
end

function vDSP_mtransD(__A, __IA, __C, __IC, __M, __N)
    @ccall libacc.vDSP_mtransD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length)::Cvoid
end

function vDSP_mmul(__A, __IA, __B, __IB, __C, __IC, __M, __N, __P)
    @ccall libacc.vDSP_mmul(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_mmulD(__A, __IA, __B, __IB, __C, __IC, __M, __N, __P)
    @ccall libacc.vDSP_mmulD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmma(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __M, __N, __P)
    @ccall libacc.vDSP_zmma(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmmaD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __M, __N, __P)
    @ccall libacc.vDSP_zmmaD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmms(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __M, __N, __P)
    @ccall libacc.vDSP_zmms(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmmsD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __M, __N, __P)
    @ccall libacc.vDSP_zmmsD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zvmmaa(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __F, __IF, __N)
    @ccall libacc.vDSP_zvmmaa(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __E::Ptr{DSPSplitComplex}, __IE::vDSP_Stride, __F::Ptr{DSPSplitComplex}, __IF::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmmaaD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __F, __IF, __N)
    @ccall libacc.vDSP_zvmmaaD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __E::Ptr{DSPDoubleSplitComplex}, __IE::vDSP_Stride, __F::Ptr{DSPDoubleSplitComplex}, __IF::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zmsm(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __M, __N, __P)
    @ccall libacc.vDSP_zmsm(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmsmD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __M, __N, __P)
    @ccall libacc.vDSP_zmsmD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmmul(__A, __IA, __B, __IB, __C, __IC, __M, __N, __P)
    @ccall libacc.vDSP_zmmul(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zmmulD(__A, __IA, __B, __IB, __C, __IC, __M, __N, __P)
    @ccall libacc.vDSP_zmmulD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __M::vDSP_Length, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_vadd(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vadd(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vaddD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vaddD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vaddi(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vaddi(__A::Ptr{Cint}, __IA::vDSP_Stride, __B::Ptr{Cint}, __IB::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvadd(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvadd(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvaddD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvaddD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvadd(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvadd(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvaddD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvaddD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsub(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vsub(__B::Ptr{Cfloat}, __IB::vDSP_Stride, __A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsubD(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vsubD(__B::Ptr{Cdouble}, __IB::vDSP_Stride, __A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvsub(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvsub(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvsubD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvsubD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmul(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmul(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmulD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmulD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvmul(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvmul(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvmulD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvmulD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vdiv(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vdiv(__B::Ptr{Cfloat}, __IB::vDSP_Stride, __A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vdivD(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vdivD(__B::Ptr{Cdouble}, __IB::vDSP_Stride, __A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vdivi(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vdivi(__B::Ptr{Cint}, __IB::vDSP_Stride, __A::Ptr{Cint}, __IA::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvdiv(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvdiv(__B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvdivD(__B, __IB, __A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvdivD(__B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvdiv(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvdiv(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvdivD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvdivD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmul(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsmul(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmulD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsmulD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsq(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vsq(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsqD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vsqD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vssq(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vssq(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vssqD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vssqD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_distancesq(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_distancesq(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_distancesqD(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_distancesqD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotpr(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_dotpr(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotprD(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_dotprD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_zdotpr(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_zdotpr(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zdotprD(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_zdotprD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zrdotpr(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_zrdotpr(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zrdotprD(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_zrdotprD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_vam(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vam(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vamD(__A, __IA, __B, __IB, __C, __IC, __D, __IDD, __N)
    @ccall libacc.vDSP_vamD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __IDD::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vma(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vma(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmaD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vmaD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvma(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_zvma(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmaD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_zvmaD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmul(__A, __IA, __B, __IB, __C, __IC, __N, __Conjugate)
    @ccall libacc.vDSP_zvmul(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length, __Conjugate::Cint)::Cvoid
end

function vDSP_zvmulD(__A, __IA, __B, __IB, __C, __IC, __N, __Conjugate)
    @ccall libacc.vDSP_zvmulD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length, __Conjugate::Cint)::Cvoid
end

function vDSP_zidotpr(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_zidotpr(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zidotprD(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_zidotprD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zvcma(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_zvcma(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvcmaD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_zvcmaD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvsub(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvsub(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zrvsubD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zrvsubD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vdpsp(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vdpsp(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vspdp(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vspdp(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vabs(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vabs(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vabsD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vabsD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vabsi(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vabsi(__A::Ptr{Cint}, __IA::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvabs(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvabs(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvabsD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvabsD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_veqvi(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_veqvi(__A::Ptr{Cint}, __IA::vDSP_Stride, __B::Ptr{Cint}, __IB::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfill(__A, __C, __IC, __N)
    @ccall libacc.vDSP_vfill(__A::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfillD(__A, __C, __IC, __N)
    @ccall libacc.vDSP_vfillD(__A::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfilli(__A, __C, __IC, __N)
    @ccall libacc.vDSP_vfilli(__A::Ptr{Cint}, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvfill(__A, __C, __IC, __N)
    @ccall libacc.vDSP_zvfill(__A::Ptr{DSPSplitComplex}, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvfillD(__A, __C, __IC, __N)
    @ccall libacc.vDSP_zvfillD(__A::Ptr{DSPDoubleSplitComplex}, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsadd(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsadd(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsaddD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsaddD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsaddi(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsaddi(__A::Ptr{Cint}, __IA::vDSP_Stride, __B::Ptr{Cint}, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsdiv(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsdiv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsdivD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsdivD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsdivi(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsdivi(__A::Ptr{Cint}, __IA::vDSP_Stride, __B::Ptr{Cint}, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zaspec(__A, __C, __N)
    @ccall libacc.vDSP_zaspec(__A::Ptr{DSPSplitComplex}, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_zaspecD(__A, __C, __N)
    @ccall libacc.vDSP_zaspecD(__A::Ptr{DSPDoubleSplitComplex}, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_blkman_window(__C, __N, __Flag)
    @ccall libacc.vDSP_blkman_window(__C::Ptr{Cfloat}, __N::vDSP_Length, __Flag::Cint)::Cvoid
end

function vDSP_blkman_windowD(__C, __N, __Flag)
    @ccall libacc.vDSP_blkman_windowD(__C::Ptr{Cdouble}, __N::vDSP_Length, __Flag::Cint)::Cvoid
end

function vDSP_zcoher(__A, __B, __C, __D, __N)
    @ccall libacc.vDSP_zcoher(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __C::Ptr{DSPSplitComplex}, __D::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_zcoherD(__A, __B, __C, __D, __N)
    @ccall libacc.vDSP_zcoherD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __C::Ptr{DSPDoubleSplitComplex}, __D::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_desamp(__A, __DF, __F, __C, __N, __P)
    @ccall libacc.vDSP_desamp(__A::Ptr{Cfloat}, __DF::vDSP_Stride, __F::Ptr{Cfloat}, __C::Ptr{Cfloat}, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_desampD(__A, __DF, __F, __C, __N, __P)
    @ccall libacc.vDSP_desampD(__A::Ptr{Cdouble}, __DF::vDSP_Stride, __F::Ptr{Cdouble}, __C::Ptr{Cdouble}, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zrdesamp(__A, __DF, __F, __C, __N, __P)
    @ccall libacc.vDSP_zrdesamp(__A::Ptr{DSPSplitComplex}, __DF::vDSP_Stride, __F::Ptr{Cfloat}, __C::Ptr{DSPSplitComplex}, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_zrdesampD(__A, __DF, __F, __C, __N, __P)
    @ccall libacc.vDSP_zrdesampD(__A::Ptr{DSPDoubleSplitComplex}, __DF::vDSP_Stride, __F::Ptr{Cdouble}, __C::Ptr{DSPDoubleSplitComplex}, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_ztrans(__A, __B, __C, __N)
    @ccall libacc.vDSP_ztrans(__A::Ptr{Cfloat}, __B::Ptr{DSPSplitComplex}, __C::Ptr{DSPSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_ztransD(__A, __B, __C, __N)
    @ccall libacc.vDSP_ztransD(__A::Ptr{Cdouble}, __B::Ptr{DSPDoubleSplitComplex}, __C::Ptr{DSPDoubleSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zcspec(__A, __B, __C, __N)
    @ccall libacc.vDSP_zcspec(__A::Ptr{DSPSplitComplex}, __B::Ptr{DSPSplitComplex}, __C::Ptr{DSPSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zcspecD(__A, __B, __C, __N)
    @ccall libacc.vDSP_zcspecD(__A::Ptr{DSPDoubleSplitComplex}, __B::Ptr{DSPDoubleSplitComplex}, __C::Ptr{DSPDoubleSplitComplex}, __N::vDSP_Length)::Cvoid
end

function vDSP_zvcmul(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvcmul(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvcmulD(__A, __IA, __B, __IB, __C, __iC, __N)
    @ccall libacc.vDSP_zvcmulD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __IB::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __iC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvconj(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvconj(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvconjD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvconjD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvzsml(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_zvzsml(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvzsmlD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_zvzsmlD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmags(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvmags(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmagsD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvmagsD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmgsa(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvmgsa(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmgsaD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_zvmgsaD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmov(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvmov(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvmovD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvmovD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvneg(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvneg(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvnegD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvnegD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvphas(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvphas(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvphasD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_zvphasD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvsma(__A, __IA, __B, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_zvsma(__A::Ptr{DSPSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPSplitComplex}, __C::Ptr{DSPSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPSplitComplex}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_zvsmaD(__A, __IA, __B, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_zvsmaD(__A::Ptr{DSPDoubleSplitComplex}, __IA::vDSP_Stride, __B::Ptr{DSPDoubleSplitComplex}, __C::Ptr{DSPDoubleSplitComplex}, __IC::vDSP_Stride, __D::Ptr{DSPDoubleSplitComplex}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_deq22(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_deq22(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_deq22D(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_deq22D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_hamm_window(__C, __N, __Flag)
    @ccall libacc.vDSP_hamm_window(__C::Ptr{Cfloat}, __N::vDSP_Length, __Flag::Cint)::Cvoid
end

function vDSP_hamm_windowD(__C, __N, __Flag)
    @ccall libacc.vDSP_hamm_windowD(__C::Ptr{Cdouble}, __N::vDSP_Length, __Flag::Cint)::Cvoid
end

function vDSP_hann_window(__C, __N, __Flag)
    @ccall libacc.vDSP_hann_window(__C::Ptr{Cfloat}, __N::vDSP_Length, __Flag::Cint)::Cvoid
end

function vDSP_hann_windowD(__C, __N, __Flag)
    @ccall libacc.vDSP_hann_windowD(__C::Ptr{Cdouble}, __N::vDSP_Length, __Flag::Cint)::Cvoid
end

function vDSP_maxmgv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_maxmgv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxmgvD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_maxmgvD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxmgvi(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_maxmgvi(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxmgviD(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_maxmgviD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_maxv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxvD(__A, __I, __C, __N)
    @ccall libacc.vDSP_maxvD(__A::Ptr{Cdouble}, __I::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxvi(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_maxvi(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_maxviD(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_maxviD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_meamgv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_meamgv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_meamgvD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_meamgvD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_meanv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_meanv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_meanvD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_meanvD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_measqv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_measqv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_measqvD(__A, __I, __C, __N)
    @ccall libacc.vDSP_measqvD(__A::Ptr{Cdouble}, __I::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_minmgv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_minmgv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_minmgvD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_minmgvD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_minmgvi(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_minmgvi(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_minmgviD(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_minmgviD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_minv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_minv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_minvD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_minvD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_minvi(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_minvi(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_minviD(__A, __IA, __C, __I, __N)
    @ccall libacc.vDSP_minviD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __I::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_mmov(__A, __C, __M, __N, __TA, __TC)
    @ccall libacc.vDSP_mmov(__A::Ptr{Cfloat}, __C::Ptr{Cfloat}, __M::vDSP_Length, __N::vDSP_Length, __TA::vDSP_Length, __TC::vDSP_Length)::Cvoid
end

function vDSP_mmovD(__A, __C, __M, __N, __TA, __TC)
    @ccall libacc.vDSP_mmovD(__A::Ptr{Cdouble}, __C::Ptr{Cdouble}, __M::vDSP_Length, __N::vDSP_Length, __TA::vDSP_Length, __TC::vDSP_Length)::Cvoid
end

function vDSP_mvessq(__A, __IA, __C, __N)
    @ccall libacc.vDSP_mvessq(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_mvessqD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_mvessqD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_nzcros(__A, __IA, __B, __C, __D, __N)
    @ccall libacc.vDSP_nzcros(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::vDSP_Length, __C::Ptr{vDSP_Length}, __D::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_nzcrosD(__A, __IA, __B, __C, __D, __N)
    @ccall libacc.vDSP_nzcrosD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::vDSP_Length, __C::Ptr{vDSP_Length}, __D::Ptr{vDSP_Length}, __N::vDSP_Length)::Cvoid
end

function vDSP_polar(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_polar(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_polarD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_polarD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_rect(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_rect(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_rectD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_rectD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_rmsqv(__A, __IA, __C, __N)
    @ccall libacc.vDSP_rmsqv(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_rmsqvD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_rmsqvD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_svdiv(__A, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_svdiv(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_svdivD(__A, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_svdivD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_sve(__A, __I, __C, __N)
    @ccall libacc.vDSP_sve(__A::Ptr{Cfloat}, __I::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_sveD(__A, __I, __C, __N)
    @ccall libacc.vDSP_sveD(__A::Ptr{Cdouble}, __I::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_svemg(__A, __IA, __C, __N)
    @ccall libacc.vDSP_svemg(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_svemgD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_svemgD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_svesq(__A, __IA, __C, __N)
    @ccall libacc.vDSP_svesq(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_svesqD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_svesqD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_sve_svesq(__A, __IA, __Sum, __SumOfSquares, __N)
    @ccall libacc.vDSP_sve_svesq(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __Sum::Ptr{Cfloat}, __SumOfSquares::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_sve_svesqD(__A, __IA, __Sum, __SumOfSquares, __N)
    @ccall libacc.vDSP_sve_svesqD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __Sum::Ptr{Cdouble}, __SumOfSquares::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_normalize(__A, __IA, __C, __IC, __Mean, __StandardDeviation, __N)
    @ccall libacc.vDSP_normalize(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __Mean::Ptr{Cfloat}, __StandardDeviation::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_normalizeD(__A, __IA, __C, __IC, __Mean, __StandardDeviation, __N)
    @ccall libacc.vDSP_normalizeD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __Mean::Ptr{Cdouble}, __StandardDeviation::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_svs(__A, __IA, __C, __N)
    @ccall libacc.vDSP_svs(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_svsD(__A, __IA, __C, __N)
    @ccall libacc.vDSP_svsD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_vaam(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vaam(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vaamD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vaamD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vasbm(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vasbm(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vasbmD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vasbmD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vasm(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vasm(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vasmD(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vasmD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vavlin(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vavlin(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vavlinD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vavlinD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vclip(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vclip(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vclipD(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vclipD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vclipc(__A, __IA, __B, __C, __D, __ID, __N, __NLow, __NHigh)
    @ccall libacc.vDSP_vclipc(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length, __NLow::Ptr{vDSP_Length}, __NHigh::Ptr{vDSP_Length})::Cvoid
end

function vDSP_vclipcD(__A, __IA, __B, __C, __D, __ID, __N, __NLow, __NHigh)
    @ccall libacc.vDSP_vclipcD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length, __NLow::Ptr{vDSP_Length}, __NHigh::Ptr{vDSP_Length})::Cvoid
end

function vDSP_vclr(__C, __IC, __N)
    @ccall libacc.vDSP_vclr(__C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vclrD(__C, __IC, __N)
    @ccall libacc.vDSP_vclrD(__C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vcmprs(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vcmprs(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vcmprsD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vcmprsD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vdbcon(__A, __IA, __B, __C, __IC, __N, __F)
    @ccall libacc.vDSP_vdbcon(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __F::Cuint)::Cvoid
end

function vDSP_vdbconD(__A, __IA, __B, __C, __IC, __N, __F)
    @ccall libacc.vDSP_vdbconD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __F::Cuint)::Cvoid
end

function vDSP_vdist(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vdist(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vdistD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vdistD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_venvlp(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_venvlp(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_venvlpD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_venvlpD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfix8(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfix8(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfix8D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfix8D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfix16(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfix16(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cshort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfix16D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfix16D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cshort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfix32(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfix32(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfix32D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfix32D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixu8(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixu8(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cuchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixu8D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixu8D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cuchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixu16(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixu16(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cushort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixu16D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixu16D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cushort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixu32(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixu32(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cuint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixu32D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixu32D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cuint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmfixu24(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsmfixu24(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{vDSP_uint24}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmfix24(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsmfix24(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{vDSP_int24}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu24(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu24(__A::Ptr{vDSP_uint24}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt24(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt24(__A::Ptr{vDSP_int24}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltsmu24(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vfltsmu24(__A::Ptr{vDSP_uint24}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltsm24(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vfltsm24(__A::Ptr{vDSP_int24}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixr8(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixr8(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixr8D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixr8D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixr16(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixr16(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cshort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixr16D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixr16D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cshort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixr32(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixr32(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixr32D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixr32D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixru8(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixru8(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cuchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixru8D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixru8D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cuchar}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixru16(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixru16(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cushort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixru16D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixru16D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cushort}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixru32(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixru32(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cuint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfixru32D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfixru32D(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cuint}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt8(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt8(__A::Ptr{Cchar}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt8D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt8D(__A::Ptr{Cchar}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt16(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt16(__A::Ptr{Cshort}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt16D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt16D(__A::Ptr{Cshort}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt32(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt32(__A::Ptr{Cint}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vflt32D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vflt32D(__A::Ptr{Cint}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu8(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu8(__A::Ptr{Cuchar}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu8D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu8D(__A::Ptr{Cuchar}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu16(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu16(__A::Ptr{Cushort}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu16D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu16D(__A::Ptr{Cushort}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu32(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu32(__A::Ptr{Cuint}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfltu32D(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfltu32D(__A::Ptr{Cuint}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfrac(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfrac(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vfracD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vfracD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgathr(__A, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vgathr(__A::Ptr{Cfloat}, __B::Ptr{vDSP_Length}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgathrD(__A, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vgathrD(__A::Ptr{Cdouble}, __B::Ptr{vDSP_Length}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgathra(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vgathra(__A::Ptr{Ptr{Cfloat}}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgathraD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vgathraD(__A::Ptr{Ptr{Cdouble}}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgen(__A, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vgen(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgenD(__A, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vgenD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vgenp(__A, __IA, __B, __IB, __C, __IC, __N, __M)
    @ccall libacc.vDSP_vgenp(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __M::vDSP_Length)::Cvoid
end

function vDSP_vgenpD(__A, __IA, __B, __IB, __C, __IC, __N, __M)
    @ccall libacc.vDSP_vgenpD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __M::vDSP_Length)::Cvoid
end

function vDSP_viclip(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_viclip(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_viclipD(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_viclipD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vindex(__A, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vindex(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vindexD(__A, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vindexD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vintb(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vintb(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vintbD(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vintbD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vlim(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vlim(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vlimD(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vlimD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vlint(__A, __B, __IB, __C, __IC, __N, __M)
    @ccall libacc.vDSP_vlint(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __M::vDSP_Length)::Cvoid
end

function vDSP_vlintD(__A, __B, __IB, __C, __IC, __N, __M)
    @ccall libacc.vDSP_vlintD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __M::vDSP_Length)::Cvoid
end

function vDSP_vmax(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmax(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmaxD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmaxD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmaxmg(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmaxmg(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmaxmgD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmaxmgD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vswmax(__A, __IA, __C, __IC, __N, __WindowLength)
    @ccall libacc.vDSP_vswmax(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __WindowLength::vDSP_Length)::Cvoid
end

function vDSP_vswmaxD(__A, __IA, __C, __IC, __N, __WindowLength)
    @ccall libacc.vDSP_vswmaxD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __WindowLength::vDSP_Length)::Cvoid
end

function vDSP_vmin(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vmin(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vminD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vminD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vminmg(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vminmg(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vminmgD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vminmgD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmma(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vmma(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmmaD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vmmaD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmmsb(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vmmsb(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmmsbD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vmmsbD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmsa(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vmsa(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmsaD(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vmsaD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmsb(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vmsb(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vmsbD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vmsbD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vnabs(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vnabs(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vnabsD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vnabsD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vneg(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vneg(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vnegD(__A, __IA, __C, __IC, __N)
    @ccall libacc.vDSP_vnegD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vpoly(__A, __IA, __B, __IB, __C, __IC, __N, __P)
    @ccall libacc.vDSP_vpoly(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_vpolyD(__A, __IA, __B, __IB, __C, __IC, __N, __P)
    @ccall libacc.vDSP_vpolyD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_vpythg(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vpythg(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vpythgD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vpythgD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vqint(__A, __B, __IB, __C, __IC, __N, __M)
    @ccall libacc.vDSP_vqint(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __M::vDSP_Length)::Cvoid
end

function vDSP_vqintD(__A, __B, __IB, __C, __IC, __N, __M)
    @ccall libacc.vDSP_vqintD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __M::vDSP_Length)::Cvoid
end

function vDSP_vramp(__A, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vramp(__A::Ptr{Cfloat}, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampD(__A, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vrampD(__A::Ptr{Cdouble}, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrsum(__A, __IA, __S, __C, __IC, __N)
    @ccall libacc.vDSP_vrsum(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __S::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrsumD(__A, __IA, __S, __C, __IC, __N)
    @ccall libacc.vDSP_vrsumD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __S::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrvrs(__C, __IC, __N)
    @ccall libacc.vDSP_vrvrs(__C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrvrsD(__C, __IC, __N)
    @ccall libacc.vDSP_vrvrsD(__C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsbm(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vsbm(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsbmD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vsbmD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsbsbm(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vsbsbm(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsbsbmD(__A, __IA, __B, __IB, __C, __IC, __D, __ID, __E, __IE, __N)
    @ccall libacc.vDSP_vsbsbmD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsbsm(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vsbsm(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsbsmD(__A, __IA, __B, __IB, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vsbsmD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsimps(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsimps(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsimpsD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vsimpsD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsma(__A, __IA, __B, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vsma(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmaD(__A, __IA, __B, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vsmaD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmsa(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vsmsa(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmsaD(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vsmsaD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmsb(__A, __IA, __B, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vsmsb(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmsbD(__A, __IA, __B, __C, __IC, __D, __ID, __N)
    @ccall libacc.vDSP_vsmsbD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmsma(__A, __IA, __B, __C, __IC, __D, __E, __IE, __N)
    @ccall libacc.vDSP_vsmsma(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __D::Ptr{Cfloat}, __E::Ptr{Cfloat}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsmsmaD(__A, __IA, __B, __C, __IC, __D, __E, __IE, __N)
    @ccall libacc.vDSP_vsmsmaD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __D::Ptr{Cdouble}, __E::Ptr{Cdouble}, __IE::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vsort(__C, __N, __Order)
    @ccall libacc.vDSP_vsort(__C::Ptr{Cfloat}, __N::vDSP_Length, __Order::Cint)::Cvoid
end

function vDSP_vsortD(__C, __N, __Order)
    @ccall libacc.vDSP_vsortD(__C::Ptr{Cdouble}, __N::vDSP_Length, __Order::Cint)::Cvoid
end

function vDSP_vsorti(__C, __I, __Temporary, __N, __Order)
    @ccall libacc.vDSP_vsorti(__C::Ptr{Cfloat}, __I::Ptr{vDSP_Length}, __Temporary::Ptr{vDSP_Length}, __N::vDSP_Length, __Order::Cint)::Cvoid
end

function vDSP_vsortiD(__C, __I, __Temporary, __N, __Order)
    @ccall libacc.vDSP_vsortiD(__C::Ptr{Cdouble}, __I::Ptr{vDSP_Length}, __Temporary::Ptr{vDSP_Length}, __N::vDSP_Length, __Order::Cint)::Cvoid
end

function vDSP_vswap(__A, __IA, __B, __IB, __N)
    @ccall libacc.vDSP_vswap(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vswapD(__A, __IA, __B, __IB, __N)
    @ccall libacc.vDSP_vswapD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vswsum(__A, __IA, __C, __IC, __N, __P)
    @ccall libacc.vDSP_vswsum(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_vswsumD(__A, __IA, __C, __IC, __N, __P)
    @ccall libacc.vDSP_vswsumD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length, __P::vDSP_Length)::Cvoid
end

function vDSP_vtabi(__A, __IA, __S1, __S2, __C, __M, __D, __ID, __N)
    @ccall libacc.vDSP_vtabi(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __S1::Ptr{Cfloat}, __S2::Ptr{Cfloat}, __C::Ptr{Cfloat}, __M::vDSP_Length, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vtabiD(__A, __IA, __S1, __S2, __C, __M, __D, __ID, __N)
    @ccall libacc.vDSP_vtabiD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __S1::Ptr{Cdouble}, __S2::Ptr{Cdouble}, __C::Ptr{Cdouble}, __M::vDSP_Length, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vthr(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vthr(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vthrD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vthrD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vthres(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vthres(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vthresD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vthresD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vthrsc(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vthrsc(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __D::Ptr{Cfloat}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vthrscD(__A, __IA, __B, __C, __D, __ID, __N)
    @ccall libacc.vDSP_vthrscD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __D::Ptr{Cdouble}, __ID::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vtmerg(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vtmerg(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vtmergD(__A, __IA, __B, __IB, __C, __IC, __N)
    @ccall libacc.vDSP_vtmergD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vtrapz(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vtrapz(__A::Ptr{Cfloat}, __IA::vDSP_Stride, __B::Ptr{Cfloat}, __C::Ptr{Cfloat}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vtrapzD(__A, __IA, __B, __C, __IC, __N)
    @ccall libacc.vDSP_vtrapzD(__A::Ptr{Cdouble}, __IA::vDSP_Stride, __B::Ptr{Cdouble}, __C::Ptr{Cdouble}, __IC::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_wiener(__L, __A, __C, __F, __P, __Flag, __Error)
    @ccall libacc.vDSP_wiener(__L::vDSP_Length, __A::Ptr{Cfloat}, __C::Ptr{Cfloat}, __F::Ptr{Cfloat}, __P::Ptr{Cfloat}, __Flag::Cint, __Error::Ptr{Cint})::Cvoid
end

function vDSP_wienerD(__L, __A, __C, __F, __P, __Flag, __Error)
    @ccall libacc.vDSP_wienerD(__L::vDSP_Length, __A::Ptr{Cdouble}, __C::Ptr{Cdouble}, __F::Ptr{Cdouble}, __P::Ptr{Cdouble}, __Flag::Cint, __Error::Ptr{Cint})::Cvoid
end

function vDSP_FFT16_copv(__Output, __Input, __Direction)
    @ccall libacc.vDSP_FFT16_copv(__Output::Ptr{Cfloat}, __Input::Ptr{Cfloat}, __Direction::FFTDirection)::Cvoid
end

function vDSP_FFT32_copv(__Output, __Input, __Direction)
    @ccall libacc.vDSP_FFT32_copv(__Output::Ptr{Cfloat}, __Input::Ptr{Cfloat}, __Direction::FFTDirection)::Cvoid
end

function vDSP_FFT16_zopv(__Or, __Oi, __Ir, __Ii, __Direction)
    @ccall libacc.vDSP_FFT16_zopv(__Or::Ptr{Cfloat}, __Oi::Ptr{Cfloat}, __Ir::Ptr{Cfloat}, __Ii::Ptr{Cfloat}, __Direction::FFTDirection)::Cvoid
end

function vDSP_FFT32_zopv(__Or, __Oi, __Ir, __Ii, __Direction)
    @ccall libacc.vDSP_FFT32_zopv(__Or::Ptr{Cfloat}, __Oi::Ptr{Cfloat}, __Ir::Ptr{Cfloat}, __Ii::Ptr{Cfloat}, __Direction::FFTDirection)::Cvoid
end

mutable struct vDSP_DFT_SetupStruct end

const vDSP_DFT_Setup = Ptr{vDSP_DFT_SetupStruct}

mutable struct vDSP_DFT_SetupStructD end

const vDSP_DFT_SetupD = Ptr{vDSP_DFT_SetupStructD}

mutable struct vDSP_DFT_Interleaved_SetupStruct end

const vDSP_DFT_Interleaved_Setup = Ptr{vDSP_DFT_Interleaved_SetupStruct}

mutable struct vDSP_DFT_Interleaved_SetupStructD end

const vDSP_DFT_Interleaved_SetupD = Ptr{vDSP_DFT_Interleaved_SetupStructD}

function vDSP_DFT_CreateSetup(__Previous, __Length)
    @ccall libacc.vDSP_DFT_CreateSetup(__Previous::vDSP_DFT_Setup, __Length::vDSP_Length)::vDSP_DFT_Setup
end

function vDSP_DFT_zop_CreateSetup(__Previous, __Length, __Direction)
    @ccall libacc.vDSP_DFT_zop_CreateSetup(__Previous::vDSP_DFT_Setup, __Length::vDSP_Length, __Direction::Cint)::vDSP_DFT_Setup
end

function vDSP_DFT_zop_CreateSetupD(__Previous, __Length, __Direction)
    @ccall libacc.vDSP_DFT_zop_CreateSetupD(__Previous::vDSP_DFT_SetupD, __Length::vDSP_Length, __Direction::Cint)::vDSP_DFT_SetupD
end

function vDSP_DFT_zrop_CreateSetup(__Previous, __Length, __Direction)
    @ccall libacc.vDSP_DFT_zrop_CreateSetup(__Previous::vDSP_DFT_Setup, __Length::vDSP_Length, __Direction::Cint)::vDSP_DFT_Setup
end

function vDSP_DFT_zrop_CreateSetupD(__Previous, __Length, __Direction)
    @ccall libacc.vDSP_DFT_zrop_CreateSetupD(__Previous::vDSP_DFT_SetupD, __Length::vDSP_Length, __Direction::Cint)::vDSP_DFT_SetupD
end

function vDSP_DFT_DestroySetup(__Setup)
    @ccall libacc.vDSP_DFT_DestroySetup(__Setup::vDSP_DFT_Setup)::Cvoid
end

function vDSP_DFT_DestroySetupD(__Setup)
    @ccall libacc.vDSP_DFT_DestroySetupD(__Setup::vDSP_DFT_SetupD)::Cvoid
end

function vDSP_DFT_zop(__Setup, __Ir, __Ii, __Is, __Or, __Oi, __Os, __Direction)
    @ccall libacc.vDSP_DFT_zop(__Setup::Ptr{vDSP_DFT_SetupStruct}, __Ir::Ptr{Cfloat}, __Ii::Ptr{Cfloat}, __Is::vDSP_Stride, __Or::Ptr{Cfloat}, __Oi::Ptr{Cfloat}, __Os::vDSP_Stride, __Direction::Cint)::Cvoid
end

function vDSP_DFT_Execute(__Setup, __Ir, __Ii, __Or, __Oi)
    @ccall libacc.vDSP_DFT_Execute(__Setup::Ptr{vDSP_DFT_SetupStruct}, __Ir::Ptr{Cfloat}, __Ii::Ptr{Cfloat}, __Or::Ptr{Cfloat}, __Oi::Ptr{Cfloat})::Cvoid
end

function vDSP_DFT_ExecuteD(__Setup, __Ir, __Ii, __Or, __Oi)
    @ccall libacc.vDSP_DFT_ExecuteD(__Setup::Ptr{vDSP_DFT_SetupStructD}, __Ir::Ptr{Cdouble}, __Ii::Ptr{Cdouble}, __Or::Ptr{Cdouble}, __Oi::Ptr{Cdouble})::Cvoid
end

function vDSP_DCT_CreateSetup(__Previous, __Length, __Type)
    @ccall libacc.vDSP_DCT_CreateSetup(__Previous::vDSP_DFT_Setup, __Length::vDSP_Length, __Type::Cint)::vDSP_DFT_Setup
end

function vDSP_DCT_Execute(__Setup, __Input, __Output)
    @ccall libacc.vDSP_DCT_Execute(__Setup::Ptr{vDSP_DFT_SetupStruct}, __Input::Ptr{Cfloat}, __Output::Ptr{Cfloat})::Cvoid
end

function vDSP_DFT_Interleaved_CreateSetup(Previous, Length, Direction, RealtoComplex)
    @ccall libacc.vDSP_DFT_Interleaved_CreateSetup(Previous::vDSP_DFT_Interleaved_Setup, Length::vDSP_Length, Direction::Cint, RealtoComplex::Cint)::vDSP_DFT_Interleaved_Setup
end

function vDSP_DFT_Interleaved_CreateSetupD(Previous, Length, Direction, RealtoComplex)
    @ccall libacc.vDSP_DFT_Interleaved_CreateSetupD(Previous::vDSP_DFT_Interleaved_SetupD, Length::vDSP_Length, Direction::Cint, RealtoComplex::Cint)::vDSP_DFT_Interleaved_SetupD
end

function vDSP_DFT_Interleaved_Execute(Setup, Iri, Ori)
    @ccall libacc.vDSP_DFT_Interleaved_Execute(Setup::vDSP_DFT_Interleaved_Setup, Iri::Ptr{DSPComplex}, Ori::Ptr{DSPComplex})::Cvoid
end

function vDSP_DFT_Interleaved_ExecuteD(Setup, Iri, Ori)
    @ccall libacc.vDSP_DFT_Interleaved_ExecuteD(Setup::vDSP_DFT_Interleaved_SetupD, Iri::Ptr{DSPDoubleComplex}, Ori::Ptr{DSPDoubleComplex})::Cvoid
end

function vDSP_DFT_Interleaved_DestroySetup(Setup)
    @ccall libacc.vDSP_DFT_Interleaved_DestroySetup(Setup::vDSP_DFT_Interleaved_Setup)::Cvoid
end

function vDSP_DFT_Interleaved_DestroySetupD(Setup)
    @ccall libacc.vDSP_DFT_Interleaved_DestroySetupD(Setup::vDSP_DFT_Interleaved_SetupD)::Cvoid
end

function vDSP_dotpr2(__A0, __IA0, __A1, __IA1, __B, __IB, __C0, __C1, __N)
    @ccall libacc.vDSP_dotpr2(__A0::Ptr{Cfloat}, __IA0::vDSP_Stride, __A1::Ptr{Cfloat}, __IA1::vDSP_Stride, __B::Ptr{Cfloat}, __IB::vDSP_Stride, __C0::Ptr{Cfloat}, __C1::Ptr{Cfloat}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotpr2D(__A0, __IA0, __A1, __IA1, __B, __IB, __C0, __C1, __N)
    @ccall libacc.vDSP_dotpr2D(__A0::Ptr{Cdouble}, __IA0::vDSP_Stride, __A1::Ptr{Cdouble}, __IA1::vDSP_Stride, __B::Ptr{Cdouble}, __IB::vDSP_Stride, __C0::Ptr{Cdouble}, __C1::Ptr{Cdouble}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotpr_s1_15(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_dotpr_s1_15(__A::Ptr{Cshort}, __IA::vDSP_Stride, __B::Ptr{Cshort}, __IB::vDSP_Stride, __C::Ptr{Cshort}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotpr2_s1_15(__A0, __IA0, __A1, __IA1, __B, __IB, __C0, __C1, __N)
    @ccall libacc.vDSP_dotpr2_s1_15(__A0::Ptr{Cshort}, __IA0::vDSP_Stride, __A1::Ptr{Cshort}, __IA1::vDSP_Stride, __B::Ptr{Cshort}, __IB::vDSP_Stride, __C0::Ptr{Cshort}, __C1::Ptr{Cshort}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotpr_s8_24(__A, __IA, __B, __IB, __C, __N)
    @ccall libacc.vDSP_dotpr_s8_24(__A::Ptr{Cint}, __IA::vDSP_Stride, __B::Ptr{Cint}, __IB::vDSP_Stride, __C::Ptr{Cint}, __N::vDSP_Length)::Cvoid
end

function vDSP_dotpr2_s8_24(__A0, __IA0, __A1, __IA1, __B, __IB, __C0, __C1, __N)
    @ccall libacc.vDSP_dotpr2_s8_24(__A0::Ptr{Cint}, __IA0::vDSP_Stride, __A1::Ptr{Cint}, __IA1::vDSP_Stride, __B::Ptr{Cint}, __IB::vDSP_Stride, __C0::Ptr{Cint}, __C1::Ptr{Cint}, __N::vDSP_Length)::Cvoid
end

function vDSP_vaddsub(__I0, __I0S, __I1, __I1S, __O0, __O0S, __O1, __O1S, __N)
    @ccall libacc.vDSP_vaddsub(__I0::Ptr{Cfloat}, __I0S::vDSP_Stride, __I1::Ptr{Cfloat}, __I1S::vDSP_Stride, __O0::Ptr{Cfloat}, __O0S::vDSP_Stride, __O1::Ptr{Cfloat}, __O1S::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vaddsubD(__I0, __I0S, __I1, __I1S, __O0, __O0S, __O1, __O1S, __N)
    @ccall libacc.vDSP_vaddsubD(__I0::Ptr{Cdouble}, __I0S::vDSP_Stride, __I1::Ptr{Cdouble}, __I1S::vDSP_Stride, __O0::Ptr{Cdouble}, __O0S::vDSP_Stride, __O1::Ptr{Cdouble}, __O1S::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmul(__I::Ptr{Cfloat}, __IS::vDSP_Stride, __Start::Ptr{Cfloat}, __Step::Ptr{Cfloat}, __O::Ptr{Cfloat}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmulD(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmulD(__I::Ptr{Cdouble}, __IS::vDSP_Stride, __Start::Ptr{Cdouble}, __Step::Ptr{Cdouble}, __O::Ptr{Cdouble}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd(__I::Ptr{Cfloat}, __IS::vDSP_Stride, __Start::Ptr{Cfloat}, __Step::Ptr{Cfloat}, __O::Ptr{Cfloat}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladdD(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmuladdD(__I::Ptr{Cdouble}, __IS::vDSP_Stride, __Start::Ptr{Cdouble}, __Step::Ptr{Cdouble}, __O::Ptr{Cdouble}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul2(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmul2(__I0::Ptr{Cfloat}, __I1::Ptr{Cfloat}, __IS::vDSP_Stride, __Start::Ptr{Cfloat}, __Step::Ptr{Cfloat}, __O0::Ptr{Cfloat}, __O1::Ptr{Cfloat}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul2D(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmul2D(__I0::Ptr{Cdouble}, __I1::Ptr{Cdouble}, __IS::vDSP_Stride, __Start::Ptr{Cdouble}, __Step::Ptr{Cdouble}, __O0::Ptr{Cdouble}, __O1::Ptr{Cdouble}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd2(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd2(__I0::Ptr{Cfloat}, __I1::Ptr{Cfloat}, __IS::vDSP_Stride, __Start::Ptr{Cfloat}, __Step::Ptr{Cfloat}, __O0::Ptr{Cfloat}, __O1::Ptr{Cfloat}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd2D(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd2D(__I0::Ptr{Cdouble}, __I1::Ptr{Cdouble}, __IS::vDSP_Stride, __Start::Ptr{Cdouble}, __Step::Ptr{Cdouble}, __O0::Ptr{Cdouble}, __O1::Ptr{Cdouble}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul_s1_15(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmul_s1_15(__I::Ptr{Cshort}, __IS::vDSP_Stride, __Start::Ptr{Cshort}, __Step::Ptr{Cshort}, __O::Ptr{Cshort}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd_s1_15(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd_s1_15(__I::Ptr{Cshort}, __IS::vDSP_Stride, __Start::Ptr{Cshort}, __Step::Ptr{Cshort}, __O::Ptr{Cshort}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul2_s1_15(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmul2_s1_15(__I0::Ptr{Cshort}, __I1::Ptr{Cshort}, __IS::vDSP_Stride, __Start::Ptr{Cshort}, __Step::Ptr{Cshort}, __O0::Ptr{Cshort}, __O1::Ptr{Cshort}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd2_s1_15(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd2_s1_15(__I0::Ptr{Cshort}, __I1::Ptr{Cshort}, __IS::vDSP_Stride, __Start::Ptr{Cshort}, __Step::Ptr{Cshort}, __O0::Ptr{Cshort}, __O1::Ptr{Cshort}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul_s8_24(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmul_s8_24(__I::Ptr{Cint}, __IS::vDSP_Stride, __Start::Ptr{Cint}, __Step::Ptr{Cint}, __O::Ptr{Cint}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd_s8_24(__I, __IS, __Start, __Step, __O, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd_s8_24(__I::Ptr{Cint}, __IS::vDSP_Stride, __Start::Ptr{Cint}, __Step::Ptr{Cint}, __O::Ptr{Cint}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmul2_s8_24(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmul2_s8_24(__I0::Ptr{Cint}, __I1::Ptr{Cint}, __IS::vDSP_Stride, __Start::Ptr{Cint}, __Step::Ptr{Cint}, __O0::Ptr{Cint}, __O1::Ptr{Cint}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

function vDSP_vrampmuladd2_s8_24(__I0, __I1, __IS, __Start, __Step, __O0, __O1, __OS, __N)
    @ccall libacc.vDSP_vrampmuladd2_s8_24(__I0::Ptr{Cint}, __I1::Ptr{Cint}, __IS::vDSP_Stride, __Start::Ptr{Cint}, __Step::Ptr{Cint}, __O0::Ptr{Cint}, __O1::Ptr{Cint}, __OS::vDSP_Stride, __N::vDSP_Length)::Cvoid
end

@enum var"##Ctag#280"::Int32 begin
    FFT_FORWARD = 1
    FFT_INVERSE = -1
end

@enum var"##Ctag#281"::UInt32 begin
    FFT_RADIX2 = 0
    FFT_RADIX3 = 1
    FFT_RADIX5 = 2
end

const COMPLEX = DSPComplex

const COMPLEX_SPLIT = DSPSplitComplex

const DOUBLE_COMPLEX = DSPDoubleComplex

const DOUBLE_COMPLEX_SPLIT = DSPDoubleSplitComplex

const __float_complex_t = ComplexF32

const __double_complex_t = ComplexF64

function vvrecf(arg1, arg2, arg3)
    @ccall libacc.vvrecf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvrec(arg1, arg2, arg3)
    @ccall libacc.vvrec(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvdivf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvdivf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvdiv(arg1, arg2, arg3, arg4)
    @ccall libacc.vvdiv(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvsqrtf(arg1, arg2, arg3)
    @ccall libacc.vvsqrtf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvsqrt(arg1, arg2, arg3)
    @ccall libacc.vvsqrt(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvcbrtf(arg1, arg2, arg3)
    @ccall libacc.vvcbrtf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvcbrt(arg1, arg2, arg3)
    @ccall libacc.vvcbrt(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvrsqrtf(arg1, arg2, arg3)
    @ccall libacc.vvrsqrtf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvrsqrt(arg1, arg2, arg3)
    @ccall libacc.vvrsqrt(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvexpf(arg1, arg2, arg3)
    @ccall libacc.vvexpf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvexp(arg1, arg2, arg3)
    @ccall libacc.vvexp(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvexp2f(arg1, arg2, arg3)
    @ccall libacc.vvexp2f(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvexp2(arg1, arg2, arg3)
    @ccall libacc.vvexp2(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvexpm1f(arg1, arg2, arg3)
    @ccall libacc.vvexpm1f(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvexpm1(arg1, arg2, arg3)
    @ccall libacc.vvexpm1(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvlogf(arg1, arg2, arg3)
    @ccall libacc.vvlogf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvlog(arg1, arg2, arg3)
    @ccall libacc.vvlog(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvlog10f(arg1, arg2, arg3)
    @ccall libacc.vvlog10f(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvlog10(arg1, arg2, arg3)
    @ccall libacc.vvlog10(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvlog1pf(arg1, arg2, arg3)
    @ccall libacc.vvlog1pf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvlog1p(arg1, arg2, arg3)
    @ccall libacc.vvlog1p(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvlog2f(arg1, arg2, arg3)
    @ccall libacc.vvlog2f(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvlog2(arg1, arg2, arg3)
    @ccall libacc.vvlog2(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvlogbf(arg1, arg2, arg3)
    @ccall libacc.vvlogbf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvlogb(arg1, arg2, arg3)
    @ccall libacc.vvlogb(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvfabsf(arg1, arg2, arg3)
    @ccall libacc.vvfabsf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvfabs(arg1, arg2, arg3)
    @ccall libacc.vvfabs(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvpowf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvpowf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvpow(arg1, arg2, arg3, arg4)
    @ccall libacc.vvpow(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvpowsf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvpowsf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvpows(arg1, arg2, arg3, arg4)
    @ccall libacc.vvpows(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvsinf(arg1, arg2, arg3)
    @ccall libacc.vvsinf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvsin(arg1, arg2, arg3)
    @ccall libacc.vvsin(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvcosf(arg1, arg2, arg3)
    @ccall libacc.vvcosf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvcos(arg1, arg2, arg3)
    @ccall libacc.vvcos(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvtanf(arg1, arg2, arg3)
    @ccall libacc.vvtanf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvtan(arg1, arg2, arg3)
    @ccall libacc.vvtan(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvasinf(arg1, arg2, arg3)
    @ccall libacc.vvasinf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvasin(arg1, arg2, arg3)
    @ccall libacc.vvasin(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvacosf(arg1, arg2, arg3)
    @ccall libacc.vvacosf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvacos(arg1, arg2, arg3)
    @ccall libacc.vvacos(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvatanf(arg1, arg2, arg3)
    @ccall libacc.vvatanf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvatan(arg1, arg2, arg3)
    @ccall libacc.vvatan(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvatan2f(arg1, arg2, arg3, arg4)
    @ccall libacc.vvatan2f(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvatan2(arg1, arg2, arg3, arg4)
    @ccall libacc.vvatan2(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvsincosf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvsincosf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvsincos(arg1, arg2, arg3, arg4)
    @ccall libacc.vvsincos(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvcosisinf(arg1, arg2, arg3)
    @ccall libacc.vvcosisinf(arg1::Ptr{__float_complex_t}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvcosisin(arg1, arg2, arg3)
    @ccall libacc.vvcosisin(arg1::Ptr{__double_complex_t}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvsinhf(arg1, arg2, arg3)
    @ccall libacc.vvsinhf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvsinh(arg1, arg2, arg3)
    @ccall libacc.vvsinh(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvcoshf(arg1, arg2, arg3)
    @ccall libacc.vvcoshf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvcosh(arg1, arg2, arg3)
    @ccall libacc.vvcosh(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvtanhf(arg1, arg2, arg3)
    @ccall libacc.vvtanhf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvtanh(arg1, arg2, arg3)
    @ccall libacc.vvtanh(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvasinhf(arg1, arg2, arg3)
    @ccall libacc.vvasinhf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvasinh(arg1, arg2, arg3)
    @ccall libacc.vvasinh(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvacoshf(arg1, arg2, arg3)
    @ccall libacc.vvacoshf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvacosh(arg1, arg2, arg3)
    @ccall libacc.vvacosh(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvatanhf(arg1, arg2, arg3)
    @ccall libacc.vvatanhf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvatanh(arg1, arg2, arg3)
    @ccall libacc.vvatanh(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvintf(arg1, arg2, arg3)
    @ccall libacc.vvintf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvint(arg1, arg2, arg3)
    @ccall libacc.vvint(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvnintf(arg1, arg2, arg3)
    @ccall libacc.vvnintf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvnint(arg1, arg2, arg3)
    @ccall libacc.vvnint(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvceilf(arg1, arg2, arg3)
    @ccall libacc.vvceilf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvceil(arg1, arg2, arg3)
    @ccall libacc.vvceil(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvfloorf(arg1, arg2, arg3)
    @ccall libacc.vvfloorf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvfloor(arg1, arg2, arg3)
    @ccall libacc.vvfloor(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvfmodf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvfmodf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvfmod(arg1, arg2, arg3, arg4)
    @ccall libacc.vvfmod(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvremainderf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvremainderf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvremainder(arg1, arg2, arg3, arg4)
    @ccall libacc.vvremainder(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvcopysignf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvcopysignf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvcopysign(arg1, arg2, arg3, arg4)
    @ccall libacc.vvcopysign(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvnextafterf(arg1, arg2, arg3, arg4)
    @ccall libacc.vvnextafterf(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cfloat}, arg4::Ptr{Cint})::Cvoid
end

function vvnextafter(arg1, arg2, arg3, arg4)
    @ccall libacc.vvnextafter(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cdouble}, arg4::Ptr{Cint})::Cvoid
end

function vvsinpif(arg1, arg2, arg3)
    @ccall libacc.vvsinpif(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvsinpi(arg1, arg2, arg3)
    @ccall libacc.vvsinpi(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvcospif(arg1, arg2, arg3)
    @ccall libacc.vvcospif(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvcospi(arg1, arg2, arg3)
    @ccall libacc.vvcospi(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

function vvtanpif(arg1, arg2, arg3)
    @ccall libacc.vvtanpif(arg1::Ptr{Cfloat}, arg2::Ptr{Cfloat}, arg3::Ptr{Cint})::Cvoid
end

function vvtanpi(arg1, arg2, arg3)
    @ccall libacc.vvtanpi(arg1::Ptr{Cdouble}, arg2::Ptr{Cdouble}, arg3::Ptr{Cint})::Cvoid
end

struct vU128
    data::NTuple{16, UInt8}
end

function Base.getproperty(x::Ptr{vU128}, f::Symbol)
    f === :v && return Ptr{vUInt32}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#351"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vU128, f::Symbol)
    r = Ref{vU128}(x)
    ptr = Base.unsafe_convert(Ptr{vU128}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vU128}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vU128, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vS128
    data::NTuple{16, UInt8}
end

function Base.getproperty(x::Ptr{vS128}, f::Symbol)
    f === :v && return Ptr{vUInt32}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#344"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vS128, f::Symbol)
    r = Ref{vS128}(x)
    ptr = Base.unsafe_convert(Ptr{vS128}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vS128}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vS128, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vU256
    data::NTuple{32, UInt8}
end

function Base.getproperty(x::Ptr{vU256}, f::Symbol)
    f === :v && return Ptr{NTuple{2, vUInt32}}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#353"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vU256, f::Symbol)
    r = Ref{vU256}(x)
    ptr = Base.unsafe_convert(Ptr{vU256}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vU256}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vU256, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vS256
    data::NTuple{32, UInt8}
end

function Base.getproperty(x::Ptr{vS256}, f::Symbol)
    f === :v && return Ptr{NTuple{2, vUInt32}}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#357"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vS256, f::Symbol)
    r = Ref{vS256}(x)
    ptr = Base.unsafe_convert(Ptr{vS256}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vS256}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vS256, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vU512
    data::NTuple{64, UInt8}
end

function Base.getproperty(x::Ptr{vU512}, f::Symbol)
    f === :v && return Ptr{NTuple{4, vUInt32}}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#359"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vU512, f::Symbol)
    r = Ref{vU512}(x)
    ptr = Base.unsafe_convert(Ptr{vU512}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vU512}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vU512, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vS512
    data::NTuple{64, UInt8}
end

function Base.getproperty(x::Ptr{vS512}, f::Symbol)
    f === :v && return Ptr{NTuple{4, vUInt32}}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#348"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vS512, f::Symbol)
    r = Ref{vS512}(x)
    ptr = Base.unsafe_convert(Ptr{vS512}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vS512}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vS512, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vU1024
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{vU1024}, f::Symbol)
    f === :v && return Ptr{NTuple{8, vUInt32}}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#355"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vU1024, f::Symbol)
    r = Ref{vU1024}(x)
    ptr = Base.unsafe_convert(Ptr{vU1024}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vU1024}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vU1024, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct vS1024
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{vS1024}, f::Symbol)
    f === :v && return Ptr{NTuple{8, vUInt32}}(x + 0)
    f === :vs && return Ptr{Cvoid}(x + 0)
    f === :s && return Ptr{var"##Ctag#346"}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::vS1024, f::Symbol)
    r = Ref{vS1024}(x)
    ptr = Base.unsafe_convert(Ptr{vS1024}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{vS1024}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::vS1024, private::Bool = false)
    (:v, :vs, :s, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

function vU256Divide(numerator, divisor, result, remainder)
    @ccall libacc.vU256Divide(numerator::Ptr{vU256}, divisor::Ptr{vU256}, result::Ptr{vU256}, remainder::Ptr{vU256})::Cvoid
end

function vS256Divide(numerator, divisor, result, remainder)
    @ccall libacc.vS256Divide(numerator::Ptr{vS256}, divisor::Ptr{vS256}, result::Ptr{vS256}, remainder::Ptr{vS256})::Cvoid
end

function vU512Divide(numerator, divisor, result, remainder)
    @ccall libacc.vU512Divide(numerator::Ptr{vU512}, divisor::Ptr{vU512}, result::Ptr{vU512}, remainder::Ptr{vU512})::Cvoid
end

function vS512Divide(numerator, divisor, result, remainder)
    @ccall libacc.vS512Divide(numerator::Ptr{vS512}, divisor::Ptr{vS512}, result::Ptr{vS512}, remainder::Ptr{vS512})::Cvoid
end

function vU1024Divide(numerator, divisor, result, remainder)
    @ccall libacc.vU1024Divide(numerator::Ptr{vU1024}, divisor::Ptr{vU1024}, result::Ptr{vU1024}, remainder::Ptr{vU1024})::Cvoid
end

function vS1024Divide(numerator, divisor, result, remainder)
    @ccall libacc.vS1024Divide(numerator::Ptr{vS1024}, divisor::Ptr{vS1024}, result::Ptr{vS1024}, remainder::Ptr{vS1024})::Cvoid
end

function vU128FullMultiply(a, b, result)
    @ccall libacc.vU128FullMultiply(a::Ptr{vU128}, b::Ptr{vU128}, result::Ptr{vU256})::Cvoid
end

function vS128FullMultiply(a, b, result)
    @ccall libacc.vS128FullMultiply(a::Ptr{vS128}, b::Ptr{vS128}, result::Ptr{vS256})::Cvoid
end

function vU256FullMultiply(a, b, result)
    @ccall libacc.vU256FullMultiply(a::Ptr{vU256}, b::Ptr{vU256}, result::Ptr{vU512})::Cvoid
end

function vS256FullMultiply(a, b, result)
    @ccall libacc.vS256FullMultiply(a::Ptr{vS256}, b::Ptr{vS256}, result::Ptr{vS512})::Cvoid
end

function vU512FullMultiply(a, b, result)
    @ccall libacc.vU512FullMultiply(a::Ptr{vU512}, b::Ptr{vU512}, result::Ptr{vU1024})::Cvoid
end

function vS512FullMultiply(a, b, result)
    @ccall libacc.vS512FullMultiply(a::Ptr{vS512}, b::Ptr{vS512}, result::Ptr{vS1024})::Cvoid
end

function vU256HalfMultiply(a, b, result)
    @ccall libacc.vU256HalfMultiply(a::Ptr{vU256}, b::Ptr{vU256}, result::Ptr{vU256})::Cvoid
end

function vS256HalfMultiply(a, b, result)
    @ccall libacc.vS256HalfMultiply(a::Ptr{vS256}, b::Ptr{vS256}, result::Ptr{vS256})::Cvoid
end

function vU512HalfMultiply(a, b, result)
    @ccall libacc.vU512HalfMultiply(a::Ptr{vU512}, b::Ptr{vU512}, result::Ptr{vU512})::Cvoid
end

function vS512HalfMultiply(a, b, result)
    @ccall libacc.vS512HalfMultiply(a::Ptr{vS512}, b::Ptr{vS512}, result::Ptr{vS512})::Cvoid
end

function vU1024HalfMultiply(a, b, result)
    @ccall libacc.vU1024HalfMultiply(a::Ptr{vU1024}, b::Ptr{vU1024}, result::Ptr{vU1024})::Cvoid
end

function vS1024HalfMultiply(a, b, result)
    @ccall libacc.vS1024HalfMultiply(a::Ptr{vS1024}, b::Ptr{vS1024}, result::Ptr{vS1024})::Cvoid
end

function vU256Sub(a, b, result)
    @ccall libacc.vU256Sub(a::Ptr{vU256}, b::Ptr{vU256}, result::Ptr{vU256})::Cvoid
end

function vS256Sub(a, b, result)
    @ccall libacc.vS256Sub(a::Ptr{vS256}, b::Ptr{vS256}, result::Ptr{vS256})::Cvoid
end

function vU256SubS(a, b, result)
    @ccall libacc.vU256SubS(a::Ptr{vU256}, b::Ptr{vU256}, result::Ptr{vU256})::Cvoid
end

function vS256SubS(a, b, result)
    @ccall libacc.vS256SubS(a::Ptr{vS256}, b::Ptr{vS256}, result::Ptr{vS256})::Cvoid
end

function vU512Sub(a, b, result)
    @ccall libacc.vU512Sub(a::Ptr{vU512}, b::Ptr{vU512}, result::Ptr{vU512})::Cvoid
end

function vS512Sub(a, b, result)
    @ccall libacc.vS512Sub(a::Ptr{vS512}, b::Ptr{vS512}, result::Ptr{vS512})::Cvoid
end

function vU512SubS(a, b, result)
    @ccall libacc.vU512SubS(a::Ptr{vU512}, b::Ptr{vU512}, result::Ptr{vU512})::Cvoid
end

function vS512SubS(a, b, result)
    @ccall libacc.vS512SubS(a::Ptr{vS512}, b::Ptr{vS512}, result::Ptr{vS512})::Cvoid
end

function vU1024Sub(a, b, result)
    @ccall libacc.vU1024Sub(a::Ptr{vU1024}, b::Ptr{vU1024}, result::Ptr{vU1024})::Cvoid
end

function vS1024Sub(a, b, result)
    @ccall libacc.vS1024Sub(a::Ptr{vS1024}, b::Ptr{vS1024}, result::Ptr{vS1024})::Cvoid
end

function vU1024SubS(a, b, result)
    @ccall libacc.vU1024SubS(a::Ptr{vU1024}, b::Ptr{vU1024}, result::Ptr{vU1024})::Cvoid
end

function vS1024SubS(a, b, result)
    @ccall libacc.vS1024SubS(a::Ptr{vS1024}, b::Ptr{vS1024}, result::Ptr{vS1024})::Cvoid
end

function vU256Neg(a, result)
    @ccall libacc.vU256Neg(a::Ptr{vU256}, result::Ptr{vU256})::Cvoid
end

function vS256Neg(a, result)
    @ccall libacc.vS256Neg(a::Ptr{vS256}, result::Ptr{vS256})::Cvoid
end

function vU512Neg(a, result)
    @ccall libacc.vU512Neg(a::Ptr{vU512}, result::Ptr{vU512})::Cvoid
end

function vS512Neg(a, result)
    @ccall libacc.vS512Neg(a::Ptr{vS512}, result::Ptr{vS512})::Cvoid
end

function vU1024Neg(a, result)
    @ccall libacc.vU1024Neg(a::Ptr{vU1024}, result::Ptr{vU1024})::Cvoid
end

function vS1024Neg(a, result)
    @ccall libacc.vS1024Neg(a::Ptr{vS1024}, result::Ptr{vS1024})::Cvoid
end

function vU256Add(a, b, result)
    @ccall libacc.vU256Add(a::Ptr{vU256}, b::Ptr{vU256}, result::Ptr{vU256})::Cvoid
end

function vS256Add(a, b, result)
    @ccall libacc.vS256Add(a::Ptr{vS256}, b::Ptr{vS256}, result::Ptr{vS256})::Cvoid
end

function vU256AddS(a, b, result)
    @ccall libacc.vU256AddS(a::Ptr{vU256}, b::Ptr{vU256}, result::Ptr{vU256})::Cvoid
end

function vS256AddS(a, b, result)
    @ccall libacc.vS256AddS(a::Ptr{vS256}, b::Ptr{vS256}, result::Ptr{vS256})::Cvoid
end

function vU512Add(a, b, result)
    @ccall libacc.vU512Add(a::Ptr{vU512}, b::Ptr{vU512}, result::Ptr{vU512})::Cvoid
end

function vS512Add(a, b, result)
    @ccall libacc.vS512Add(a::Ptr{vS512}, b::Ptr{vS512}, result::Ptr{vS512})::Cvoid
end

function vU512AddS(a, b, result)
    @ccall libacc.vU512AddS(a::Ptr{vU512}, b::Ptr{vU512}, result::Ptr{vU512})::Cvoid
end

function vS512AddS(a, b, result)
    @ccall libacc.vS512AddS(a::Ptr{vS512}, b::Ptr{vS512}, result::Ptr{vS512})::Cvoid
end

function vU1024Add(a, b, result)
    @ccall libacc.vU1024Add(a::Ptr{vU1024}, b::Ptr{vU1024}, result::Ptr{vU1024})::Cvoid
end

function vS1024Add(a, b, result)
    @ccall libacc.vS1024Add(a::Ptr{vS1024}, b::Ptr{vS1024}, result::Ptr{vS1024})::Cvoid
end

function vU1024AddS(a, b, result)
    @ccall libacc.vU1024AddS(a::Ptr{vU1024}, b::Ptr{vU1024}, result::Ptr{vU1024})::Cvoid
end

function vS1024AddS(a, b, result)
    @ccall libacc.vS1024AddS(a::Ptr{vS1024}, b::Ptr{vS1024}, result::Ptr{vS1024})::Cvoid
end

function vU256Mod(numerator, divisor, remainder)
    @ccall libacc.vU256Mod(numerator::Ptr{vU256}, divisor::Ptr{vU256}, remainder::Ptr{vU256})::Cvoid
end

function vS256Mod(numerator, divisor, remainder)
    @ccall libacc.vS256Mod(numerator::Ptr{vS256}, divisor::Ptr{vS256}, remainder::Ptr{vS256})::Cvoid
end

function vU512Mod(numerator, divisor, remainder)
    @ccall libacc.vU512Mod(numerator::Ptr{vU512}, divisor::Ptr{vU512}, remainder::Ptr{vU512})::Cvoid
end

function vS512Mod(numerator, divisor, remainder)
    @ccall libacc.vS512Mod(numerator::Ptr{vS512}, divisor::Ptr{vS512}, remainder::Ptr{vS512})::Cvoid
end

function vU1024Mod(numerator, divisor, remainder)
    @ccall libacc.vU1024Mod(numerator::Ptr{vU1024}, divisor::Ptr{vU1024}, remainder::Ptr{vU1024})::Cvoid
end

function vS1024Mod(numerator, divisor, remainder)
    @ccall libacc.vS1024Mod(numerator::Ptr{vS1024}, divisor::Ptr{vS1024}, remainder::Ptr{vS1024})::Cvoid
end

function vLL256Shift(a, shiftAmount, result)
    @ccall libacc.vLL256Shift(a::Ptr{vU256}, shiftAmount::UInt32, result::Ptr{vU256})::Cvoid
end

function vLL512Shift(a, shiftAmount, result)
    @ccall libacc.vLL512Shift(a::Ptr{vU512}, shiftAmount::UInt32, result::Ptr{vU512})::Cvoid
end

function vLL1024Shift(a, shiftAmount, result)
    @ccall libacc.vLL1024Shift(a::Ptr{vU1024}, shiftAmount::UInt32, result::Ptr{vU1024})::Cvoid
end

function vLR256Shift(a, shiftAmount, result)
    @ccall libacc.vLR256Shift(a::Ptr{vU256}, shiftAmount::UInt32, result::Ptr{vU256})::Cvoid
end

function vLR512Shift(a, shiftAmount, result)
    @ccall libacc.vLR512Shift(a::Ptr{vU512}, shiftAmount::UInt32, result::Ptr{vU512})::Cvoid
end

function vLR1024Shift(a, shiftAmount, result)
    @ccall libacc.vLR1024Shift(a::Ptr{vU1024}, shiftAmount::UInt32, result::Ptr{vU1024})::Cvoid
end

function vA256Shift(a, shiftAmount, result)
    @ccall libacc.vA256Shift(a::Ptr{vS256}, shiftAmount::UInt32, result::Ptr{vS256})::Cvoid
end

function vA512Shift(a, shiftAmount, result)
    @ccall libacc.vA512Shift(a::Ptr{vS512}, shiftAmount::UInt32, result::Ptr{vS512})::Cvoid
end

function vA1024Shift(a, shiftAmount, result)
    @ccall libacc.vA1024Shift(a::Ptr{vS1024}, shiftAmount::UInt32, result::Ptr{vS1024})::Cvoid
end

function vL256Rotate(a, rotateAmount, result)
    @ccall libacc.vL256Rotate(a::Ptr{vU256}, rotateAmount::UInt32, result::Ptr{vU256})::Cvoid
end

function vL512Rotate(a, rotateAmount, result)
    @ccall libacc.vL512Rotate(a::Ptr{vU512}, rotateAmount::UInt32, result::Ptr{vU512})::Cvoid
end

function vL1024Rotate(a, rotateAmount, result)
    @ccall libacc.vL1024Rotate(a::Ptr{vU1024}, rotateAmount::UInt32, result::Ptr{vU1024})::Cvoid
end

function vR256Rotate(a, rotateAmount, result)
    @ccall libacc.vR256Rotate(a::Ptr{vU256}, rotateAmount::UInt32, result::Ptr{vU256})::Cvoid
end

function vR512Rotate(a, rotateAmount, result)
    @ccall libacc.vR512Rotate(a::Ptr{vU512}, rotateAmount::UInt32, result::Ptr{vU512})::Cvoid
end

function vR1024Rotate(a, rotateAmount, result)
    @ccall libacc.vR1024Rotate(a::Ptr{vU1024}, rotateAmount::UInt32, result::Ptr{vU1024})::Cvoid
end

function _SparseTrap()
    @ccall libacc._SparseTrap()::Cvoid
end

const SparseTriangle_t = Cuchar

const SparseKind_t = Cuint

struct SparseAttributes_t
    data::NTuple{4, UInt8}
end

function Base.getproperty(x::Ptr{SparseAttributes_t}, f::Symbol)
    f === :transpose && return (Ptr{Bool}(x + 0), 0, 1)
    f === :triangle && return (Ptr{SparseTriangle_t}(x + 0), 1, 1)
    f === :kind && return (Ptr{SparseKind_t}(x + 0), 2, 2)
    f === :_reserved && return (Ptr{Cuint}(x + 0), 4, 11)
    f === :_allocatedBySparse && return (Ptr{Bool}(x + 0), 15, 1)
    return getfield(x, f)
end

function Base.getproperty(x::SparseAttributes_t, f::Symbol)
    r = Ref{SparseAttributes_t}(x)
    ptr = Base.unsafe_convert(Ptr{SparseAttributes_t}, r)
    fptr = getproperty(ptr, f)
    begin
        if fptr isa Ptr
            return GC.@preserve(r, unsafe_load(fptr))
        else
            (baseptr, offset, width) = fptr
            ty = eltype(baseptr)
            baseptr32 = convert(Ptr{UInt32}, baseptr)
            u64 = GC.@preserve(r, unsafe_load(baseptr32))
            if offset + width > 32
                u64 |= GC.@preserve(r, unsafe_load(baseptr32 + 4)) << 32
            end
            u64 = u64 >> offset & (1 << width - 1)
            return u64 % ty
        end
    end
end

function Base.setproperty!(x::Ptr{SparseAttributes_t}, f::Symbol, v)
    fptr = getproperty(x, f)
    if fptr isa Ptr
        unsafe_store!(getproperty(x, f), v)
    else
        (baseptr, offset, width) = fptr
        baseptr32 = convert(Ptr{UInt32}, baseptr)
        u64 = unsafe_load(baseptr32)
        straddle = offset + width > 32
        if straddle
            u64 |= unsafe_load(baseptr32 + 4) << 32
        end
        mask = 1 << width - 1
        u64 &= ~(mask << offset)
        u64 |= (unsigned(v) & mask) << offset
        unsafe_store!(baseptr32, u64 & typemax(UInt32))
        if straddle
            unsafe_store!(baseptr32 + 4, u64 >> 32)
        end
    end
end

function Base.propertynames(x::SparseAttributes_t, private::Bool = false)
    (:transpose, :triangle, :kind, :_reserved, :_allocatedBySparse, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

const __SPARSE_double_complex = ComplexF64

struct SparseAttributesComplex_t
    data::NTuple{4, UInt8}
end

function Base.getproperty(x::Ptr{SparseAttributesComplex_t}, f::Symbol)
    f === :transpose && return (Ptr{Bool}(x + 0), 0, 1)
    f === :triangle && return (Ptr{SparseTriangle_t}(x + 0), 1, 1)
    f === :kind && return (Ptr{SparseKind_t}(x + 0), 2, 3)
    f === :conjugate_transpose && return (Ptr{Bool}(x + 0), 5, 1)
    f === :_reserved && return (Ptr{Cuint}(x + 0), 6, 9)
    f === :_allocatedBySparse && return (Ptr{Bool}(x + 0), 15, 1)
    return getfield(x, f)
end

function Base.getproperty(x::SparseAttributesComplex_t, f::Symbol)
    r = Ref{SparseAttributesComplex_t}(x)
    ptr = Base.unsafe_convert(Ptr{SparseAttributesComplex_t}, r)
    fptr = getproperty(ptr, f)
    begin
        if fptr isa Ptr
            return GC.@preserve(r, unsafe_load(fptr))
        else
            (baseptr, offset, width) = fptr
            ty = eltype(baseptr)
            baseptr32 = convert(Ptr{UInt32}, baseptr)
            u64 = GC.@preserve(r, unsafe_load(baseptr32))
            if offset + width > 32
                u64 |= GC.@preserve(r, unsafe_load(baseptr32 + 4)) << 32
            end
            u64 = u64 >> offset & (1 << width - 1)
            return u64 % ty
        end
    end
end

function Base.setproperty!(x::Ptr{SparseAttributesComplex_t}, f::Symbol, v)
    fptr = getproperty(x, f)
    if fptr isa Ptr
        unsafe_store!(getproperty(x, f), v)
    else
        (baseptr, offset, width) = fptr
        baseptr32 = convert(Ptr{UInt32}, baseptr)
        u64 = unsafe_load(baseptr32)
        straddle = offset + width > 32
        if straddle
            u64 |= unsafe_load(baseptr32 + 4) << 32
        end
        mask = 1 << width - 1
        u64 &= ~(mask << offset)
        u64 |= (unsigned(v) & mask) << offset
        unsafe_store!(baseptr32, u64 & typemax(UInt32))
        if straddle
            unsafe_store!(baseptr32 + 4, u64 >> 32)
        end
    end
end

function Base.propertynames(x::SparseAttributesComplex_t, private::Bool = false)
    (:transpose, :triangle, :kind, :conjugate_transpose, :_reserved, :_allocatedBySparse, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

const __SPARSE_float_complex = ComplexF32

@enum CBLAS_ORDER::UInt32 begin
    CblasRowMajor = 101
    CblasColMajor = 102
end

@enum CBLAS_TRANSPOSE::UInt32 begin
    CblasNoTrans = 111
    CblasTrans = 112
    CblasConjTrans = 113
    AtlasConj = 114
end

@enum CBLAS_UPLO::UInt32 begin
    CblasUpper = 121
    CblasLower = 122
end

@enum CBLAS_DIAG::UInt32 begin
    CblasNonUnit = 131
    CblasUnit = 132
end

@enum CBLAS_SIDE::UInt32 begin
    CblasLeft = 141
    CblasRight = 142
end

# typedef void ( * BLASParamErrorProc ) ( const char * funcName , const char * paramName , const int * paramPos , const int * paramValue )
const BLASParamErrorProc = Ptr{Cvoid}

function SetBLASParamErrorProc(__ErrorProc)
    @ccall libacc.SetBLASParamErrorProc(__ErrorProc::BLASParamErrorProc)::Cvoid
end

@enum var"##Ctag#282"::UInt32 begin
    SparseOrdinary = 0
    SparseTriangular = 1
    SparseUnitTriangular = 2
    SparseSymmetric = 3
    SparseHermitian = 7
end

@enum SparseTriangle::UInt32 begin
    SparseUpperTriangle = 0
    SparseLowerTriangle = 1
end

struct SparseMatrixStructure
    data::NTuple{32, UInt8}
end

function Base.getproperty(x::Ptr{SparseMatrixStructure}, f::Symbol)
    f === :rowCount && return Ptr{Cint}(x + 0)
    f === :columnCount && return Ptr{Cint}(x + 4)
    f === :columnStarts && return Ptr{Ptr{Clong}}(x + 8)
    f === :rowIndices && return Ptr{Ptr{Cint}}(x + 16)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 24)
    f === :blockSize && return Ptr{UInt8}(x + 28)
    return getfield(x, f)
end

function Base.getproperty(x::SparseMatrixStructure, f::Symbol)
    r = Ref{SparseMatrixStructure}(x)
    ptr = Base.unsafe_convert(Ptr{SparseMatrixStructure}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseMatrixStructure}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseMatrixStructure, private::Bool = false)
    (:rowCount, :columnCount, :columnStarts, :rowIndices, :attributes, :blockSize, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseMatrixStructureComplex
    data::NTuple{32, UInt8}
end

function Base.getproperty(x::Ptr{SparseMatrixStructureComplex}, f::Symbol)
    f === :rowCount && return Ptr{Cint}(x + 0)
    f === :columnCount && return Ptr{Cint}(x + 4)
    f === :columnStarts && return Ptr{Ptr{Clong}}(x + 8)
    f === :rowIndices && return Ptr{Ptr{Cint}}(x + 16)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 24)
    f === :blockSize && return Ptr{UInt8}(x + 28)
    return getfield(x, f)
end

function Base.getproperty(x::SparseMatrixStructureComplex, f::Symbol)
    r = Ref{SparseMatrixStructureComplex}(x)
    ptr = Base.unsafe_convert(Ptr{SparseMatrixStructureComplex}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseMatrixStructureComplex}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseMatrixStructureComplex, private::Bool = false)
    (:rowCount, :columnCount, :columnStarts, :rowIndices, :attributes, :blockSize, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseMatrix_Double
    data::NTuple{40, UInt8}
end

function Base.getproperty(x::Ptr{SparseMatrix_Double}, f::Symbol)
    f === :structure && return Ptr{SparseMatrixStructure}(x + 0)
    f === :data && return Ptr{Ptr{Cdouble}}(x + 32)
    return getfield(x, f)
end

function Base.getproperty(x::SparseMatrix_Double, f::Symbol)
    r = Ref{SparseMatrix_Double}(x)
    ptr = Base.unsafe_convert(Ptr{SparseMatrix_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseMatrix_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseMatrix_Double, private::Bool = false)
    (:structure, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseMatrix_Float
    data::NTuple{40, UInt8}
end

function Base.getproperty(x::Ptr{SparseMatrix_Float}, f::Symbol)
    f === :structure && return Ptr{SparseMatrixStructure}(x + 0)
    f === :data && return Ptr{Ptr{Cfloat}}(x + 32)
    return getfield(x, f)
end

function Base.getproperty(x::SparseMatrix_Float, f::Symbol)
    r = Ref{SparseMatrix_Float}(x)
    ptr = Base.unsafe_convert(Ptr{SparseMatrix_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseMatrix_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseMatrix_Float, private::Bool = false)
    (:structure, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseMatrix_Complex_Double
    data::NTuple{40, UInt8}
end

function Base.getproperty(x::Ptr{SparseMatrix_Complex_Double}, f::Symbol)
    f === :structure && return Ptr{SparseMatrixStructureComplex}(x + 0)
    f === :data && return Ptr{Ptr{__SPARSE_double_complex}}(x + 32)
    return getfield(x, f)
end

function Base.getproperty(x::SparseMatrix_Complex_Double, f::Symbol)
    r = Ref{SparseMatrix_Complex_Double}(x)
    ptr = Base.unsafe_convert(Ptr{SparseMatrix_Complex_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseMatrix_Complex_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseMatrix_Complex_Double, private::Bool = false)
    (:structure, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseMatrix_Complex_Float
    data::NTuple{40, UInt8}
end

function Base.getproperty(x::Ptr{SparseMatrix_Complex_Float}, f::Symbol)
    f === :structure && return Ptr{SparseMatrixStructureComplex}(x + 0)
    f === :data && return Ptr{Ptr{__SPARSE_float_complex}}(x + 32)
    return getfield(x, f)
end

function Base.getproperty(x::SparseMatrix_Complex_Float, f::Symbol)
    r = Ref{SparseMatrix_Complex_Float}(x)
    ptr = Base.unsafe_convert(Ptr{SparseMatrix_Complex_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseMatrix_Complex_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseMatrix_Complex_Float, private::Bool = false)
    (:structure, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct DenseVector_Double
    count::Cint
    data::Ptr{Cdouble}
end

struct DenseVector_Float
    count::Cint
    data::Ptr{Cfloat}
end

struct DenseVector_Complex_Double
    count::Cint
    data::Ptr{__SPARSE_double_complex}
end

struct DenseVector_Complex_Float
    count::Cint
    data::Ptr{__SPARSE_float_complex}
end

struct DenseMatrix_Double
    data::NTuple{24, UInt8}
end

function Base.getproperty(x::Ptr{DenseMatrix_Double}, f::Symbol)
    f === :rowCount && return Ptr{Cint}(x + 0)
    f === :columnCount && return Ptr{Cint}(x + 4)
    f === :columnStride && return Ptr{Cint}(x + 8)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 12)
    f === :data && return Ptr{Ptr{Cdouble}}(x + 16)
    return getfield(x, f)
end

function Base.getproperty(x::DenseMatrix_Double, f::Symbol)
    r = Ref{DenseMatrix_Double}(x)
    ptr = Base.unsafe_convert(Ptr{DenseMatrix_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{DenseMatrix_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::DenseMatrix_Double, private::Bool = false)
    (:rowCount, :columnCount, :columnStride, :attributes, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct DenseMatrix_Float
    data::NTuple{24, UInt8}
end

function Base.getproperty(x::Ptr{DenseMatrix_Float}, f::Symbol)
    f === :rowCount && return Ptr{Cint}(x + 0)
    f === :columnCount && return Ptr{Cint}(x + 4)
    f === :columnStride && return Ptr{Cint}(x + 8)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 12)
    f === :data && return Ptr{Ptr{Cfloat}}(x + 16)
    return getfield(x, f)
end

function Base.getproperty(x::DenseMatrix_Float, f::Symbol)
    r = Ref{DenseMatrix_Float}(x)
    ptr = Base.unsafe_convert(Ptr{DenseMatrix_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{DenseMatrix_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::DenseMatrix_Float, private::Bool = false)
    (:rowCount, :columnCount, :columnStride, :attributes, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct DenseMatrix_Complex_Double
    data::NTuple{24, UInt8}
end

function Base.getproperty(x::Ptr{DenseMatrix_Complex_Double}, f::Symbol)
    f === :rowCount && return Ptr{Cint}(x + 0)
    f === :columnCount && return Ptr{Cint}(x + 4)
    f === :columnStride && return Ptr{Cint}(x + 8)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 12)
    f === :data && return Ptr{Ptr{__SPARSE_double_complex}}(x + 16)
    return getfield(x, f)
end

function Base.getproperty(x::DenseMatrix_Complex_Double, f::Symbol)
    r = Ref{DenseMatrix_Complex_Double}(x)
    ptr = Base.unsafe_convert(Ptr{DenseMatrix_Complex_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{DenseMatrix_Complex_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::DenseMatrix_Complex_Double, private::Bool = false)
    (:rowCount, :columnCount, :columnStride, :attributes, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct DenseMatrix_Complex_Float
    data::NTuple{24, UInt8}
end

function Base.getproperty(x::Ptr{DenseMatrix_Complex_Float}, f::Symbol)
    f === :rowCount && return Ptr{Cint}(x + 0)
    f === :columnCount && return Ptr{Cint}(x + 4)
    f === :columnStride && return Ptr{Cint}(x + 8)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 12)
    f === :data && return Ptr{Ptr{__SPARSE_float_complex}}(x + 16)
    return getfield(x, f)
end

function Base.getproperty(x::DenseMatrix_Complex_Float, f::Symbol)
    r = Ref{DenseMatrix_Complex_Float}(x)
    ptr = Base.unsafe_convert(Ptr{DenseMatrix_Complex_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{DenseMatrix_Complex_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::DenseMatrix_Complex_Float, private::Bool = false)
    (:rowCount, :columnCount, :columnStride, :attributes, :data, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

const SparseStatus_t = Cint

@enum var"##Ctag#283"::Int32 begin
    SparseStatusOK = 0
    SparseFactorizationFailed = -1
    SparseMatrixIsSingular = -2
    SparseInternalError = -3
    SparseParameterError = -4
    SparseStatusReleased = -2147483647
end

const SparseFactorization_t = UInt8

@enum var"##Ctag#284"::UInt32 begin
    SparseFactorizationCholesky = 0
    SparseFactorizationLDLT = 1
    SparseFactorizationLDLTUnpivoted = 2
    SparseFactorizationLDLTSBK = 3
    SparseFactorizationLDLTTPP = 4
    SparseFactorizationQR = 40
    SparseFactorizationCholeskyAtA = 41
    SparseFactorizationLU = 80
    SparseFactorizationLUUnpivoted = 81
    SparseFactorizationLUSPP = 82
    SparseFactorizationLUTPP = 83
end

const SparseControl_t = UInt32

@enum var"##Ctag#285"::UInt32 begin
    SparseDefaultControl = 0
end

const SparseOrder_t = UInt8

@enum var"##Ctag#286"::UInt32 begin
    SparseOrderDefault = 0
    SparseOrderUser = 1
    SparseOrderAMD = 2
    SparseOrderMetis = 3
    SparseOrderCOLAMD = 4
    SparseOrderMTMetis = 5
end

const SparseScaling_t = UInt8

@enum var"##Ctag#287"::UInt32 begin
    SparseScalingDefault = 0
    SparseScalingUser = 1
    SparseScalingEquilibriationInf = 2
    SparseScalingHungarianScalingOnly = 3
    SparseScalingHungarianScalingAndOrdering = 4
end

struct SparseSymbolicFactorOptions
    control::SparseControl_t
    orderMethod::SparseOrder_t
    order::Ptr{Cint}
    ignoreRowsAndColumns::Ptr{Cint}
    malloc::Ptr{Cvoid}
    free::Ptr{Cvoid}
    reportError::Ptr{Cvoid}
end

struct SparseNumericFactorOptions
    control::SparseControl_t
    scalingMethod::SparseScaling_t
    scaling::Ptr{Cvoid}
    pivotTolerance::Cdouble
    zeroTolerance::Cdouble
end

struct SparseOpaqueSymbolicFactorization
    data::NTuple{64, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueSymbolicFactorization}, f::Symbol)
    f === :status && return Ptr{SparseStatus_t}(x + 0)
    f === :rowCount && return Ptr{Cint}(x + 4)
    f === :columnCount && return Ptr{Cint}(x + 8)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 12)
    f === :blockSize && return Ptr{UInt8}(x + 16)
    f === :type && return Ptr{SparseFactorization_t}(x + 17)
    f === :factorization && return Ptr{Ptr{Cvoid}}(x + 24)
    f === :workspaceSize_Float && return Ptr{Csize_t}(x + 32)
    f === :workspaceSize_Double && return Ptr{Csize_t}(x + 40)
    f === :factorSize_Float && return Ptr{Csize_t}(x + 48)
    f === :factorSize_Double && return Ptr{Csize_t}(x + 56)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueSymbolicFactorization, f::Symbol)
    r = Ref{SparseOpaqueSymbolicFactorization}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueSymbolicFactorization}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueSymbolicFactorization}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueSymbolicFactorization, private::Bool = false)
    (:status, :rowCount, :columnCount, :attributes, :blockSize, :type, :factorization, :workspaceSize_Float, :workspaceSize_Double, :factorSize_Float, :factorSize_Double, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueFactorization_Double
    data::NTuple{104, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueFactorization_Double}, f::Symbol)
    f === :status && return Ptr{SparseStatus_t}(x + 0)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 4)
    f === :symbolicFactorization && return Ptr{SparseOpaqueSymbolicFactorization}(x + 8)
    f === :userFactorStorage && return Ptr{Bool}(x + 72)
    f === :numericFactorization && return Ptr{Ptr{Cvoid}}(x + 80)
    f === :solveWorkspaceRequiredStatic && return Ptr{Csize_t}(x + 88)
    f === :solveWorkspaceRequiredPerRHS && return Ptr{Csize_t}(x + 96)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueFactorization_Double, f::Symbol)
    r = Ref{SparseOpaqueFactorization_Double}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueFactorization_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueFactorization_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueFactorization_Double, private::Bool = false)
    (:status, :attributes, :symbolicFactorization, :userFactorStorage, :numericFactorization, :solveWorkspaceRequiredStatic, :solveWorkspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueFactorization_Float
    data::NTuple{104, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueFactorization_Float}, f::Symbol)
    f === :status && return Ptr{SparseStatus_t}(x + 0)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 4)
    f === :symbolicFactorization && return Ptr{SparseOpaqueSymbolicFactorization}(x + 8)
    f === :userFactorStorage && return Ptr{Bool}(x + 72)
    f === :numericFactorization && return Ptr{Ptr{Cvoid}}(x + 80)
    f === :solveWorkspaceRequiredStatic && return Ptr{Csize_t}(x + 88)
    f === :solveWorkspaceRequiredPerRHS && return Ptr{Csize_t}(x + 96)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueFactorization_Float, f::Symbol)
    r = Ref{SparseOpaqueFactorization_Float}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueFactorization_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueFactorization_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueFactorization_Float, private::Bool = false)
    (:status, :attributes, :symbolicFactorization, :userFactorStorage, :numericFactorization, :solveWorkspaceRequiredStatic, :solveWorkspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueFactorization_Complex_Double
    data::NTuple{104, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueFactorization_Complex_Double}, f::Symbol)
    f === :status && return Ptr{SparseStatus_t}(x + 0)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 4)
    f === :symbolicFactorization && return Ptr{SparseOpaqueSymbolicFactorization}(x + 8)
    f === :userFactorStorage && return Ptr{Bool}(x + 72)
    f === :numericFactorization && return Ptr{Ptr{Cvoid}}(x + 80)
    f === :solveWorkspaceRequiredStatic && return Ptr{Csize_t}(x + 88)
    f === :solveWorkspaceRequiredPerRHS && return Ptr{Csize_t}(x + 96)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueFactorization_Complex_Double, f::Symbol)
    r = Ref{SparseOpaqueFactorization_Complex_Double}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueFactorization_Complex_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueFactorization_Complex_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueFactorization_Complex_Double, private::Bool = false)
    (:status, :attributes, :symbolicFactorization, :userFactorStorage, :numericFactorization, :solveWorkspaceRequiredStatic, :solveWorkspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueFactorization_Complex_Float
    data::NTuple{104, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueFactorization_Complex_Float}, f::Symbol)
    f === :status && return Ptr{SparseStatus_t}(x + 0)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 4)
    f === :symbolicFactorization && return Ptr{SparseOpaqueSymbolicFactorization}(x + 8)
    f === :userFactorStorage && return Ptr{Bool}(x + 72)
    f === :numericFactorization && return Ptr{Ptr{Cvoid}}(x + 80)
    f === :solveWorkspaceRequiredStatic && return Ptr{Csize_t}(x + 88)
    f === :solveWorkspaceRequiredPerRHS && return Ptr{Csize_t}(x + 96)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueFactorization_Complex_Float, f::Symbol)
    r = Ref{SparseOpaqueFactorization_Complex_Float}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueFactorization_Complex_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueFactorization_Complex_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueFactorization_Complex_Float, private::Bool = false)
    (:status, :attributes, :symbolicFactorization, :userFactorStorage, :numericFactorization, :solveWorkspaceRequiredStatic, :solveWorkspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

const SparseSubfactor_t = UInt8

@enum var"##Ctag#288"::UInt32 begin
    SparseSubfactorInvalid = 0
    SparseSubfactorP = 1
    SparseSubfactorS = 2
    SparseSubfactorL = 3
    SparseSubfactorD = 4
    SparseSubfactorPLPS = 5
    SparseSubfactorQ = 6
    SparseSubfactorR = 7
    SparseSubfactorRP = 8
    SparseSubfactorSr = 9
    SparseSubfactorSc = 10
end

struct SparseOpaqueSubfactor_Double
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueSubfactor_Double}, f::Symbol)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 0)
    f === :contents && return Ptr{SparseSubfactor_t}(x + 4)
    f === :factor && return Ptr{SparseOpaqueFactorization_Double}(x + 8)
    f === :workspaceRequiredStatic && return Ptr{Csize_t}(x + 112)
    f === :workspaceRequiredPerRHS && return Ptr{Csize_t}(x + 120)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueSubfactor_Double, f::Symbol)
    r = Ref{SparseOpaqueSubfactor_Double}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueSubfactor_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueSubfactor_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueSubfactor_Double, private::Bool = false)
    (:attributes, :contents, :factor, :workspaceRequiredStatic, :workspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueSubfactor_Float
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueSubfactor_Float}, f::Symbol)
    f === :attributes && return Ptr{SparseAttributes_t}(x + 0)
    f === :contents && return Ptr{SparseSubfactor_t}(x + 4)
    f === :factor && return Ptr{SparseOpaqueFactorization_Float}(x + 8)
    f === :workspaceRequiredStatic && return Ptr{Csize_t}(x + 112)
    f === :workspaceRequiredPerRHS && return Ptr{Csize_t}(x + 120)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueSubfactor_Float, f::Symbol)
    r = Ref{SparseOpaqueSubfactor_Float}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueSubfactor_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueSubfactor_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueSubfactor_Float, private::Bool = false)
    (:attributes, :contents, :factor, :workspaceRequiredStatic, :workspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueSubfactor_Complex_Double
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueSubfactor_Complex_Double}, f::Symbol)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 0)
    f === :contents && return Ptr{SparseSubfactor_t}(x + 4)
    f === :factor && return Ptr{SparseOpaqueFactorization_Complex_Double}(x + 8)
    f === :workspaceRequiredStatic && return Ptr{Csize_t}(x + 112)
    f === :workspaceRequiredPerRHS && return Ptr{Csize_t}(x + 120)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueSubfactor_Complex_Double, f::Symbol)
    r = Ref{SparseOpaqueSubfactor_Complex_Double}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueSubfactor_Complex_Double}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueSubfactor_Complex_Double}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueSubfactor_Complex_Double, private::Bool = false)
    (:attributes, :contents, :factor, :workspaceRequiredStatic, :workspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseOpaqueSubfactor_Complex_Float
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{SparseOpaqueSubfactor_Complex_Float}, f::Symbol)
    f === :attributes && return Ptr{SparseAttributesComplex_t}(x + 0)
    f === :contents && return Ptr{SparseSubfactor_t}(x + 4)
    f === :factor && return Ptr{SparseOpaqueFactorization_Complex_Float}(x + 8)
    f === :workspaceRequiredStatic && return Ptr{Csize_t}(x + 112)
    f === :workspaceRequiredPerRHS && return Ptr{Csize_t}(x + 120)
    return getfield(x, f)
end

function Base.getproperty(x::SparseOpaqueSubfactor_Complex_Float, f::Symbol)
    r = Ref{SparseOpaqueSubfactor_Complex_Float}(x)
    ptr = Base.unsafe_convert(Ptr{SparseOpaqueSubfactor_Complex_Float}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseOpaqueSubfactor_Complex_Float}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseOpaqueSubfactor_Complex_Float, private::Bool = false)
    (:attributes, :contents, :factor, :workspaceRequiredStatic, :workspaceRequiredPerRHS, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

const SparseUpdate_t = UInt8

@enum var"##Ctag#289"::UInt32 begin
    SparseUpdatePartialRefactor = 0
end

const SparsePreconditioner_t = Cint

@enum var"##Ctag#290"::UInt32 begin
    SparsePreconditionerNone = 0
    SparsePreconditionerUser = 1
    SparsePreconditionerDiagonal = 2
    SparsePreconditionerDiagScaling = 3
end

struct SparseOpaquePreconditioner_Double
    type::SparsePreconditioner_t
    mem::Ptr{Cvoid}
    apply::Ptr{Cvoid}
end

struct SparseOpaquePreconditioner_Float
    type::SparsePreconditioner_t
    mem::Ptr{Cvoid}
    apply::Ptr{Cvoid}
end

struct SparseOpaquePreconditioner_Complex_Double
    type::SparsePreconditioner_t
    mem::Ptr{Cvoid}
    apply::Ptr{Cvoid}
end

struct SparseOpaquePreconditioner_Complex_Float
    type::SparsePreconditioner_t
    mem::Ptr{Cvoid}
    apply::Ptr{Cvoid}
end

const SparseIterativeStatus_t = Cint

@enum var"##Ctag#291"::Int32 begin
    SparseIterativeConverged = 0
    SparseIterativeMaxIterations = 1
    SparseIterativeParameterError = -1
    SparseIterativeIllConditioned = -2
    SparseIterativeInternalError = -99
end

struct _SparseIterativeMethodBaseOptions
    reportError::Ptr{Cvoid}
end

struct SparseCGOptions
    reportError::Ptr{Cvoid}
    maxIterations::Cint
    atol::Cdouble
    rtol::Cdouble
    reportStatus::Ptr{Cvoid}
end

const SparseGMRESVariant_t = UInt8

@enum var"##Ctag#292"::UInt32 begin
    SparseVariantDQGMRES = 0
    SparseVariantGMRES = 1
    SparseVariantFGMRES = 2
end

struct SparseGMRESOptions
    reportError::Ptr{Cvoid}
    variant::SparseGMRESVariant_t
    nvec::Cint
    maxIterations::Cint
    atol::Cdouble
    rtol::Cdouble
    reportStatus::Ptr{Cvoid}
end

const SparseLSMRConvergenceTest_t = Cint

@enum var"##Ctag#293"::UInt32 begin
    SparseLSMRCTDefault = 0
    SparseLSMRCTFongSaunders = 1
end

struct SparseLSMROptions
    reportError::Ptr{Cvoid}
    lambda::Cdouble
    nvec::Cint
    convergenceTest::SparseLSMRConvergenceTest_t
    atol::Cdouble
    rtol::Cdouble
    btol::Cdouble
    conditionLimit::Cdouble
    maxIterations::Cint
    reportStatus::Ptr{Cvoid}
end

struct var"##Ctag#349"
    data::NTuple{256, UInt8}
end

function Base.getproperty(x::Ptr{var"##Ctag#349"}, f::Symbol)
    f === :base && return Ptr{_SparseIterativeMethodBaseOptions}(x + 0)
    f === :cg && return Ptr{SparseCGOptions}(x + 0)
    f === :gmres && return Ptr{SparseGMRESOptions}(x + 0)
    f === :lsmr && return Ptr{SparseLSMROptions}(x + 0)
    f === :padding && return Ptr{NTuple{256, Cchar}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#349", f::Symbol)
    r = Ref{var"##Ctag#349"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#349"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#349"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::var"##Ctag#349", private::Bool = false)
    (:base, :cg, :gmres, :lsmr, :padding, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct SparseIterativeMethod
    data::NTuple{264, UInt8}
end

function Base.getproperty(x::Ptr{SparseIterativeMethod}, f::Symbol)
    f === :method && return Ptr{Cint}(x + 0)
    f === :options && return Ptr{var"##Ctag#349"}(x + 8)
    return getfield(x, f)
end

function Base.getproperty(x::SparseIterativeMethod, f::Symbol)
    r = Ref{SparseIterativeMethod}(x)
    ptr = Base.unsafe_convert(Ptr{SparseIterativeMethod}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{SparseIterativeMethod}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::SparseIterativeMethod, private::Bool = false)
    (:method, :options, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

const _SparseIterativeMethod_t = Cint

@enum var"##Ctag#298"::UInt32 begin
    _SparseMethodCG = 0
    _SparseMethodGMRES = 1
    _SparseMethodLSMR = 2
end

function _SparseSymbolicFactorSymmetric(factorType, Matrix, options)
    @ccall libacc._SparseSymbolicFactorSymmetric(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrixStructure}, options::Ptr{SparseSymbolicFactorOptions})::SparseOpaqueSymbolicFactorization
end

function _SparseSymbolicFactorQR(factorType, Matrix, options)
    @ccall libacc._SparseSymbolicFactorQR(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrixStructure}, options::Ptr{SparseSymbolicFactorOptions})::SparseOpaqueSymbolicFactorization
end

function _SparseSymbolicFactorLU(factorType, Matrix, options)
    @ccall libacc._SparseSymbolicFactorLU(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrixStructure}, options::Ptr{SparseSymbolicFactorOptions})::SparseOpaqueSymbolicFactorization
end

function _SparseRetainSymbolic(symbolicFactor)
    @ccall libacc._SparseRetainSymbolic(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization})::Cvoid
end

function _SparseDestroyOpaqueSymbolic(toFree)
    @ccall libacc._SparseDestroyOpaqueSymbolic(toFree::Ptr{SparseOpaqueSymbolicFactorization})::Cvoid
end

function _SparseGetOptionsFromSymbolicFactor(factor)
    @ccall libacc._SparseGetOptionsFromSymbolicFactor(factor::Ptr{SparseOpaqueSymbolicFactorization})::SparseSymbolicFactorOptions
end

function _SparseFromKindComplex(K)
    @ccall libacc._SparseFromKindComplex(K::SparseKind_t)::SparseKind_t
end

function _SparseToKindComplex(K)
    @ccall libacc._SparseToKindComplex(K::SparseKind_t)::SparseKind_t
end

function _SparseFromAttributeComplex(K)
    @ccall libacc._SparseFromAttributeComplex(K::SparseAttributesComplex_t)::SparseAttributes_t
end

function _SparseToAttributeComplex(K)
    @ccall libacc._SparseToAttributeComplex(K::SparseAttributes_t)::SparseAttributesComplex_t
end

function _SparseFromStructureComplex(K)
    @ccall libacc._SparseFromStructureComplex(K::SparseMatrixStructureComplex)::SparseMatrixStructure
end

function _SparseToStructureComplex(K)
    @ccall libacc._SparseToStructureComplex(K::SparseMatrixStructure)::SparseMatrixStructureComplex
end

function _SparseConvertFromCoordinate_Double(m, n, nBlock, blockSize, attributes, row, col, val, storage, workspace)
    @ccall libacc._SparseConvertFromCoordinate_Double(m::Cint, n::Cint, nBlock::Clong, blockSize::UInt8, attributes::SparseAttributes_t, row::Ptr{Cint}, col::Ptr{Cint}, val::Ptr{Cdouble}, storage::Ptr{Cchar}, workspace::Ptr{Cint})::SparseMatrix_Double
end

function _SparseNumericFactorSymmetric_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorSymmetric_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Double
end

function _SparseNumericFactorQR_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorQR_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Double
end

function _SparseNumericFactorLU_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorLU_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Double
end

function _SparseFactorSymmetric_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorSymmetric_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Double
end

function _SparseFactorQR_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorQR_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Double
end

function _SparseFactorLU_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorLU_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Double
end

function _SparseRefactorSymmetric_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorSymmetric_Double(Matrix::Ptr{SparseMatrix_Double}, Factorization::Ptr{SparseOpaqueFactorization_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorQR_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorQR_Double(Matrix::Ptr{SparseMatrix_Double}, Factorization::Ptr{SparseOpaqueFactorization_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorLU_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorLU_Double(Matrix::Ptr{SparseMatrix_Double}, Factorization::Ptr{SparseOpaqueFactorization_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseUpdatePartialRefactorLU_Double(Opaque, updateCount, updatedIndices, newMatrix)
    @ccall libacc._SparseUpdatePartialRefactorLU_Double(Opaque::Ptr{SparseOpaqueFactorization_Double}, updateCount::Cint, updatedIndices::Ptr{Cint}, newMatrix::SparseMatrix_Double)::Cvoid
end

function _SparseMultiplySubfactor_Double(Subfactor, x, y, workspace)
    @ccall libacc._SparseMultiplySubfactor_Double(Subfactor::Ptr{SparseOpaqueSubfactor_Double}, x::Ptr{DenseMatrix_Double}, y::Ptr{DenseMatrix_Double}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveSubfactor_Double(Subfactor, b, x, workspace)
    @ccall libacc._SparseSolveSubfactor_Double(Subfactor::Ptr{SparseOpaqueSubfactor_Double}, b::Ptr{DenseMatrix_Double}, x::Ptr{DenseMatrix_Double}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveOpaque_Double(Factored, RHS, Soln, workspace)
    @ccall libacc._SparseSolveOpaque_Double(Factored::Ptr{SparseOpaqueFactorization_Double}, RHS::Ptr{DenseMatrix_Double}, Soln::Ptr{DenseMatrix_Double}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseDestroyOpaqueNumeric_Double(toFree)
    @ccall libacc._SparseDestroyOpaqueNumeric_Double(toFree::Ptr{SparseOpaqueFactorization_Double})::Cvoid
end

function _SparseRetainNumeric_Double(numericFactor)
    @ccall libacc._SparseRetainNumeric_Double(numericFactor::Ptr{SparseOpaqueFactorization_Double})::Cvoid
end

function _SparseGetOptionsFromNumericFactor_Double(factor)
    @ccall libacc._SparseGetOptionsFromNumericFactor_Double(factor::Ptr{SparseOpaqueFactorization_Double})::SparseNumericFactorOptions
end

function _SparseGetWorkspaceRequired_Double(Subfactor, Factor, workStatic, workPerRHS)
    @ccall libacc._SparseGetWorkspaceRequired_Double(Subfactor::SparseSubfactor_t, Factor::SparseOpaqueFactorization_Double, workStatic::Ptr{Csize_t}, workPerRHS::Ptr{Csize_t})::Cvoid
end

function _SparseGetIterativeStateSize_Double(method, preconditioner, m, n, nrhs)
    @ccall libacc._SparseGetIterativeStateSize_Double(method::Ptr{SparseIterativeMethod}, preconditioner::Bool, m::Cint, n::Cint, nrhs::Cint)::Csize_t
end

function _SparseCreatePreconditioner_Double(type, A)
    @ccall libacc._SparseCreatePreconditioner_Double(type::SparsePreconditioner_t, A::Ptr{SparseMatrix_Double})::SparseOpaquePreconditioner_Double
end

function _SparseReleaseOpaquePreconditioner_Double(toFree)
    @ccall libacc._SparseReleaseOpaquePreconditioner_Double(toFree::Ptr{SparseOpaquePreconditioner_Double})::Cvoid
end

function _SparseSpMV_Double(alpha, A, x, accumulate, y)
    @ccall libacc._SparseSpMV_Double(alpha::Cdouble, A::SparseMatrix_Double, x::DenseMatrix_Double, accumulate::Bool, y::DenseMatrix_Double)::Cvoid
end

function _SparseConvertFromCoordinate_Float(m, n, nBlock, blockSize, attributes, row, col, val, storage, workspace)
    @ccall libacc._SparseConvertFromCoordinate_Float(m::Cint, n::Cint, nBlock::Clong, blockSize::UInt8, attributes::SparseAttributes_t, row::Ptr{Cint}, col::Ptr{Cint}, val::Ptr{Cfloat}, storage::Ptr{Cchar}, workspace::Ptr{Cint})::SparseMatrix_Float
end

function _SparseNumericFactorSymmetric_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorSymmetric_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Float
end

function _SparseNumericFactorQR_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorQR_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Float
end

function _SparseNumericFactorLU_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorLU_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Float
end

function _SparseFactorSymmetric_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorSymmetric_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Float
end

function _SparseFactorQR_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorQR_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Float
end

function _SparseFactorLU_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorLU_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Float
end

function _SparseRefactorSymmetric_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorSymmetric_Float(Matrix::Ptr{SparseMatrix_Float}, Factorization::Ptr{SparseOpaqueFactorization_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorQR_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorQR_Float(Matrix::Ptr{SparseMatrix_Float}, Factorization::Ptr{SparseOpaqueFactorization_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorLU_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorLU_Float(Matrix::Ptr{SparseMatrix_Float}, Factorization::Ptr{SparseOpaqueFactorization_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseUpdatePartialRefactorLU_Float(Opaque, updateCount, updatedIndices, newMatrix)
    @ccall libacc._SparseUpdatePartialRefactorLU_Float(Opaque::Ptr{SparseOpaqueFactorization_Float}, updateCount::Cint, updatedIndices::Ptr{Cint}, newMatrix::SparseMatrix_Float)::Cvoid
end

function _SparseMultiplySubfactor_Float(Subfactor, x, y, workspace)
    @ccall libacc._SparseMultiplySubfactor_Float(Subfactor::Ptr{SparseOpaqueSubfactor_Float}, x::Ptr{DenseMatrix_Float}, y::Ptr{DenseMatrix_Float}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveSubfactor_Float(Subfactor, b, x, workspace)
    @ccall libacc._SparseSolveSubfactor_Float(Subfactor::Ptr{SparseOpaqueSubfactor_Float}, b::Ptr{DenseMatrix_Float}, x::Ptr{DenseMatrix_Float}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveOpaque_Float(Factored, RHS, Soln, workspace)
    @ccall libacc._SparseSolveOpaque_Float(Factored::Ptr{SparseOpaqueFactorization_Float}, RHS::Ptr{DenseMatrix_Float}, Soln::Ptr{DenseMatrix_Float}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseDestroyOpaqueNumeric_Float(toFree)
    @ccall libacc._SparseDestroyOpaqueNumeric_Float(toFree::Ptr{SparseOpaqueFactorization_Float})::Cvoid
end

function _SparseRetainNumeric_Float(numericFactor)
    @ccall libacc._SparseRetainNumeric_Float(numericFactor::Ptr{SparseOpaqueFactorization_Float})::Cvoid
end

function _SparseGetOptionsFromNumericFactor_Float(factor)
    @ccall libacc._SparseGetOptionsFromNumericFactor_Float(factor::Ptr{SparseOpaqueFactorization_Float})::SparseNumericFactorOptions
end

function _SparseGetWorkspaceRequired_Float(Subfactor, Factor, workStatic, workPerRHS)
    @ccall libacc._SparseGetWorkspaceRequired_Float(Subfactor::SparseSubfactor_t, Factor::SparseOpaqueFactorization_Float, workStatic::Ptr{Csize_t}, workPerRHS::Ptr{Csize_t})::Cvoid
end

function _SparseGetIterativeStateSize_Float(method, preconditioner, m, n, nrhs)
    @ccall libacc._SparseGetIterativeStateSize_Float(method::Ptr{SparseIterativeMethod}, preconditioner::Bool, m::Cint, n::Cint, nrhs::Cint)::Csize_t
end

function _SparseCreatePreconditioner_Float(type, A)
    @ccall libacc._SparseCreatePreconditioner_Float(type::SparsePreconditioner_t, A::Ptr{SparseMatrix_Float})::SparseOpaquePreconditioner_Float
end

function _SparseReleaseOpaquePreconditioner_Float(toFree)
    @ccall libacc._SparseReleaseOpaquePreconditioner_Float(toFree::Ptr{SparseOpaquePreconditioner_Float})::Cvoid
end

function _SparseSpMV_Float(alpha, A, x, accumulate, y)
    @ccall libacc._SparseSpMV_Float(alpha::Cfloat, A::SparseMatrix_Float, x::DenseMatrix_Float, accumulate::Bool, y::DenseMatrix_Float)::Cvoid
end

function _SparseConvertFromCoordinate_Complex_Double(m, n, nBlock, blockSize, attributes, row, col, val, storage, workspace)
    @ccall libacc._SparseConvertFromCoordinate_Complex_Double(m::Cint, n::Cint, nBlock::Clong, blockSize::UInt8, attributes::SparseAttributesComplex_t, row::Ptr{Cint}, col::Ptr{Cint}, val::Ptr{__SPARSE_double_complex}, storage::Ptr{Cchar}, workspace::Ptr{Cint})::SparseMatrix_Complex_Double
end

function _SparseNumericFactorSymmetric_Complex_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorSymmetric_Complex_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Double
end

function _SparseNumericFactorHermitian_Complex_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorHermitian_Complex_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Double
end

function _SparseNumericFactorQR_Complex_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorQR_Complex_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Double
end

function _SparseNumericFactorLU_Complex_Double(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorLU_Complex_Double(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Double}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Double
end

function _SparseFactorSymmetric_Complex_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorSymmetric_Complex_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Double
end

function _SparseFactorHermitian_Complex_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorHermitian_Complex_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Double
end

function _SparseFactorQR_Complex_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorQR_Complex_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Double
end

function _SparseFactorLU_Complex_Double(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorLU_Complex_Double(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Double}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Double
end

function _SparseRefactorSymmetric_Complex_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorSymmetric_Complex_Double(Matrix::Ptr{SparseMatrix_Complex_Double}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorHermitian_Complex_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorHermitian_Complex_Double(Matrix::Ptr{SparseMatrix_Complex_Double}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorQR_Complex_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorQR_Complex_Double(Matrix::Ptr{SparseMatrix_Complex_Double}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorLU_Complex_Double(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorLU_Complex_Double(Matrix::Ptr{SparseMatrix_Complex_Double}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Double}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseUpdatePartialRefactorLU_Complex_Double(Opaque, updateCount, updatedIndices, newMatrix)
    @ccall libacc._SparseUpdatePartialRefactorLU_Complex_Double(Opaque::Ptr{SparseOpaqueFactorization_Complex_Double}, updateCount::Cint, updatedIndices::Ptr{Cint}, newMatrix::SparseMatrix_Complex_Double)::Cvoid
end

function _SparseMultiplySubfactor_Complex_Double(Subfactor, x, y, workspace)
    @ccall libacc._SparseMultiplySubfactor_Complex_Double(Subfactor::Ptr{SparseOpaqueSubfactor_Complex_Double}, x::Ptr{DenseMatrix_Complex_Double}, y::Ptr{DenseMatrix_Complex_Double}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveSubfactor_Complex_Double(Subfactor, b, x, workspace)
    @ccall libacc._SparseSolveSubfactor_Complex_Double(Subfactor::Ptr{SparseOpaqueSubfactor_Complex_Double}, b::Ptr{DenseMatrix_Complex_Double}, x::Ptr{DenseMatrix_Complex_Double}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveOpaque_Complex_Double(Factored, RHS, Soln, workspace)
    @ccall libacc._SparseSolveOpaque_Complex_Double(Factored::Ptr{SparseOpaqueFactorization_Complex_Double}, RHS::Ptr{DenseMatrix_Complex_Double}, Soln::Ptr{DenseMatrix_Complex_Double}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseDestroyOpaqueNumeric_Complex_Double(toFree)
    @ccall libacc._SparseDestroyOpaqueNumeric_Complex_Double(toFree::Ptr{SparseOpaqueFactorization_Complex_Double})::Cvoid
end

function _SparseRetainNumeric_Complex_Double(numericFactor)
    @ccall libacc._SparseRetainNumeric_Complex_Double(numericFactor::Ptr{SparseOpaqueFactorization_Complex_Double})::Cvoid
end

function _SparseGetOptionsFromNumericFactor_Complex_Double(factor)
    @ccall libacc._SparseGetOptionsFromNumericFactor_Complex_Double(factor::Ptr{SparseOpaqueFactorization_Complex_Double})::SparseNumericFactorOptions
end

function _SparseGetWorkspaceRequired_Complex_Double(Subfactor, Factor, workStatic, workPerRHS)
    @ccall libacc._SparseGetWorkspaceRequired_Complex_Double(Subfactor::SparseSubfactor_t, Factor::SparseOpaqueFactorization_Complex_Double, workStatic::Ptr{Csize_t}, workPerRHS::Ptr{Csize_t})::Cvoid
end

function _SparseGetIterativeStateSize_Complex_Double(method, preconditioner, m, n, nrhs)
    @ccall libacc._SparseGetIterativeStateSize_Complex_Double(method::Ptr{SparseIterativeMethod}, preconditioner::Bool, m::Cint, n::Cint, nrhs::Cint)::Csize_t
end

function _SparseCreatePreconditioner_Complex_Double(type, A)
    @ccall libacc._SparseCreatePreconditioner_Complex_Double(type::SparsePreconditioner_t, A::Ptr{SparseMatrix_Complex_Double})::SparseOpaquePreconditioner_Complex_Double
end

function _SparseReleaseOpaquePreconditioner_Complex_Double(toFree)
    @ccall libacc._SparseReleaseOpaquePreconditioner_Complex_Double(toFree::Ptr{SparseOpaquePreconditioner_Complex_Double})::Cvoid
end

function _SparseSpMV_Complex_Double(alpha, A, x, accumulate, y)
    @ccall libacc._SparseSpMV_Complex_Double(alpha::__SPARSE_double_complex, A::SparseMatrix_Complex_Double, x::DenseMatrix_Complex_Double, accumulate::Bool, y::DenseMatrix_Complex_Double)::Cvoid
end

function _SparseConvertFromCoordinate_Complex_Float(m, n, nBlock, blockSize, attributes, row, col, val, storage, workspace)
    @ccall libacc._SparseConvertFromCoordinate_Complex_Float(m::Cint, n::Cint, nBlock::Clong, blockSize::UInt8, attributes::SparseAttributesComplex_t, row::Ptr{Cint}, col::Ptr{Cint}, val::Ptr{__SPARSE_float_complex}, storage::Ptr{Cchar}, workspace::Ptr{Cint})::SparseMatrix_Complex_Float
end

function _SparseNumericFactorSymmetric_Complex_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorSymmetric_Complex_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Float
end

function _SparseNumericFactorHermitian_Complex_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorHermitian_Complex_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Float
end

function _SparseNumericFactorQR_Complex_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorQR_Complex_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Float
end

function _SparseNumericFactorLU_Complex_Float(symbolicFactor, Matrix, options, factorStorage, workspace)
    @ccall libacc._SparseNumericFactorLU_Complex_Float(symbolicFactor::Ptr{SparseOpaqueSymbolicFactorization}, Matrix::Ptr{SparseMatrix_Complex_Float}, options::Ptr{SparseNumericFactorOptions}, factorStorage::Ptr{Cvoid}, workspace::Ptr{Cvoid})::SparseOpaqueFactorization_Complex_Float
end

function _SparseFactorSymmetric_Complex_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorSymmetric_Complex_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Float
end

function _SparseFactorHermitian_Complex_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorHermitian_Complex_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Float
end

function _SparseFactorQR_Complex_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorQR_Complex_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Float
end

function _SparseFactorLU_Complex_Float(factorType, Matrix, sfoptions, nfoptions)
    @ccall libacc._SparseFactorLU_Complex_Float(factorType::SparseFactorization_t, Matrix::Ptr{SparseMatrix_Complex_Float}, sfoptions::Ptr{SparseSymbolicFactorOptions}, nfoptions::Ptr{SparseNumericFactorOptions})::SparseOpaqueFactorization_Complex_Float
end

function _SparseRefactorSymmetric_Complex_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorSymmetric_Complex_Float(Matrix::Ptr{SparseMatrix_Complex_Float}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorHermitian_Complex_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorHermitian_Complex_Float(Matrix::Ptr{SparseMatrix_Complex_Float}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorQR_Complex_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorQR_Complex_Float(Matrix::Ptr{SparseMatrix_Complex_Float}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseRefactorLU_Complex_Float(Matrix, Factorization, nfoptions, workspace)
    @ccall libacc._SparseRefactorLU_Complex_Float(Matrix::Ptr{SparseMatrix_Complex_Float}, Factorization::Ptr{SparseOpaqueFactorization_Complex_Float}, nfoptions::Ptr{SparseNumericFactorOptions}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseUpdatePartialRefactorLU_Complex_Float(Opaque, updateCount, updatedIndices, newMatrix)
    @ccall libacc._SparseUpdatePartialRefactorLU_Complex_Float(Opaque::Ptr{SparseOpaqueFactorization_Complex_Float}, updateCount::Cint, updatedIndices::Ptr{Cint}, newMatrix::SparseMatrix_Complex_Float)::Cvoid
end

function _SparseMultiplySubfactor_Complex_Float(Subfactor, x, y, workspace)
    @ccall libacc._SparseMultiplySubfactor_Complex_Float(Subfactor::Ptr{SparseOpaqueSubfactor_Complex_Float}, x::Ptr{DenseMatrix_Complex_Float}, y::Ptr{DenseMatrix_Complex_Float}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveSubfactor_Complex_Float(Subfactor, b, x, workspace)
    @ccall libacc._SparseSolveSubfactor_Complex_Float(Subfactor::Ptr{SparseOpaqueSubfactor_Complex_Float}, b::Ptr{DenseMatrix_Complex_Float}, x::Ptr{DenseMatrix_Complex_Float}, workspace::Ptr{Cchar})::Cvoid
end

function _SparseSolveOpaque_Complex_Float(Factored, RHS, Soln, workspace)
    @ccall libacc._SparseSolveOpaque_Complex_Float(Factored::Ptr{SparseOpaqueFactorization_Complex_Float}, RHS::Ptr{DenseMatrix_Complex_Float}, Soln::Ptr{DenseMatrix_Complex_Float}, workspace::Ptr{Cvoid})::Cvoid
end

function _SparseDestroyOpaqueNumeric_Complex_Float(toFree)
    @ccall libacc._SparseDestroyOpaqueNumeric_Complex_Float(toFree::Ptr{SparseOpaqueFactorization_Complex_Float})::Cvoid
end

function _SparseRetainNumeric_Complex_Float(numericFactor)
    @ccall libacc._SparseRetainNumeric_Complex_Float(numericFactor::Ptr{SparseOpaqueFactorization_Complex_Float})::Cvoid
end

function _SparseGetOptionsFromNumericFactor_Complex_Float(factor)
    @ccall libacc._SparseGetOptionsFromNumericFactor_Complex_Float(factor::Ptr{SparseOpaqueFactorization_Complex_Float})::SparseNumericFactorOptions
end

function _SparseGetWorkspaceRequired_Complex_Float(Subfactor, Factor, workStatic, workPerRHS)
    @ccall libacc._SparseGetWorkspaceRequired_Complex_Float(Subfactor::SparseSubfactor_t, Factor::SparseOpaqueFactorization_Complex_Float, workStatic::Ptr{Csize_t}, workPerRHS::Ptr{Csize_t})::Cvoid
end

function _SparseGetIterativeStateSize_Complex_Float(method, preconditioner, m, n, nrhs)
    @ccall libacc._SparseGetIterativeStateSize_Complex_Float(method::Ptr{SparseIterativeMethod}, preconditioner::Bool, m::Cint, n::Cint, nrhs::Cint)::Csize_t
end

function _SparseCreatePreconditioner_Complex_Float(type, A)
    @ccall libacc._SparseCreatePreconditioner_Complex_Float(type::SparsePreconditioner_t, A::Ptr{SparseMatrix_Complex_Float})::SparseOpaquePreconditioner_Complex_Float
end

function _SparseReleaseOpaquePreconditioner_Complex_Float(toFree)
    @ccall libacc._SparseReleaseOpaquePreconditioner_Complex_Float(toFree::Ptr{SparseOpaquePreconditioner_Complex_Float})::Cvoid
end

function _SparseSpMV_Complex_Float(alpha, A, x, accumulate, y)
    @ccall libacc._SparseSpMV_Complex_Float(alpha::__SPARSE_float_complex, A::SparseMatrix_Complex_Float, x::DenseMatrix_Complex_Float, accumulate::Bool, y::DenseMatrix_Complex_Float)::Cvoid
end

@enum var"##Ctag#299"::UInt32 begin
    BNNSDataTypeFloatBit = 65536
    BNNSDataTypeFloat16 = 65552
    BNNSDataTypeFloat32 = 65568
    BNNSDataTypeBFloat16 = 98320
    BNNSDataTypeIntBit = 131072
    BNNSDataTypeInt1 = 131073
    BNNSDataTypeInt2 = 131074
    BNNSDataTypeInt4 = 131076
    BNNSDataTypeInt8 = 131080
    BNNSDataTypeInt16 = 131088
    BNNSDataTypeInt32 = 131104
    BNNSDataTypeInt64 = 131136
    BNNSDataTypeUIntBit = 262144
    BNNSDataTypeUInt1 = 262145
    BNNSDataTypeUInt2 = 262146
    BNNSDataTypeUInt3 = 262147
    BNNSDataTypeUInt4 = 262148
    BNNSDataTypeUInt6 = 262150
    BNNSDataTypeUInt8 = 262152
    BNNSDataTypeUInt16 = 262160
    BNNSDataTypeUInt32 = 262176
    BNNSDataTypeUInt64 = 262208
    BNNSDataTypeIndexedBit = 524288
    BNNSDataTypeIndexed1 = 524289
    BNNSDataTypeIndexed2 = 524290
    BNNSDataTypeIndexed4 = 524292
    BNNSDataTypeIndexed8 = 524296
    BNNSDataTypeMiscellaneousBit = 1048576
    BNNSDataTypeBoolean = 1048584
end

const BNNSDataType = UInt32

@enum var"##Ctag#300"::UInt32 begin
    BNNSPoolingFunctionMax = 0
    BNNSPoolingFunctionAverageCountIncludePadding = 1
    BNNSPoolingFunctionAverageCountExcludePadding = 2
    BNNSPoolingFunctionUnMax = 3
    BNNSPoolingFunctionL2Norm = 4
    # BNNSPoolingFunctionAverage = 1
end

const BNNSPoolingFunction = UInt32

@enum var"##Ctag#301"::UInt32 begin
    BNNSActivationFunctionIdentity = 0
    BNNSActivationFunctionRectifiedLinear = 1
    BNNSActivationFunctionLeakyRectifiedLinear = 2
    BNNSActivationFunctionSigmoid = 3
    BNNSActivationFunctionTanh = 4
    BNNSActivationFunctionScaledTanh = 5
    BNNSActivationFunctionAbs = 6
    BNNSActivationFunctionLinear = 7
    BNNSActivationFunctionClamp = 8
    BNNSActivationFunctionIntegerLinearSaturate = 9
    BNNSActivationFunctionIntegerLinearSaturatePerChannel = 10
    BNNSActivationFunctionSoftmax = 11
    BNNSActivationFunctionGELUApproximation = 12
    BNNSActivationFunctionGumbel = 13
    BNNSActivationFunctionGumbelMax = 14
    BNNSActivationFunctionHardSigmoid = 15
    BNNSActivationFunctionSoftplus = 16
    BNNSActivationFunctionSoftsign = 17
    BNNSActivationFunctionELU = 18
    BNNSActivationFunctionClampedLeakyRectifiedLinear = 19
    BNNSActivationFunctionLinearWithBias = 20
    BNNSActivationFunctionLogSoftmax = 21
    BNNSActivationFunctionLogSigmoid = 22
    BNNSActivationFunctionSELU = 23
    BNNSActivationFunctionCELU = 24
    BNNSActivationFunctionHardShrink = 25
    BNNSActivationFunctionSoftShrink = 26
    BNNSActivationFunctionTanhShrink = 27
    BNNSActivationFunctionThreshold = 28
    BNNSActivationFunctionPReLUPerChannel = 29
    BNNSActivationFunctionGELUApproximation2 = 30
    # BNNSActivationFunctionHardSwish = 30
    BNNSActivationFunctionSiLU = 31
    BNNSActivationFunctionReLU6 = 32
    BNNSActivationFunctionErf = 33
    BNNSActivationFunctionGELU = 34
    BNNSActivationFunctionGELUApproximationSigmoid = 35
end

const BNNSActivationFunction = UInt32

@enum var"##Ctag#302"::UInt32 begin
    BNNSFlagsUseClientPtr = 1
end

const BNNSFlags = UInt32

@enum var"##Ctag#303"::UInt32 begin
    BNNSLossFunctionSoftmaxCrossEntropy = 1
    BNNSLossFunctionSigmoidCrossEntropy = 2
    BNNSLossFunctionMeanSquareError = 3
    BNNSLossFunctionHuber = 4
    BNNSLossFunctionYolo = 5
    BNNSLossFunctionLog = 6
    BNNSLossFunctionCosineDistance = 7
    BNNSLossFunctionHinge = 8
    BNNSLossFunctionMeanAbsoluteError = 9
    BNNSLossFunctionCategoricalCrossEntropy = 10
end

const BNNSLossFunction = UInt32

@enum var"##Ctag#304"::UInt32 begin
    BNNSLossReductionNone = 0
    BNNSLossReductionSum = 1
    BNNSLossReductionWeightedMean = 2
    BNNSLossReductionMean = 3
    BNNSLossReductionNonZeroWeightMean = 4
end

const BNNSLossReductionFunction = UInt32

@enum var"##Ctag#305"::UInt32 begin
    BNNSArithmeticAdd = 0
    BNNSArithmeticSubtract = 1
    BNNSArithmeticMultiply = 2
    BNNSArithmeticDivide = 3
    BNNSArithmeticSquareRoot = 4
    BNNSArithmeticReciprocalSquareRoot = 5
    BNNSArithmeticCeil = 6
    BNNSArithmeticFloor = 7
    BNNSArithmeticRound = 8
    BNNSArithmeticSin = 9
    BNNSArithmeticCos = 10
    BNNSArithmeticTan = 11
    BNNSArithmeticAsin = 12
    BNNSArithmeticAcos = 13
    BNNSArithmeticAtan = 14
    BNNSArithmeticSinh = 15
    BNNSArithmeticCosh = 16
    BNNSArithmeticTanh = 17
    BNNSArithmeticAsinh = 18
    BNNSArithmeticAcosh = 19
    BNNSArithmeticAtanh = 20
    BNNSArithmeticPow = 21
    BNNSArithmeticExp = 22
    BNNSArithmeticExp2 = 23
    BNNSArithmeticLog = 24
    BNNSArithmeticLog2 = 25
    BNNSArithmeticMultiplyNoNaN = 26
    BNNSArithmeticDivideNoNaN = 27
    BNNSArithmeticMultiplyAdd = 28
    BNNSArithmeticMinimum = 29
    BNNSArithmeticMaximum = 30
    BNNSArithmeticSelect = 31
    BNNSArithmeticAbs = 32
    BNNSArithmeticSign = 33
    BNNSArithmeticNegate = 34
    BNNSArithmeticReciprocal = 35
    BNNSArithmeticSquare = 36
    BNNSArithmeticFloorDivide = 37
    BNNSArithmeticTruncDivide = 38
    BNNSArithmeticTruncRemainder = 39
    BNNSArithmeticErf = 40
end

const BNNSArithmeticFunction = UInt32

@enum var"##Ctag#306"::UInt32 begin
    BNNSConstant = 0
    BNNSSample = 1
    BNNSParameter = 2
end

const BNNSDescriptorType = UInt32

@enum var"##Ctag#307"::UInt32 begin
    BNNSOptimizerFunctionSGDMomentum = 1
    BNNSOptimizerFunctionAdam = 2
    BNNSOptimizerFunctionRMSProp = 3
    BNNSOptimizerFunctionAdamW = 4
    BNNSOptimizerFunctionAdamAMSGrad = 5
    BNNSOptimizerFunctionAdamWAMSGrad = 6
    BNNSOptimizerFunctionSGDMomentumWithClipping = 7
    BNNSOptimizerFunctionAdamWithClipping = 8
    BNNSOptimizerFunctionRMSPropWithClipping = 9
    BNNSOptimizerFunctionAdamWWithClipping = 10
    BNNSOptimizerFunctionAdamAMSGradWithClipping = 11
    BNNSOptimizerFunctionAdamWAMSGradWithClipping = 12
end

const BNNSOptimizerFunction = UInt32

@enum var"##Ctag#308"::UInt32 begin
    BNNSOptimizerRegularizationNone = 0
    BNNSOptimizerRegularizationL1 = 1
    BNNSOptimizerRegularizationL2 = 2
end

const BNNSOptimizerRegularizationFunction = UInt32

@enum var"##Ctag#309"::UInt32 begin
    BNNSSGDMomentumVariant0 = 0
    BNNSSGDMomentumVariant1 = 1
    BNNSSGDMomentumVariant2 = 2
end

const BNNSOptimizerSGDMomentumVariant = UInt32

@enum var"##Ctag#310"::UInt32 begin
    BNNSOptimizerClippingNone = 0
    BNNSOptimizerClippingByValue = 1
    BNNSOptimizerClippingByNorm = 2
    BNNSOptimizerClippingByGlobalNorm = 3
end

const BNNSOptimizerClippingFunction = UInt32

@enum var"##Ctag#311"::UInt32 begin
    BNNSL2Norm = 1
end

const BNNSNormType = UInt32

@enum var"##Ctag#312"::UInt32 begin
    BNNSConvolution = 0
    BNNSFullyConnected = 1
    BNNSBatchNorm = 2
    BNNSInstanceNorm = 3
    BNNSLayerNorm = 4
    BNNSGroupNorm = 5
    BNNSTransposedConvolution = 6
    BNNSQuantization = 7
    BNNSArithmetic = 8
end

const BNNSFilterType = UInt32

@enum var"##Ctag#313"::UInt32 begin
    BNNSReduceFunctionMax = 0
    BNNSReduceFunctionMin = 1
    BNNSReduceFunctionArgMax = 2
    BNNSReduceFunctionArgMin = 3
    BNNSReduceFunctionMean = 4
    BNNSReduceFunctionMeanNonZero = 5
    BNNSReduceFunctionSum = 6
    BNNSReduceFunctionSumSquare = 7
    BNNSReduceFunctionSumLog = 8
    BNNSReduceFunctionL1Norm = 9
    BNNSReduceFunctionLogicalOr = 10
    BNNSReduceFunctionLogicalAnd = 11
    BNNSReduceFunctionL2Norm = 12
    BNNSReduceFunctionLogSumExp = 13
    BNNSReduceFunctionProduct = 14
    BNNSReduceFunctionNone = 15
    BNNSReduceFunctionLogSum = 16
    # BNNSReduceFunctionAny = 10
    # BNNSReduceFunctionAll = 11
end

const BNNSReduceFunction = UInt32

@enum var"##Ctag#314"::UInt32 begin
    BNNSLayerFlagsLSTMBidirectional = 1
    BNNSLayerFlagsLSTMDefaultActivations = 2
end

const BNNSLayerFlags = UInt32

@enum var"##Ctag#315"::UInt32 begin
    BNNSDataLayoutVector = 65536
    BNNSDataLayout1DLastMajor = 98304
    BNNSDataLayout1DFirstMajor = 98305
    BNNSDataLayoutRowMajorMatrix = 131072
    BNNSDataLayoutColumnMajorMatrix = 131073
    BNNSDataLayout2DLastMajor = 163840
    BNNSDataLayout2DFirstMajor = 163841
    BNNSDataLayoutFullyConnectedSparse = 135169
    BNNSDataLayoutImageCHW = 196608
    BNNSDataLayoutSNE = 196609
    BNNSDataLayoutNSE = 196610
    BNNSDataLayoutMHA_DHK = 196611
    BNNSDataLayout3DLastMajor = 229376
    BNNSDataLayout3DFirstMajor = 229377
    BNNSDataLayoutConvolutionWeightsOIHW = 262144
    BNNSDataLayoutConvolutionWeightsOIHrWr = 262145
    BNNSDataLayoutConvolutionWeightsIOHrWr = 262146
    BNNSDataLayoutConvolutionWeightsOIHW_Pack32 = 262160
    BNNSDataLayout4DLastMajor = 294912
    BNNSDataLayout4DFirstMajor = 294913
    BNNSDataLayout5DLastMajor = 360448
    BNNSDataLayout5DFirstMajor = 360449
    BNNSDataLayout6DLastMajor = 425984
    BNNSDataLayout6DFirstMajor = 425985
    BNNSDataLayout7DLastMajor = 491520
    BNNSDataLayout7DFirstMajor = 491521
    BNNSDataLayout8DLastMajor = 557056
    BNNSDataLayout8DFirstMajor = 557057
end

const BNNSDataLayout = UInt32

@enum var"##Ctag#316"::UInt32 begin
    BNNSInterpolationMethodNearest = 0
    BNNSInterpolationMethodLinear = 1
end

const BNNSInterpolationMethod = UInt32

@enum var"##Ctag#317"::UInt32 begin
    BNNSLinearSamplingDefault = 0
    BNNSLinearSamplingAlignCorners = 1
    BNNSLinearSamplingUnalignCorners = 2
    BNNSLinearSamplingStrictAlignCorners = 3
    BNNSLinearSamplingOffsetCorners = 4
end

const BNNSLinearSamplingMode = UInt32

@enum var"##Ctag#318"::UInt32 begin
    BNNSCornersHeightFirst = 0
    BNNSCornersWidthFirst = 1
    BNNSCenterSizeHeightFirst = 2
    BNNSCenterSizeWidthFirst = 3
end

const BNNSBoxCoordinateMode = UInt32

@enum var"##Ctag#319"::UInt32 begin
    BNNSPaddingModeConstant = 0
    BNNSPaddingModeReflect = 1
    BNNSPaddingModeSymmetric = 2
end

const BNNSPaddingMode = UInt32

@enum var"##Ctag#320"::UInt32 begin
    BNNSRelationalOperatorEqual = 0
    BNNSRelationalOperatorLess = 1
    BNNSRelationalOperatorLessEqual = 2
    BNNSRelationalOperatorGreater = 3
    BNNSRelationalOperatorGreaterEqual = 4
    BNNSRelationalOperatorNotEqual = 5
    BNNSRelationalOperatorLogicalAND = 6
    BNNSRelationalOperatorLogicalOR = 7
    BNNSRelationalOperatorLogicalNOT = 8
    BNNSRelationalOperatorLogicalNAND = 9
    BNNSRelationalOperatorLogicalNOR = 10
    BNNSRelationalOperatorLogicalXOR = 11
end

const BNNSRelationalOperator = UInt32

@enum var"##Ctag#321"::UInt32 begin
    BNNSPointerSpecifierAlpha = 0
    BNNSPointerSpecifierBeta = 1
end

const BNNSPointerSpecifier = UInt32

@enum var"##Ctag#322"::UInt32 begin
    BNNSNDArrayFlagBackpropSet = 0
    BNNSNDArrayFlagBackpropAccumulate = 1
end

const BNNSNDArrayFlags = UInt32

@enum var"##Ctag#323"::UInt32 begin
    BNNSEmbeddingFlagScaleGradientByFrequency = 1
end

const BNNSEmbeddingFlags = UInt32

@enum var"##Ctag#324"::UInt32 begin
    BNNSQuantizerFunctionQuantize = 0
    BNNSQuantizerFunctionDequantize = 1
end

const BNNSQuantizerFunction = UInt32

@enum var"##Ctag#325"::UInt32 begin
    BNNSRandomGeneratorMethodAES_CTR = 0
end

const BNNSRandomGeneratorMethod = UInt32

@enum var"##Ctag#326"::UInt32 begin
    BNNSSparsityTypeUnstructured = 0
end

const BNNSSparsityType = UInt32

@enum var"##Ctag#327"::UInt32 begin
    BNNSTargetSystemGeneric = 0
end

const BNNSTargetSystem = UInt32

@enum var"##Ctag#328"::UInt32 begin
    BNNSShuffleTypePixelShuffleNCHW = 0
    BNNSShuffleTypePixelUnshuffleNCHW = 1
    BNNSShuffleTypeDepthToSpaceNCHW = 2
    BNNSShuffleTypeSpaceToDepthNCHW = 3
end

const BNNSShuffleType = UInt32

# typedef int ( * BNNSAlloc ) ( void * _Nullable * _Nullable memptr , size_t alignment , size_t size )
const BNNSAlloc = Ptr{Cvoid}

# typedef void ( * BNNSFree ) ( void * _Null_unspecified ptr )
const BNNSFree = Ptr{Cvoid}

struct BNNSActivation
    _function::BNNSActivationFunction
    alpha::Cfloat
    beta::Cfloat
    iscale::Int32
    ioffset::Int32
    ishift::Int32
    iscale_per_channel::Ptr{Int32}
    ioffset_per_channel::Ptr{Int32}
    ishift_per_channel::Ptr{Int32}
end

struct BNNSNDArrayDescriptor
    flags::BNNSNDArrayFlags
    layout::BNNSDataLayout
    size::NTuple{8, Csize_t}
    stride::NTuple{8, Csize_t}
    data::Ptr{Cvoid}
    data_type::BNNSDataType
    table_data::Ptr{Cvoid}
    table_data_type::BNNSDataType
    data_scale::Cfloat
    data_bias::Cfloat
end

struct BNNSTensor
    data_type::BNNSDataType
    rank::UInt8
    shape::NTuple{8, Cssize_t}
    stride::NTuple{8, Cssize_t}
    data::Ptr{Cvoid}
    data_size_in_bytes::Csize_t
    name::Ptr{Cchar}
end

struct BNNSLSTMGateDescriptor
    iw_desc::NTuple{2, BNNSNDArrayDescriptor}
    hw_desc::BNNSNDArrayDescriptor
    cw_desc::BNNSNDArrayDescriptor
    b_desc::BNNSNDArrayDescriptor
    activation::BNNSActivation
end

struct BNNSLSTMDataDescriptor
    data_desc::BNNSNDArrayDescriptor
    hidden_desc::BNNSNDArrayDescriptor
    cell_state_desc::BNNSNDArrayDescriptor
end

struct BNNSArithmeticUnary
    in::BNNSNDArrayDescriptor
    in_type::BNNSDescriptorType
    out::BNNSNDArrayDescriptor
    out_type::BNNSDescriptorType
end

struct BNNSArithmeticBinary
    in1::BNNSNDArrayDescriptor
    in1_type::BNNSDescriptorType
    in2::BNNSNDArrayDescriptor
    in2_type::BNNSDescriptorType
    out::BNNSNDArrayDescriptor
    out_type::BNNSDescriptorType
end

struct BNNSArithmeticTernary
    in1::BNNSNDArrayDescriptor
    in1_type::BNNSDescriptorType
    in2::BNNSNDArrayDescriptor
    in2_type::BNNSDescriptorType
    in3::BNNSNDArrayDescriptor
    in3_type::BNNSDescriptorType
    out::BNNSNDArrayDescriptor
    out_type::BNNSDescriptorType
end

struct BNNSMHAProjectionParameters
    target_desc::BNNSNDArrayDescriptor
    weights::BNNSNDArrayDescriptor
    bias::BNNSNDArrayDescriptor
end

struct BNNSLayerParametersConvolution
    i_desc::BNNSNDArrayDescriptor
    w_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    bias::BNNSNDArrayDescriptor
    activation::BNNSActivation
    x_stride::Csize_t
    y_stride::Csize_t
    x_dilation_stride::Csize_t
    y_dilation_stride::Csize_t
    x_padding::Csize_t
    y_padding::Csize_t
    groups::Csize_t
    pad::NTuple{4, Csize_t}
end

struct BNNSLayerParametersFullyConnected
    i_desc::BNNSNDArrayDescriptor
    w_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    bias::BNNSNDArrayDescriptor
    activation::BNNSActivation
end

struct BNNSLayerParametersPooling
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    bias::BNNSNDArrayDescriptor
    activation::BNNSActivation
    pooling_function::BNNSPoolingFunction
    k_width::Csize_t
    k_height::Csize_t
    x_stride::Csize_t
    y_stride::Csize_t
    x_dilation_stride::Csize_t
    y_dilation_stride::Csize_t
    x_padding::Csize_t
    y_padding::Csize_t
    pad::NTuple{4, Csize_t}
end

struct BNNSLayerParametersActivation
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    activation::BNNSActivation
    axis_flags::UInt32
end

struct BNNSLayerParametersLossBase
    _function::BNNSLossFunction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    reduction::BNNSLossReductionFunction
end

struct BNNSLayerParametersLossSoftmaxCrossEntropy
    _function::BNNSLossFunction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    reduction::BNNSLossReductionFunction
    label_smooth::Cfloat
end

struct BNNSLayerParametersLossSigmoidCrossEntropy
    _function::BNNSLossFunction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    reduction::BNNSLossReductionFunction
    label_smooth::Cfloat
end

struct BNNSLayerParametersLossHuber
    _function::BNNSLossFunction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    reduction::BNNSLossReductionFunction
    huber_delta::Cfloat
end

struct BNNSLayerParametersLossYolo
    _function::BNNSLossFunction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    reduction::BNNSLossReductionFunction
    huber_delta::Cfloat
    number_of_grid_columns::Csize_t
    number_of_grid_rows::Csize_t
    number_of_anchor_boxes::Csize_t
    anchor_box_size::Csize_t
    rescore::Bool
    scale_xy::Cfloat
    scale_wh::Cfloat
    scale_object::Cfloat
    scale_no_object::Cfloat
    scale_classification::Cfloat
    object_minimum_iou::Cfloat
    no_object_maximum_iou::Cfloat
    anchors_data::Ptr{Cfloat}
end

struct BNNSOptimizerSGDMomentumFields
    learning_rate::Cfloat
    momentum::Cfloat
    gradient_scale::Cfloat
    regularization_scale::Cfloat
    clip_gradients::Bool
    clip_gradients_min::Cfloat
    clip_gradients_max::Cfloat
    nesterov::Bool
    regularization_func::BNNSOptimizerRegularizationFunction
    sgd_momentum_variant::BNNSOptimizerSGDMomentumVariant
end

struct BNNSOptimizerSGDMomentumWithClippingFields
    learning_rate::Cfloat
    momentum::Cfloat
    gradient_scale::Cfloat
    regularization_scale::Cfloat
    nesterov::Bool
    regularization_func::BNNSOptimizerRegularizationFunction
    sgd_momentum_variant::BNNSOptimizerSGDMomentumVariant
    clipping_func::BNNSOptimizerClippingFunction
    clip_gradients_min::Cfloat
    clip_gradients_max::Cfloat
    clip_gradients_max_norm::Cfloat
    clip_gradients_use_norm::Cfloat
end

struct BNNSOptimizerAdamFields
    learning_rate::Cfloat
    beta1::Cfloat
    beta2::Cfloat
    time_step::Cfloat
    epsilon::Cfloat
    gradient_scale::Cfloat
    regularization_scale::Cfloat
    clip_gradients::Bool
    clip_gradients_min::Cfloat
    clip_gradients_max::Cfloat
    regularization_func::BNNSOptimizerRegularizationFunction
end

struct BNNSOptimizerAdamWithClippingFields
    learning_rate::Cfloat
    beta1::Cfloat
    beta2::Cfloat
    time_step::Cfloat
    epsilon::Cfloat
    gradient_scale::Cfloat
    regularization_scale::Cfloat
    regularization_func::BNNSOptimizerRegularizationFunction
    clipping_func::BNNSOptimizerClippingFunction
    clip_gradients_min::Cfloat
    clip_gradients_max::Cfloat
    clip_gradients_max_norm::Cfloat
    clip_gradients_use_norm::Cfloat
end

struct BNNSOptimizerRMSPropFields
    learning_rate::Cfloat
    alpha::Cfloat
    epsilon::Cfloat
    centered::Bool
    momentum::Cfloat
    gradient_scale::Cfloat
    regularization_scale::Cfloat
    clip_gradients::Bool
    clip_gradients_min::Cfloat
    clip_gradients_max::Cfloat
    regularization_func::BNNSOptimizerRegularizationFunction
end

struct BNNSOptimizerRMSPropWithClippingFields
    learning_rate::Cfloat
    alpha::Cfloat
    epsilon::Cfloat
    centered::Bool
    momentum::Cfloat
    gradient_scale::Cfloat
    regularization_scale::Cfloat
    regularization_func::BNNSOptimizerRegularizationFunction
    clipping_func::BNNSOptimizerClippingFunction
    clip_gradients_min::Cfloat
    clip_gradients_max::Cfloat
    clip_gradients_max_norm::Cfloat
    clip_gradients_use_norm::Cfloat
end

struct BNNSLayerParametersNormalization
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    beta_desc::BNNSNDArrayDescriptor
    gamma_desc::BNNSNDArrayDescriptor
    moving_mean_desc::BNNSNDArrayDescriptor
    moving_variance_desc::BNNSNDArrayDescriptor
    momentum::Cfloat
    epsilon::Cfloat
    activation::BNNSActivation
    num_groups::Csize_t
    normalization_axis::Csize_t
end

struct BNNSLayerParametersDropout
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    rate::Cfloat
    seed::UInt32
    control::UInt8
end

struct BNNSLayerParametersLSTM
    data::NTuple{5040, UInt8}
end

function Base.getproperty(x::Ptr{BNNSLayerParametersLSTM}, f::Symbol)
    f === :input_size && return Ptr{Csize_t}(x + 0)
    f === :hidden_size && return Ptr{Csize_t}(x + 8)
    f === :batch_size && return Ptr{Csize_t}(x + 16)
    f === :num_layers && return Ptr{Csize_t}(x + 24)
    f === :seq_len && return Ptr{Csize_t}(x + 32)
    f === :dropout && return Ptr{Cfloat}(x + 40)
    f === :lstm_flags && return Ptr{UInt32}(x + 44)
    f === :sequence_descriptor && return Ptr{BNNSNDArrayDescriptor}(x + 48)
    f === :input_descriptor && return Ptr{BNNSLSTMDataDescriptor}(x + 224)
    f === :output_descriptor && return Ptr{BNNSLSTMDataDescriptor}(x + 752)
    f === :input_gate && return Ptr{BNNSLSTMGateDescriptor}(x + 1280)
    f === :forget_gate && return Ptr{BNNSLSTMGateDescriptor}(x + 2208)
    f === :candidate_gate && return Ptr{BNNSLSTMGateDescriptor}(x + 3136)
    f === :output_gate && return Ptr{BNNSLSTMGateDescriptor}(x + 4064)
    f === :hidden_activation && return Ptr{BNNSActivation}(x + 4992)
    return getfield(x, f)
end

function Base.getproperty(x::BNNSLayerParametersLSTM, f::Symbol)
    r = Ref{BNNSLayerParametersLSTM}(x)
    ptr = Base.unsafe_convert(Ptr{BNNSLayerParametersLSTM}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{BNNSLayerParametersLSTM}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::BNNSLayerParametersLSTM, private::Bool = false)
    (:input_size, :hidden_size, :batch_size, :num_layers, :seq_len, :dropout, :lstm_flags, :sequence_descriptor, :input_descriptor, :output_descriptor, :input_gate, :forget_gate, :candidate_gate, :output_gate, :hidden_activation, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct BNNSLayerParametersArithmetic
    arithmetic_function::BNNSArithmeticFunction
    arithmetic_function_fields::Ptr{Cvoid}
    activation::BNNSActivation
end

struct BNNSLayerParametersPermute
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    permutation::NTuple{8, Csize_t}
end

struct BNNSLayerParametersTensorContraction
    operation::Ptr{Cchar}
    alpha::Cfloat
    beta::Cfloat
    iA_desc::BNNSNDArrayDescriptor
    iB_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
end

struct BNNSLayerParametersGram
    alpha::Cfloat
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
end

struct BNNSLayerParametersResize
    method::BNNSInterpolationMethod
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    align_corners::Bool
end

struct BNNSLayerParametersCropResize
    normalized_coordinates::Bool
    spatial_scale::Cfloat
    extrapolation_value::Cfloat
    sampling_mode::BNNSLinearSamplingMode
    box_coordinate_mode::BNNSBoxCoordinateMode
    method::BNNSInterpolationMethod
end

struct BNNSLayerParametersBroadcastMatMul
    alpha::Cfloat
    beta::Cfloat
    transA::Bool
    transB::Bool
    quadratic::Bool
    a_is_weights::Bool
    b_is_weights::Bool
    iA_desc::BNNSNDArrayDescriptor
    iB_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
end

struct BNNSLayerParametersMultiheadAttention
    query::BNNSMHAProjectionParameters
    key::BNNSMHAProjectionParameters
    value::BNNSMHAProjectionParameters
    add_zero_attn::Bool
    key_attn_bias::BNNSNDArrayDescriptor
    value_attn_bias::BNNSNDArrayDescriptor
    output::BNNSMHAProjectionParameters
    dropout::Cfloat
    seed::UInt32
end

struct BNNSLayerParametersReduction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    w_desc::BNNSNDArrayDescriptor
    reduce_func::BNNSReduceFunction
    epsilon::Cfloat
end

struct BNNSLayerParametersPadding
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    padding_size::NTuple{8, NTuple{2, Csize_t}}
    padding_mode::BNNSPaddingMode
    padding_value::UInt32
end

struct BNNSLayerParametersEmbedding
    flags::BNNSEmbeddingFlags
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    dictionary::BNNSNDArrayDescriptor
    padding_idx::Csize_t
    max_norm::Cfloat
    norm_type::Cfloat
end

struct BNNSLayerParametersQuantization
    axis_mask::Csize_t
    _function::BNNSQuantizerFunction
    i_desc::BNNSNDArrayDescriptor
    o_desc::BNNSNDArrayDescriptor
    scale::BNNSNDArrayDescriptor
    bias::BNNSNDArrayDescriptor
end

struct BNNSSparsityParameters
    flags::UInt64
    sparsity_ratio::NTuple{2, UInt32}
    sparsity_type::BNNSSparsityType
    target_system::BNNSTargetSystem
end

struct BNNSImageStackDescriptor
    width::Csize_t
    height::Csize_t
    channels::Csize_t
    row_stride::Csize_t
    image_stride::Csize_t
    data_type::BNNSDataType
    data_scale::Cfloat
    data_bias::Cfloat
end

struct BNNSVectorDescriptor
    size::Csize_t
    data_type::BNNSDataType
    data_scale::Cfloat
    data_bias::Cfloat
end

struct BNNSLayerData
    data::Ptr{Cvoid}
    data_type::BNNSDataType
    data_scale::Cfloat
    data_bias::Cfloat
    data_table::Ptr{Cfloat}
end

struct BNNSConvolutionLayerParameters
    data::NTuple{176, UInt8}
end

function Base.getproperty(x::Ptr{BNNSConvolutionLayerParameters}, f::Symbol)
    f === :x_stride && return Ptr{Csize_t}(x + 0)
    f === :y_stride && return Ptr{Csize_t}(x + 8)
    f === :x_padding && return Ptr{Csize_t}(x + 16)
    f === :y_padding && return Ptr{Csize_t}(x + 24)
    f === :k_width && return Ptr{Csize_t}(x + 32)
    f === :k_height && return Ptr{Csize_t}(x + 40)
    f === :in_channels && return Ptr{Csize_t}(x + 48)
    f === :out_channels && return Ptr{Csize_t}(x + 56)
    f === :weights && return Ptr{BNNSLayerData}(x + 64)
    f === :bias && return Ptr{BNNSLayerData}(x + 96)
    f === :activation && return Ptr{BNNSActivation}(x + 128)
    return getfield(x, f)
end

function Base.getproperty(x::BNNSConvolutionLayerParameters, f::Symbol)
    r = Ref{BNNSConvolutionLayerParameters}(x)
    ptr = Base.unsafe_convert(Ptr{BNNSConvolutionLayerParameters}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{BNNSConvolutionLayerParameters}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::BNNSConvolutionLayerParameters, private::Bool = false)
    (:x_stride, :y_stride, :x_padding, :y_padding, :k_width, :k_height, :in_channels, :out_channels, :weights, :bias, :activation, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct BNNSFullyConnectedLayerParameters
    data::NTuple{128, UInt8}
end

function Base.getproperty(x::Ptr{BNNSFullyConnectedLayerParameters}, f::Symbol)
    f === :in_size && return Ptr{Csize_t}(x + 0)
    f === :out_size && return Ptr{Csize_t}(x + 8)
    f === :weights && return Ptr{BNNSLayerData}(x + 16)
    f === :bias && return Ptr{BNNSLayerData}(x + 48)
    f === :activation && return Ptr{BNNSActivation}(x + 80)
    return getfield(x, f)
end

function Base.getproperty(x::BNNSFullyConnectedLayerParameters, f::Symbol)
    r = Ref{BNNSFullyConnectedLayerParameters}(x)
    ptr = Base.unsafe_convert(Ptr{BNNSFullyConnectedLayerParameters}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{BNNSFullyConnectedLayerParameters}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::BNNSFullyConnectedLayerParameters, private::Bool = false)
    (:in_size, :out_size, :weights, :bias, :activation, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct BNNSPoolingLayerParameters
    data::NTuple{152, UInt8}
end

function Base.getproperty(x::Ptr{BNNSPoolingLayerParameters}, f::Symbol)
    f === :x_stride && return Ptr{Csize_t}(x + 0)
    f === :y_stride && return Ptr{Csize_t}(x + 8)
    f === :x_padding && return Ptr{Csize_t}(x + 16)
    f === :y_padding && return Ptr{Csize_t}(x + 24)
    f === :k_width && return Ptr{Csize_t}(x + 32)
    f === :k_height && return Ptr{Csize_t}(x + 40)
    f === :in_channels && return Ptr{Csize_t}(x + 48)
    f === :out_channels && return Ptr{Csize_t}(x + 56)
    f === :pooling_function && return Ptr{BNNSPoolingFunction}(x + 64)
    f === :bias && return Ptr{BNNSLayerData}(x + 72)
    f === :activation && return Ptr{BNNSActivation}(x + 104)
    return getfield(x, f)
end

function Base.getproperty(x::BNNSPoolingLayerParameters, f::Symbol)
    r = Ref{BNNSPoolingLayerParameters}(x)
    ptr = Base.unsafe_convert(Ptr{BNNSPoolingLayerParameters}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{BNNSPoolingLayerParameters}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

function Base.propertynames(x::BNNSPoolingLayerParameters, private::Bool = false)
    (:x_stride, :y_stride, :x_padding, :y_padding, :k_width, :k_height, :in_channels, :out_channels, :pooling_function, :bias, :activation, if private
            fieldnames(typeof(x))
        else
            ()
        end...)
end

struct BNNSFilterParameters
    flags::UInt32
    n_threads::Csize_t
    alloc_memory::BNNSAlloc
    free_memory::BNNSFree
end

# typedef int ( * bnns_graph_realloc_fn_t ) ( void * _Nullable user_memory_context , size_t user_memory_context_size , void * _Nullable * _Nonnull memptr , size_t alignment , size_t size )
const bnns_graph_realloc_fn_t = Ptr{Cvoid}

# typedef void ( * bnns_graph_free_all_fn_t ) ( void * _Nullable user_memory_context , size_t user_memory_context_size )
const bnns_graph_free_all_fn_t = Ptr{Cvoid}

@enum var"##Ctag#329"::UInt32 begin
    BNNSGraphMessageLevelInfo = 1
    BNNSGraphMessageLevelUnsupported = 2
    BNNSGraphMessageLevelWarning = 4
    BNNSGraphMessageLevelError = 8
end

const BNNSGraphMessageLevel = UInt32

# typedef void ( * bnns_graph_execute_message_fn_t ) ( BNNSGraphMessageLevel msg_level , char const * _Nonnull error_msg , char const * _Nullable op_info , bnns_user_message_data_t * _Nullable additional_logging_arguments )
const bnns_graph_execute_message_fn_t = Ptr{Cvoid}

# typedef void ( * bnns_graph_compile_message_fn_t ) ( BNNSGraphMessageLevel msg_level , char const * _Nonnull error_msg , char const * _Nullable source_location , bnns_user_message_data_t * _Nullable additional_logging_arguments )
const bnns_graph_compile_message_fn_t = Ptr{Cvoid}

function BNNSGraphCompileOptionsMakeDefault()
    @ccall libacc.BNNSGraphCompileOptionsMakeDefault()::bnns_graph_compile_options_t
end

function BNNSGraphCompileOptionsDestroy(options)
    @ccall libacc.BNNSGraphCompileOptionsDestroy(options::bnns_graph_compile_options_t)::Cvoid
end

function BNNSGraphCompileOptionsSetTargetSingleThread(options, value)
    @ccall libacc.BNNSGraphCompileOptionsSetTargetSingleThread(options::bnns_graph_compile_options_t, value::Bool)::Cvoid
end

function BNNSGraphCompileOptionsGetTargetSingleThread(options)
    @ccall libacc.BNNSGraphCompileOptionsGetTargetSingleThread(options::bnns_graph_compile_options_t)::Bool
end

function BNNSGraphCompileOptionsSetGenerateDebugInfo(options, value)
    @ccall libacc.BNNSGraphCompileOptionsSetGenerateDebugInfo(options::bnns_graph_compile_options_t, value::Bool)::Cvoid
end

function BNNSGraphCompileOptionsGetGenerateDebugInfo(options)
    @ccall libacc.BNNSGraphCompileOptionsGetGenerateDebugInfo(options::bnns_graph_compile_options_t)::Bool
end

@enum var"##Ctag#330"::UInt32 begin
    BNNSGraphOptimizationPreferencePerformance = 0
    BNNSGraphOptimizationPreferenceIRSize = 1
end

const BNNSGraphOptimizationPreference = UInt32

function BNNSGraphCompileOptionsSetOptimizationPreference(options, preference)
    @ccall libacc.BNNSGraphCompileOptionsSetOptimizationPreference(options::bnns_graph_compile_options_t, preference::BNNSGraphOptimizationPreference)::Cvoid
end

function BNNSGraphCompileOptionsGetOptimizationPreference(options)
    @ccall libacc.BNNSGraphCompileOptionsGetOptimizationPreference(options::bnns_graph_compile_options_t)::BNNSGraphOptimizationPreference
end

function BNNSGraphCompileOptionsSetMessageLogCallback(options, log_callback, additional_logging_arguments)
    @ccall libacc.BNNSGraphCompileOptionsSetMessageLogCallback(options::bnns_graph_compile_options_t, log_callback::bnns_graph_compile_message_fn_t, additional_logging_arguments::Ptr{bnns_user_message_data_t})::Cvoid
end

function BNNSGraphCompileOptionsSetMessageLogMask(options, log_level_mask)
    @ccall libacc.BNNSGraphCompileOptionsSetMessageLogMask(options::bnns_graph_compile_options_t, log_level_mask::UInt32)::Cvoid
end

function BNNSGraphCompileOptionsSetOutputPath(options, path)
    @ccall libacc.BNNSGraphCompileOptionsSetOutputPath(options::bnns_graph_compile_options_t, path::Ptr{Cchar})::Cvoid
end

function BNNSGraphCompileOptionsGetOutputPath(options)
    @ccall libacc.BNNSGraphCompileOptionsGetOutputPath(options::bnns_graph_compile_options_t)::Ptr{Cchar}
end

function BNNSGraphCompileOptionsSetOutputFD(options, fd)
    @ccall libacc.BNNSGraphCompileOptionsSetOutputFD(options::bnns_graph_compile_options_t, fd::Cint)::Cvoid
end

function BNNSGraphCompileOptionsGetOutputFD(options)
    @ccall libacc.BNNSGraphCompileOptionsGetOutputFD(options::bnns_graph_compile_options_t)::Cint
end

function BNNSGraphCompileFromFile(filename, _function, options)
    @ccall libacc.BNNSGraphCompileFromFile(filename::Ptr{Cchar}, _function::Ptr{Cchar}, options::bnns_graph_compile_options_t)::bnns_graph_t
end

function BNNSGraphGetInputCount(graph, _function)
    @ccall libacc.BNNSGraphGetInputCount(graph::bnns_graph_t, _function::Ptr{Cchar})::Csize_t
end

function BNNSGraphGetOutputCount(graph, _function)
    @ccall libacc.BNNSGraphGetOutputCount(graph::bnns_graph_t, _function::Ptr{Cchar})::Csize_t
end

function BNNSGraphGetArgumentCount(graph, _function)
    @ccall libacc.BNNSGraphGetArgumentCount(graph::bnns_graph_t, _function::Ptr{Cchar})::Csize_t
end

function BNNSGraphGetFunctionCount(graph)
    @ccall libacc.BNNSGraphGetFunctionCount(graph::bnns_graph_t)::Csize_t
end

function BNNSGraphGetInputNames(graph, _function, input_names_count, input_names)
    @ccall libacc.BNNSGraphGetInputNames(graph::bnns_graph_t, _function::Ptr{Cchar}, input_names_count::Csize_t, input_names::Ptr{Ptr{Cchar}})::Cint
end

function BNNSGraphGetOutputNames(graph, _function, output_names_count, output_names)
    @ccall libacc.BNNSGraphGetOutputNames(graph::bnns_graph_t, _function::Ptr{Cchar}, output_names_count::Csize_t, output_names::Ptr{Ptr{Cchar}})::Cint
end

function BNNSGraphGetArgumentNames(graph, _function, argument_names_count, argument_names)
    @ccall libacc.BNNSGraphGetArgumentNames(graph::bnns_graph_t, _function::Ptr{Cchar}, argument_names_count::Csize_t, argument_names::Ptr{Ptr{Cchar}})::Cint
end

function BNNSGraphGetFunctionNames(graph, function_name_count, function_names)
    @ccall libacc.BNNSGraphGetFunctionNames(graph::bnns_graph_t, function_name_count::Csize_t, function_names::Ptr{Ptr{Cchar}})::Cint
end

@enum var"##Ctag#331"::UInt32 begin
    BNNSGraphArgumentIntentIn = 1
    BNNSGraphArgumentIntentOut = 2
    BNNSGraphArgumentIntentInOut = 3
end

const BNNSGraphArgumentIntent = UInt32

function BNNSGraphGetArgumentIntents(graph, _function, argument_intents_count, argument_intents)
    @ccall libacc.BNNSGraphGetArgumentIntents(graph::bnns_graph_t, _function::Ptr{Cchar}, argument_intents_count::Csize_t, argument_intents::Ptr{BNNSGraphArgumentIntent})::Cint
end

function BNNSGraphGetArgumentPosition(graph, _function, argument)
    @ccall libacc.BNNSGraphGetArgumentPosition(graph::bnns_graph_t, _function::Ptr{Cchar}, argument::Ptr{Cchar})::Csize_t
end

function BNNSGraphGetArgumentInterleaveFactors(graph, _function, argument_count, argument_interleave, argument_interleave_counts)
    @ccall libacc.BNNSGraphGetArgumentInterleaveFactors(graph::bnns_graph_t, _function::Ptr{Cchar}, argument_count::Csize_t, argument_interleave::Ptr{Ptr{UInt16}}, argument_interleave_counts::Ptr{Csize_t})::Cint
end

@enum var"##Ctag#332"::UInt32 begin
    BNNSGraphArgumentTypePointer = 0
    BNNSGraphArgumentTypeTensor = 2
end

const BNNSGraphArgumentType = UInt32

function BNNSGraphContextMake(graph)
    @ccall libacc.BNNSGraphContextMake(graph::bnns_graph_t)::bnns_graph_context_t
end

function BNNSGraphContextMakeStreaming(graph, _function, initial_states_count, initial_states)
    @ccall libacc.BNNSGraphContextMakeStreaming(graph::bnns_graph_t, _function::Ptr{Cchar}, initial_states_count::Csize_t, initial_states::Ptr{BNNSTensor})::bnns_graph_context_t
end

function BNNSGraphContextDestroy(context)
    @ccall libacc.BNNSGraphContextDestroy(context::bnns_graph_context_t)::Cvoid
end

function BNNSGraphContextSetDynamicShapes(context, _function, shapes_count, shapes)
    @ccall libacc.BNNSGraphContextSetDynamicShapes(context::bnns_graph_context_t, _function::Ptr{Cchar}, shapes_count::Csize_t, shapes::Ptr{bnns_graph_shape_t})::Cint
end

function BNNSGraphContextSetArgumentType(context, argument_type)
    @ccall libacc.BNNSGraphContextSetArgumentType(context::bnns_graph_context_t, argument_type::BNNSGraphArgumentType)::Cint
end

function BNNSGraphContextEnableNanAndInfChecks(context, enable_check_for_nans_inf)
    @ccall libacc.BNNSGraphContextEnableNanAndInfChecks(context::bnns_graph_context_t, enable_check_for_nans_inf::Bool)::Cvoid
end

function BNNSGraphContextSetStreamingAdvanceCount(context, advance_count)
    @ccall libacc.BNNSGraphContextSetStreamingAdvanceCount(context::bnns_graph_context_t, advance_count::Csize_t)::Cint
end

function BNNSGraphContextExecute(context, _function, argument_count, arguments, workspace_size, workspace)
    @ccall libacc.BNNSGraphContextExecute(context::bnns_graph_context_t, _function::Ptr{Cchar}, argument_count::Csize_t, arguments::Ptr{bnns_graph_argument_t}, workspace_size::Csize_t, workspace::Ptr{Cchar})::Cint
end

function BNNSGraphContextGetWorkspaceSize(context, _function)
    @ccall libacc.BNNSGraphContextGetWorkspaceSize(context::bnns_graph_context_t, _function::Ptr{Cchar})::Csize_t
end

function BNNSGraphContextGetTensor(context, _function, argument, fill_known_dynamic_shapes, tensor)
    @ccall libacc.BNNSGraphContextGetTensor(context::bnns_graph_context_t, _function::Ptr{Cchar}, argument::Ptr{Cchar}, fill_known_dynamic_shapes::Bool, tensor::Ptr{BNNSTensor})::Cint
end

function BNNSGraphTensorFillStrides(graph, _function, argument, tensor)
    @ccall libacc.BNNSGraphTensorFillStrides(graph::bnns_graph_t, _function::Ptr{Cchar}, argument::Ptr{Cchar}, tensor::Ptr{BNNSTensor})::Cint
end

function BNNSGraphContextSetMessageLogCallback(context, log_callback_fn, additional_logging_arguments)
    @ccall libacc.BNNSGraphContextSetMessageLogCallback(context::bnns_graph_context_t, log_callback_fn::bnns_graph_execute_message_fn_t, additional_logging_arguments::Ptr{bnns_user_message_data_t})::Cint
end

function BNNSGraphContextSetMessageLogMask(context, log_level_mask)
    @ccall libacc.BNNSGraphContextSetMessageLogMask(context::bnns_graph_context_t, log_level_mask::UInt32)::Cint
end

const BNNSFilter = Ptr{Cvoid}

function BNNSFilterCreateLayerConvolution(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerConvolution(layer_params::Ptr{BNNSLayerParametersConvolution}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerTransposedConvolution(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerTransposedConvolution(layer_params::Ptr{BNNSLayerParametersConvolution}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerFullyConnected(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerFullyConnected(layer_params::Ptr{BNNSLayerParametersFullyConnected}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerPooling(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerPooling(layer_params::Ptr{BNNSLayerParametersPooling}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerActivation(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerActivation(layer_params::Ptr{BNNSLayerParametersActivation}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerLoss(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerLoss(layer_params::Ptr{Cvoid}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerNormalization(normType, layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerNormalization(normType::BNNSFilterType, layer_params::Ptr{BNNSLayerParametersNormalization}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerArithmetic(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerArithmetic(layer_params::Ptr{BNNSLayerParametersArithmetic}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerPermute(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerPermute(layer_params::Ptr{BNNSLayerParametersPermute}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerDropout(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerDropout(layer_params::Ptr{BNNSLayerParametersDropout}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerPadding(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerPadding(layer_params::Ptr{BNNSLayerParametersPadding}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerBroadcastMatMul(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerBroadcastMatMul(layer_params::Ptr{BNNSLayerParametersBroadcastMatMul}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerTensorContraction(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerTensorContraction(layer_params::Ptr{BNNSLayerParametersTensorContraction}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerGram(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerGram(layer_params::Ptr{BNNSLayerParametersGram}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerResize(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerResize(layer_params::Ptr{BNNSLayerParametersResize}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerMultiheadAttention(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerMultiheadAttention(layer_params::Ptr{BNNSLayerParametersMultiheadAttention}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerReduction(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerReduction(layer_params::Ptr{BNNSLayerParametersReduction}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateFusedLayer(number_of_fused_filters, filter_type, layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateFusedLayer(number_of_fused_filters::Csize_t, filter_type::Ptr{BNNSFilterType}, layer_params::Ptr{Ptr{Cvoid}}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateLayerEmbedding(layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateLayerEmbedding(layer_params::Ptr{BNNSLayerParametersEmbedding}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterApply(filter, in, out)
    @ccall libacc.BNNSFilterApply(filter::Ptr{Cvoid}, in::Ptr{Cvoid}, out::Ptr{Cvoid})::Cint
end

function BNNSFilterApplyBatch(filter, batch_size, in, in_stride, out, out_stride)
    @ccall libacc.BNNSFilterApplyBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t)::Cint
end

function BNNSPoolingFilterApplyBatch(filter, batch_size, in, in_stride, out, out_stride, indices, idx_stride)
    @ccall libacc.BNNSPoolingFilterApplyBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, indices::Ptr{Csize_t}, idx_stride::Csize_t)::Cint
end

function BNNSPoolingFilterApplyBatchEx(filter, batch_size, in, in_stride, out, out_stride, indices_data_type, indices, idx_stride)
    @ccall libacc.BNNSPoolingFilterApplyBatchEx(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, indices_data_type::BNNSDataType, indices::Ptr{Cvoid}, idx_stride::Csize_t)::Cint
end

function BNNSFilterApplyTwoInput(filter, inA, inB, out)
    @ccall libacc.BNNSFilterApplyTwoInput(filter::Ptr{Cvoid}, inA::Ptr{Cvoid}, inB::Ptr{Cvoid}, out::Ptr{Cvoid})::Cint
end

function BNNSFilterApplyTwoInputBatch(filter, batch_size, inA, inA_stride, inB, inB_stride, out, out_stride)
    @ccall libacc.BNNSFilterApplyTwoInputBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, inA::Ptr{Cvoid}, inA_stride::Csize_t, inB::Ptr{Cvoid}, inB_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t)::Cint
end

function BNNSNormalizationFilterApplyBatch(filter, batch_size, in, in_stride, out, out_stride, training)
    @ccall libacc.BNNSNormalizationFilterApplyBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, training::Bool)::Cint
end

function BNNSFusedFilterApplyBatch(filter, batch_size, in, in_stride, out, out_stride, training)
    @ccall libacc.BNNSFusedFilterApplyBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, training::Bool)::Cint
end

function BNNSFusedFilterApplyMultiInputBatch(filter, batch_size, number_of_inputs, in, in_stride, out, out_stride, training)
    @ccall libacc.BNNSFusedFilterApplyMultiInputBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, number_of_inputs::Csize_t, in::Ptr{Ptr{Cvoid}}, in_stride::Ptr{Csize_t}, out::Ptr{Cvoid}, out_stride::Csize_t, training::Bool)::Cint
end

function BNNSArithmeticFilterApplyBatch(filter, batch_size, number_of_inputs, in, in_stride, out, out_stride)
    @ccall libacc.BNNSArithmeticFilterApplyBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, number_of_inputs::Csize_t, in::Ptr{Ptr{Cvoid}}, in_stride::Ptr{Csize_t}, out::Ptr{Cvoid}, out_stride::Csize_t)::Cint
end

function BNNSApplyMultiheadAttention(F, batch_size, query, query_stride, key, key_stride, key_mask, key_mask_stride, value, value_stride, output, output_stride, add_to_attention, backprop_cache_size, backprop_cache, workspace_size, workspace)
    @ccall libacc.BNNSApplyMultiheadAttention(F::Ptr{Cvoid}, batch_size::Csize_t, query::Ptr{Cvoid}, query_stride::Csize_t, key::Ptr{Cvoid}, key_stride::Csize_t, key_mask::Ptr{BNNSNDArrayDescriptor}, key_mask_stride::Csize_t, value::Ptr{Cvoid}, value_stride::Csize_t, output::Ptr{Cvoid}, output_stride::Csize_t, add_to_attention::Ptr{BNNSNDArrayDescriptor}, backprop_cache_size::Ptr{Csize_t}, backprop_cache::Ptr{Cvoid}, workspace_size::Ptr{Csize_t}, workspace::Ptr{Cvoid})::Cint
end

function BNNSDirectApplyQuantizer(layer_params, filter_params, batch_size, input_stride, output_stride)
    @ccall libacc.BNNSDirectApplyQuantizer(layer_params::Ptr{BNNSLayerParametersQuantization}, filter_params::Ptr{BNNSFilterParameters}, batch_size::Csize_t, input_stride::Csize_t, output_stride::Csize_t)::Cint
end

function BNNSFilterDestroy(filter)
    @ccall libacc.BNNSFilterDestroy(filter::Ptr{Cvoid})::Cvoid
end

function BNNSOptimizerStep(_function, OptimizerAlgFields, number_of_parameters, parameters, gradients, accumulators, filter_params)
    @ccall libacc.BNNSOptimizerStep(_function::BNNSOptimizerFunction, OptimizerAlgFields::Ptr{Cvoid}, number_of_parameters::Csize_t, parameters::Ptr{Ptr{BNNSNDArrayDescriptor}}, gradients::Ptr{Ptr{BNNSNDArrayDescriptor}}, accumulators::Ptr{Ptr{BNNSNDArrayDescriptor}}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSClipByValue(dest, src, min_val, max_val)
    @ccall libacc.BNNSClipByValue(dest::Ptr{BNNSNDArrayDescriptor}, src::Ptr{BNNSNDArrayDescriptor}, min_val::Cfloat, max_val::Cfloat)::Cint
end

function BNNSClipByNorm(dest, src, max_norm, axis_flags)
    @ccall libacc.BNNSClipByNorm(dest::Ptr{BNNSNDArrayDescriptor}, src::Ptr{BNNSNDArrayDescriptor}, max_norm::Cfloat, axis_flags::UInt32)::Cint
end

function BNNSClipByGlobalNorm(dest, src, count, max_norm, use_norm)
    @ccall libacc.BNNSClipByGlobalNorm(dest::Ptr{Ptr{BNNSNDArrayDescriptor}}, src::Ptr{Ptr{BNNSNDArrayDescriptor}}, count::Csize_t, max_norm::Cfloat, use_norm::Cfloat)::Cint
end

function BNNSComputeNorm(dest, src, norm_type, axis_flags)
    @ccall libacc.BNNSComputeNorm(dest::Ptr{BNNSNDArrayDescriptor}, src::Ptr{BNNSNDArrayDescriptor}, norm_type::BNNSNormType, axis_flags::UInt32)::Cint
end

function BNNSComputeNormBackward(in, in_delta, out, out_delta, norm_type, axis_flags)
    @ccall libacc.BNNSComputeNormBackward(in::Ptr{Cvoid}, in_delta::Ptr{BNNSNDArrayDescriptor}, out::Ptr{Cvoid}, out_delta::Ptr{BNNSNDArrayDescriptor}, norm_type::BNNSNormType, axis_flags::UInt32)::Cint
end

function BNNSFilterApplyBackwardBatch(filter, batch_size, in, in_stride, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride, weights_delta, bias_delta)
    @ccall libacc.BNNSFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, weights_delta::Ptr{BNNSNDArrayDescriptor}, bias_delta::Ptr{BNNSNDArrayDescriptor})::Cint
end

function BNNSPoolingFilterApplyBackwardBatch(filter, batch_size, in, in_stride, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride, bias_delta, indices, idx_stride)
    @ccall libacc.BNNSPoolingFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, bias_delta::Ptr{BNNSNDArrayDescriptor}, indices::Ptr{Csize_t}, idx_stride::Csize_t)::Cint
end

function BNNSPoolingFilterApplyBackwardBatchEx(filter, batch_size, in, in_stride, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride, bias_delta, indices_data_type, indices, idx_stride)
    @ccall libacc.BNNSPoolingFilterApplyBackwardBatchEx(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, bias_delta::Ptr{BNNSNDArrayDescriptor}, indices_data_type::BNNSDataType, indices::Ptr{Cvoid}, idx_stride::Csize_t)::Cint
end

function BNNSFilterApplyBackwardTwoInputBatch(filter, batch_size, inA, inA_stride, inA_delta, inA_delta_stride, inB, inB_stride, inB_delta, inB_delta_stride, out, out_stride, out_delta, out_delta_stride, weights_delta, bias_delta)
    @ccall libacc.BNNSFilterApplyBackwardTwoInputBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, inA::Ptr{Cvoid}, inA_stride::Csize_t, inA_delta::Ptr{BNNSNDArrayDescriptor}, inA_delta_stride::Csize_t, inB::Ptr{Cvoid}, inB_stride::Csize_t, inB_delta::Ptr{BNNSNDArrayDescriptor}, inB_delta_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, weights_delta::Ptr{BNNSNDArrayDescriptor}, bias_delta::Ptr{BNNSNDArrayDescriptor})::Cint
end

function BNNSNormalizationFilterApplyBackwardBatch(filter, batch_size, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride, beta_delta, gamma_delta)
    @ccall libacc.BNNSNormalizationFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, beta_delta::Ptr{BNNSNDArrayDescriptor}, gamma_delta::Ptr{BNNSNDArrayDescriptor})::Cint
end

function BNNSFusedFilterApplyBackwardBatch(filter, batch_size, in, in_stride, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride, delta_parameters)
    @ccall libacc.BNNSFusedFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, delta_parameters::Ptr{Ptr{BNNSNDArrayDescriptor}})::Cint
end

function BNNSFusedFilterApplyBackwardMultiInputBatch(filter, batch_size, number_of_inputs, in, in_stride, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride, delta_parameters)
    @ccall libacc.BNNSFusedFilterApplyBackwardMultiInputBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, number_of_inputs::Csize_t, in::Ptr{Ptr{Cvoid}}, in_stride::Ptr{Csize_t}, in_delta::Ptr{Ptr{BNNSNDArrayDescriptor}}, in_delta_stride::Ptr{Csize_t}, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t, delta_parameters::Ptr{Ptr{BNNSNDArrayDescriptor}})::Cint
end

function BNNSArithmeticFilterApplyBackwardBatch(filter, batch_size, number_of_inputs, in, in_stride, in_delta, in_delta_stride, out, out_stride, out_delta, out_delta_stride)
    @ccall libacc.BNNSArithmeticFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, number_of_inputs::Csize_t, in::Ptr{Ptr{Cvoid}}, in_stride::Ptr{Csize_t}, in_delta::Ptr{Ptr{BNNSNDArrayDescriptor}}, in_delta_stride::Ptr{Csize_t}, out::Ptr{Cvoid}, out_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t)::Cint
end

function BNNSPermuteFilterApplyBackwardBatch(filter, batch_size, in_delta, in_delta_stride, out_delta, out_delta_stride)
    @ccall libacc.BNNSPermuteFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t)::Cint
end

function BNNSLossFilterApplyBatch(filter, batch_size, in, in_stride, labels, labels_stride, weights, weights_size, out, in_delta, in_delta_stride)
    @ccall libacc.BNNSLossFilterApplyBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, labels::Ptr{Cvoid}, labels_stride::Csize_t, weights::Ptr{Cvoid}, weights_size::Csize_t, out::Ptr{Cvoid}, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t)::Cint
end

function BNNSLossFilterApplyBackwardBatch(filter, batch_size, in, in_stride, in_delta, in_delta_stride, labels, labels_stride, weights, weights_size, out_delta, out_delta_stride)
    @ccall libacc.BNNSLossFilterApplyBackwardBatch(filter::Ptr{Cvoid}, batch_size::Csize_t, in::Ptr{Cvoid}, in_stride::Csize_t, in_delta::Ptr{BNNSNDArrayDescriptor}, in_delta_stride::Csize_t, labels::Ptr{Cvoid}, labels_stride::Csize_t, weights::Ptr{Cvoid}, weights_size::Csize_t, out_delta::Ptr{BNNSNDArrayDescriptor}, out_delta_stride::Csize_t)::Cint
end

function BNNSApplyMultiheadAttentionBackward(F, batch_size, query, query_stride, query_param_delta, key, key_stride, key_mask, key_mask_stride, key_param_delta, value, value_stride, value_param_delta, add_to_attention, key_attn_bias_delta, value_attn_bias_delta, output, output_stride, output_param_delta, backprop_cache_size, backprop_cache, workspace_size, workspace)
    @ccall libacc.BNNSApplyMultiheadAttentionBackward(F::Ptr{Cvoid}, batch_size::Csize_t, query::Ptr{Cvoid}, query_stride::Csize_t, query_param_delta::Ptr{BNNSMHAProjectionParameters}, key::Ptr{Cvoid}, key_stride::Csize_t, key_mask::Ptr{BNNSNDArrayDescriptor}, key_mask_stride::Csize_t, key_param_delta::Ptr{BNNSMHAProjectionParameters}, value::Ptr{Cvoid}, value_stride::Csize_t, value_param_delta::Ptr{BNNSMHAProjectionParameters}, add_to_attention::Ptr{BNNSNDArrayDescriptor}, key_attn_bias_delta::Ptr{BNNSNDArrayDescriptor}, value_attn_bias_delta::Ptr{BNNSNDArrayDescriptor}, output::Ptr{Cvoid}, output_stride::Csize_t, output_param_delta::Ptr{BNNSMHAProjectionParameters}, backprop_cache_size::Csize_t, backprop_cache::Ptr{Cvoid}, workspace_size::Ptr{Csize_t}, workspace::Ptr{Cvoid})::Cint
end

function BNNSComputeLSTMTrainingCacheCapacity(layer_params)
    @ccall libacc.BNNSComputeLSTMTrainingCacheCapacity(layer_params::Ptr{BNNSLayerParametersLSTM})::Csize_t
end

function BNNSDirectApplyLSTMBatchTrainingCaching(layer_params, filter_params, training_cache_ptr, training_cache_capacity)
    @ccall libacc.BNNSDirectApplyLSTMBatchTrainingCaching(layer_params::Ptr{BNNSLayerParametersLSTM}, filter_params::Ptr{BNNSFilterParameters}, training_cache_ptr::Ptr{Cvoid}, training_cache_capacity::Csize_t)::Cint
end

function BNNSDirectApplyActivationBatch(layer_params, filter_params, batch_size, in_stride, out_stride)
    @ccall libacc.BNNSDirectApplyActivationBatch(layer_params::Ptr{BNNSLayerParametersActivation}, filter_params::Ptr{BNNSFilterParameters}, batch_size::Csize_t, in_stride::Csize_t, out_stride::Csize_t)::Cint
end

function BNNSCopy(dest, src, filter_params)
    @ccall libacc.BNNSCopy(dest::Ptr{BNNSNDArrayDescriptor}, src::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSMatMulWorkspaceSize(transA, transB, alpha, inputA, inputB, output, filter_params)
    @ccall libacc.BNNSMatMulWorkspaceSize(transA::Bool, transB::Bool, alpha::Cfloat, inputA::Ptr{BNNSNDArrayDescriptor}, inputB::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cssize_t
end

function BNNSMatMul(transA, transB, alpha, inputA, inputB, output, workspace, filter_params)
    @ccall libacc.BNNSMatMul(transA::Bool, transB::Bool, alpha::Cfloat, inputA::Ptr{BNNSNDArrayDescriptor}, inputB::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, workspace::Ptr{Cvoid}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSDirectApplyBroadcastMatMul(transA, transB, alpha, inputA, inputB, output, filter_params)
    @ccall libacc.BNNSDirectApplyBroadcastMatMul(transA::Bool, transB::Bool, alpha::Cfloat, inputA::Ptr{BNNSNDArrayDescriptor}, inputB::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cvoid
end

function BNNSTranspose(dest, src, axis0, axis1, filter_params)
    @ccall libacc.BNNSTranspose(dest::Ptr{BNNSNDArrayDescriptor}, src::Ptr{BNNSNDArrayDescriptor}, axis0::Csize_t, axis1::Csize_t, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSDirectApplyReduction(layer_params, filter_params)
    @ccall libacc.BNNSDirectApplyReduction(layer_params::Ptr{BNNSLayerParametersReduction}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSCompareTensor(in0, in1, op, out)
    @ccall libacc.BNNSCompareTensor(in0::Ptr{BNNSNDArrayDescriptor}, in1::Ptr{BNNSNDArrayDescriptor}, op::BNNSRelationalOperator, out::Ptr{BNNSNDArrayDescriptor})::Cint
end

function BNNSTile(input, output, filter_params)
    @ccall libacc.BNNSTile(input::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSDirectApplyTopK(K, axis, batch_size, input, input_batch_stride, best_values, best_values_batch_stride, best_indices, best_indices_batch_stride, filter_params)
    @ccall libacc.BNNSDirectApplyTopK(K::Csize_t, axis::Csize_t, batch_size::Csize_t, input::Ptr{BNNSNDArrayDescriptor}, input_batch_stride::Csize_t, best_values::Ptr{BNNSNDArrayDescriptor}, best_values_batch_stride::Csize_t, best_indices::Ptr{BNNSNDArrayDescriptor}, best_indices_batch_stride::Csize_t, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSDirectApplyInTopK(K, axis, batch_size, input, input_batch_stride, test_indices, test_indices_batch_stride, output, output_batch_stride, filter_params)
    @ccall libacc.BNNSDirectApplyInTopK(K::Csize_t, axis::Csize_t, batch_size::Csize_t, input::Ptr{BNNSNDArrayDescriptor}, input_batch_stride::Csize_t, test_indices::Ptr{BNNSNDArrayDescriptor}, test_indices_batch_stride::Csize_t, output::Ptr{BNNSNDArrayDescriptor}, output_batch_stride::Csize_t, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSGather(axis, input, indices, output, filter_params)
    @ccall libacc.BNNSGather(axis::Csize_t, input::Ptr{BNNSNDArrayDescriptor}, indices::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSScatter(axis, op, input, indices, output, filter_params)
    @ccall libacc.BNNSScatter(axis::Csize_t, op::BNNSReduceFunction, input::Ptr{BNNSNDArrayDescriptor}, indices::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSGatherND(input, indices, output, filter_params)
    @ccall libacc.BNNSGatherND(input::Ptr{BNNSNDArrayDescriptor}, indices::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSScatterND(op, input, indices, output, filter_params)
    @ccall libacc.BNNSScatterND(op::BNNSReduceFunction, input::Ptr{BNNSNDArrayDescriptor}, indices::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSShuffle(type, input, output, filter_params)
    @ccall libacc.BNNSShuffle(type::BNNSShuffleType, input::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSBandPart(num_lower, num_upper, input, output, filter_params)
    @ccall libacc.BNNSBandPart(num_lower::Cint, num_upper::Cint, input::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSCropResize(layer_params, input, roi, output, filter_params)
    @ccall libacc.BNNSCropResize(layer_params::Ptr{BNNSLayerParametersCropResize}, input::Ptr{BNNSNDArrayDescriptor}, roi::Ptr{BNNSNDArrayDescriptor}, output::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSDirectApplyLSTMBatchBackward(layer_params, layer_delta_params, filter_params, training_cache_ptr, training_cache_capacity)
    @ccall libacc.BNNSDirectApplyLSTMBatchBackward(layer_params::Ptr{BNNSLayerParametersLSTM}, layer_delta_params::Ptr{BNNSLayerParametersLSTM}, filter_params::Ptr{BNNSFilterParameters}, training_cache_ptr::Ptr{Cvoid}, training_cache_capacity::Csize_t)::Cint
end

function BNNSTileBackward(in_delta, out_delta, filter_params)
    @ccall libacc.BNNSTileBackward(in_delta::Ptr{BNNSNDArrayDescriptor}, out_delta::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSCropResizeBackward(layer_params, in_delta, roi, out_delta, filter_params)
    @ccall libacc.BNNSCropResizeBackward(layer_params::Ptr{BNNSLayerParametersCropResize}, in_delta::Ptr{BNNSNDArrayDescriptor}, roi::Ptr{BNNSNDArrayDescriptor}, out_delta::Ptr{BNNSNDArrayDescriptor}, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSGetPointer(filter, target)
    @ccall libacc.BNNSGetPointer(filter::Ptr{Cvoid}, target::BNNSPointerSpecifier)::BNNSNDArrayDescriptor
end

function BNNSNDArrayGetDataSize(array)
    @ccall libacc.BNNSNDArrayGetDataSize(array::Ptr{BNNSNDArrayDescriptor})::Csize_t
end

function BNNSTensorGetAllocationSize(tensor)
    @ccall libacc.BNNSTensorGetAllocationSize(tensor::Ptr{BNNSTensor})::Csize_t
end

function BNNSDataLayoutGetRank(layout)
    @ccall libacc.BNNSDataLayoutGetRank(layout::BNNSDataLayout)::Csize_t
end

function BNNSNDArrayFullyConnectedSparsifySparseCOO(in_dense_shape, in_indices, in_values, out, sparse_params, batch_size, workspace, workspace_size, filter_params)
    @ccall libacc.BNNSNDArrayFullyConnectedSparsifySparseCOO(in_dense_shape::Ptr{BNNSNDArrayDescriptor}, in_indices::Ptr{BNNSNDArrayDescriptor}, in_values::Ptr{BNNSNDArrayDescriptor}, out::Ptr{BNNSNDArrayDescriptor}, sparse_params::Ptr{BNNSSparsityParameters}, batch_size::Csize_t, workspace::Ptr{Cvoid}, workspace_size::Csize_t, filter_params::Ptr{BNNSFilterParameters})::Cint
end

function BNNSNDArrayFullyConnectedSparsifySparseCSR(in_dense_shape, in_column_indices, in_row_starts, in_values, out, sparse_params, batch_size, workspace, workspace_size, filter_params)
    @ccall libacc.BNNSNDArrayFullyConnectedSparsifySparseCSR(in_dense_shape::Ptr{BNNSNDArrayDescriptor}, in_column_indices::Ptr{BNNSNDArrayDescriptor}, in_row_starts::Ptr{BNNSNDArrayDescriptor}, in_values::Ptr{BNNSNDArrayDescriptor}, out::Ptr{BNNSNDArrayDescriptor}, sparse_params::Ptr{BNNSSparsityParameters}, batch_size::Csize_t, workspace::Ptr{Cvoid}, workspace_size::Csize_t, filter_params::Ptr{BNNSFilterParameters})::Cint
end

const BNNSRandomGenerator = Ptr{Cvoid}

function BNNSCreateRandomGenerator(method, filter_params)
    @ccall libacc.BNNSCreateRandomGenerator(method::BNNSRandomGeneratorMethod, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSCreateRandomGeneratorWithSeed(method, seed, filter_params)
    @ccall libacc.BNNSCreateRandomGeneratorWithSeed(method::BNNSRandomGeneratorMethod, seed::UInt64, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSDestroyRandomGenerator(generator)
    @ccall libacc.BNNSDestroyRandomGenerator(generator::Ptr{Cvoid})::Cvoid
end

function BNNSRandomGeneratorStateSize(generator)
    @ccall libacc.BNNSRandomGeneratorStateSize(generator::Ptr{Cvoid})::Csize_t
end

function BNNSRandomGeneratorGetState(generator, state_size, state)
    @ccall libacc.BNNSRandomGeneratorGetState(generator::Ptr{Cvoid}, state_size::Csize_t, state::Ptr{Cvoid})::Cint
end

function BNNSRandomGeneratorSetState(generator, state_size, state)
    @ccall libacc.BNNSRandomGeneratorSetState(generator::Ptr{Cvoid}, state_size::Csize_t, state::Ptr{Cvoid})::Cint
end

function BNNSRandomFillUniformFloat(generator, desc, a, b)
    @ccall libacc.BNNSRandomFillUniformFloat(generator::Ptr{Cvoid}, desc::Ptr{BNNSNDArrayDescriptor}, a::Cfloat, b::Cfloat)::Cint
end

function BNNSRandomFillUniformInt(generator, desc, a, b)
    @ccall libacc.BNNSRandomFillUniformInt(generator::Ptr{Cvoid}, desc::Ptr{BNNSNDArrayDescriptor}, a::Int64, b::Int64)::Cint
end

function BNNSRandomFillNormalFloat(generator, desc, mean, stddev)
    @ccall libacc.BNNSRandomFillNormalFloat(generator::Ptr{Cvoid}, desc::Ptr{BNNSNDArrayDescriptor}, mean::Cfloat, stddev::Cfloat)::Cint
end

function BNNSRandomFillCategoricalFloat(generator, desc, probabilities, log_probabilities)
    @ccall libacc.BNNSRandomFillCategoricalFloat(generator::Ptr{Cvoid}, desc::Ptr{BNNSNDArrayDescriptor}, probabilities::Ptr{BNNSNDArrayDescriptor}, log_probabilities::Bool)::Cint
end

const BNNSNearestNeighbors = Ptr{Cvoid}

function BNNSCreateNearestNeighbors(max_n_samples, n_features, n_neighbors, data_type, filter_params)
    @ccall libacc.BNNSCreateNearestNeighbors(max_n_samples::Cuint, n_features::Cuint, n_neighbors::Cuint, data_type::BNNSDataType, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSDestroyNearestNeighbors(knn)
    @ccall libacc.BNNSDestroyNearestNeighbors(knn::Ptr{Cvoid})::Cvoid
end

function BNNSNearestNeighborsLoad(knn, n_new_samples, data_ptr)
    @ccall libacc.BNNSNearestNeighborsLoad(knn::Ptr{Cvoid}, n_new_samples::Cuint, data_ptr::Ptr{Cvoid})::Cint
end

function BNNSNearestNeighborsGetInfo(knn, sample_number, indices, distances)
    @ccall libacc.BNNSNearestNeighborsGetInfo(knn::Ptr{Cvoid}, sample_number::Cint, indices::Ptr{Cint}, distances::Ptr{Cvoid})::Cint
end

function BNNSFilterCreateConvolutionLayer(in_desc, out_desc, layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateConvolutionLayer(in_desc::Ptr{BNNSImageStackDescriptor}, out_desc::Ptr{BNNSImageStackDescriptor}, layer_params::Ptr{BNNSConvolutionLayerParameters}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateFullyConnectedLayer(in_desc, out_desc, layer_params, filter_params)
    @ccall libacc.BNNSFilterCreateFullyConnectedLayer(in_desc::Ptr{BNNSVectorDescriptor}, out_desc::Ptr{BNNSVectorDescriptor}, layer_params::Ptr{BNNSFullyConnectedLayerParameters}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreatePoolingLayer(in_desc, out_desc, layer_params, filter_params)
    @ccall libacc.BNNSFilterCreatePoolingLayer(in_desc::Ptr{BNNSImageStackDescriptor}, out_desc::Ptr{BNNSImageStackDescriptor}, layer_params::Ptr{BNNSPoolingLayerParameters}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

function BNNSFilterCreateVectorActivationLayer(in_desc, out_desc, activation, filter_params)
    @ccall libacc.BNNSFilterCreateVectorActivationLayer(in_desc::Ptr{BNNSVectorDescriptor}, out_desc::Ptr{BNNSVectorDescriptor}, activation::Ptr{BNNSActivation}, filter_params::Ptr{BNNSFilterParameters})::Ptr{Cvoid}
end

struct var"##Ctag#333"
    data::Ptr{Cvoid}
    size::Cint
end
function Base.getproperty(x::Ptr{var"##Ctag#333"}, f::Symbol)
    f === :data && return Ptr{Ptr{Cvoid}}(x + 0)
    f === :size && return Ptr{Cint}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#333", f::Symbol)
    r = Ref{var"##Ctag#333"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#333"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#333"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#334"
    data::Ptr{Cvoid}
    size::Cint
end
function Base.getproperty(x::Ptr{var"##Ctag#334"}, f::Symbol)
    f === :data && return Ptr{Ptr{Cvoid}}(x + 0)
    f === :size && return Ptr{Cint}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#334", f::Symbol)
    r = Ref{var"##Ctag#334"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#334"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#334"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#335"
    data::Ptr{Cvoid}
    size::Cint
end
function Base.getproperty(x::Ptr{var"##Ctag#335"}, f::Symbol)
    f === :data && return Ptr{Ptr{Cvoid}}(x + 0)
    f === :size && return Ptr{Cint}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#335", f::Symbol)
    r = Ref{var"##Ctag#335"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#335"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#335"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#336"
    rank::Cint
    shape::Ptr{Cint}
end
function Base.getproperty(x::Ptr{var"##Ctag#336"}, f::Symbol)
    f === :rank && return Ptr{Cint}(x + 0)
    f === :shape && return Ptr{Ptr{Cint}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#336", f::Symbol)
    r = Ref{var"##Ctag#336"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#336"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#336"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#338"
    size::Cint
    data::Ptr{Cvoid}
end
function Base.getproperty(x::Ptr{var"##Ctag#338"}, f::Symbol)
    f === :size && return Ptr{Cint}(x + 0)
    f === :data && return Ptr{Ptr{Cvoid}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#338", f::Symbol)
    r = Ref{var"##Ctag#338"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#338"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#338"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#342"
    data_ptr_size::Cint
end
function Base.getproperty(x::Ptr{var"##Ctag#342"}, f::Symbol)
    f === :tensor && return Ptr{Ptr{Cint}}(x + 0)
    f === :descriptor && return Ptr{Ptr{Cint}}(x + 0)
    f === :data_ptr && return Ptr{Ptr{Cvoid}}(x + 0)
    f === :data_ptr_size && return Ptr{Cint}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#342", f::Symbol)
    r = Ref{var"##Ctag#342"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#342"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#342"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


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

struct var"##Ctag#344"
    LSW::Int32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#344"}, f::Symbol)
    f === :LSW && return Ptr{Int32}(x + 0)
    f === :d3 && return Ptr{UInt32}(x + 4)
    f === :d2 && return Ptr{UInt32}(x + 8)
    f === :MSW && return Ptr{UInt32}(x + 12)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#344", f::Symbol)
    r = Ref{var"##Ctag#344"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#344"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#344"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#346"
    LSW::Int32
    d31::UInt32
    d30::UInt32
    d29::UInt32
    d28::UInt32
    d27::UInt32
    d26::UInt32
    d25::UInt32
    d24::UInt32
    d23::UInt32
    d22::UInt32
    d21::UInt32
    d20::UInt32
    d19::UInt32
    d18::UInt32
    d17::UInt32
    d16::UInt32
    d15::UInt32
    d14::UInt32
    d13::UInt32
    d12::UInt32
    d11::UInt32
    d10::UInt32
    d9::UInt32
    d8::UInt32
    d7::UInt32
    d6::UInt32
    d5::UInt32
    d4::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#346"}, f::Symbol)
    f === :LSW && return Ptr{Int32}(x + 0)
    f === :d31 && return Ptr{UInt32}(x + 4)
    f === :d30 && return Ptr{UInt32}(x + 8)
    f === :d29 && return Ptr{UInt32}(x + 12)
    f === :d28 && return Ptr{UInt32}(x + 16)
    f === :d27 && return Ptr{UInt32}(x + 20)
    f === :d26 && return Ptr{UInt32}(x + 24)
    f === :d25 && return Ptr{UInt32}(x + 28)
    f === :d24 && return Ptr{UInt32}(x + 32)
    f === :d23 && return Ptr{UInt32}(x + 36)
    f === :d22 && return Ptr{UInt32}(x + 40)
    f === :d21 && return Ptr{UInt32}(x + 44)
    f === :d20 && return Ptr{UInt32}(x + 48)
    f === :d19 && return Ptr{UInt32}(x + 52)
    f === :d18 && return Ptr{UInt32}(x + 56)
    f === :d17 && return Ptr{UInt32}(x + 60)
    f === :d16 && return Ptr{UInt32}(x + 64)
    f === :d15 && return Ptr{UInt32}(x + 68)
    f === :d14 && return Ptr{UInt32}(x + 72)
    f === :d13 && return Ptr{UInt32}(x + 76)
    f === :d12 && return Ptr{UInt32}(x + 80)
    f === :d11 && return Ptr{UInt32}(x + 84)
    f === :d10 && return Ptr{UInt32}(x + 88)
    f === :d9 && return Ptr{UInt32}(x + 92)
    f === :d8 && return Ptr{UInt32}(x + 96)
    f === :d7 && return Ptr{UInt32}(x + 100)
    f === :d6 && return Ptr{UInt32}(x + 104)
    f === :d5 && return Ptr{UInt32}(x + 108)
    f === :d4 && return Ptr{UInt32}(x + 112)
    f === :d3 && return Ptr{UInt32}(x + 116)
    f === :d2 && return Ptr{UInt32}(x + 120)
    f === :MSW && return Ptr{UInt32}(x + 124)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#346", f::Symbol)
    r = Ref{var"##Ctag#346"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#346"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#346"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#348"
    LSW::Int32
    d15::UInt32
    d14::UInt32
    d13::UInt32
    d12::UInt32
    d11::UInt32
    d10::UInt32
    d9::UInt32
    d8::UInt32
    d7::UInt32
    d6::UInt32
    d5::UInt32
    d4::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#348"}, f::Symbol)
    f === :LSW && return Ptr{Int32}(x + 0)
    f === :d15 && return Ptr{UInt32}(x + 4)
    f === :d14 && return Ptr{UInt32}(x + 8)
    f === :d13 && return Ptr{UInt32}(x + 12)
    f === :d12 && return Ptr{UInt32}(x + 16)
    f === :d11 && return Ptr{UInt32}(x + 20)
    f === :d10 && return Ptr{UInt32}(x + 24)
    f === :d9 && return Ptr{UInt32}(x + 28)
    f === :d8 && return Ptr{UInt32}(x + 32)
    f === :d7 && return Ptr{UInt32}(x + 36)
    f === :d6 && return Ptr{UInt32}(x + 40)
    f === :d5 && return Ptr{UInt32}(x + 44)
    f === :d4 && return Ptr{UInt32}(x + 48)
    f === :d3 && return Ptr{UInt32}(x + 52)
    f === :d2 && return Ptr{UInt32}(x + 56)
    f === :MSW && return Ptr{UInt32}(x + 60)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#348", f::Symbol)
    r = Ref{var"##Ctag#348"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#348"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#348"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#351"
    LSW::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#351"}, f::Symbol)
    f === :LSW && return Ptr{UInt32}(x + 0)
    f === :d3 && return Ptr{UInt32}(x + 4)
    f === :d2 && return Ptr{UInt32}(x + 8)
    f === :MSW && return Ptr{UInt32}(x + 12)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#351", f::Symbol)
    r = Ref{var"##Ctag#351"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#351"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#351"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#353"
    LSW::UInt32
    d7::UInt32
    d6::UInt32
    d5::UInt32
    d4::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#353"}, f::Symbol)
    f === :LSW && return Ptr{UInt32}(x + 0)
    f === :d7 && return Ptr{UInt32}(x + 4)
    f === :d6 && return Ptr{UInt32}(x + 8)
    f === :d5 && return Ptr{UInt32}(x + 12)
    f === :d4 && return Ptr{UInt32}(x + 16)
    f === :d3 && return Ptr{UInt32}(x + 20)
    f === :d2 && return Ptr{UInt32}(x + 24)
    f === :MSW && return Ptr{UInt32}(x + 28)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#353", f::Symbol)
    r = Ref{var"##Ctag#353"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#353"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#353"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#355"
    LSW::UInt32
    d31::UInt32
    d30::UInt32
    d29::UInt32
    d28::UInt32
    d27::UInt32
    d26::UInt32
    d25::UInt32
    d24::UInt32
    d23::UInt32
    d22::UInt32
    d21::UInt32
    d20::UInt32
    d19::UInt32
    d18::UInt32
    d17::UInt32
    d16::UInt32
    d15::UInt32
    d14::UInt32
    d13::UInt32
    d12::UInt32
    d11::UInt32
    d10::UInt32
    d9::UInt32
    d8::UInt32
    d7::UInt32
    d6::UInt32
    d5::UInt32
    d4::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#355"}, f::Symbol)
    f === :LSW && return Ptr{UInt32}(x + 0)
    f === :d31 && return Ptr{UInt32}(x + 4)
    f === :d30 && return Ptr{UInt32}(x + 8)
    f === :d29 && return Ptr{UInt32}(x + 12)
    f === :d28 && return Ptr{UInt32}(x + 16)
    f === :d27 && return Ptr{UInt32}(x + 20)
    f === :d26 && return Ptr{UInt32}(x + 24)
    f === :d25 && return Ptr{UInt32}(x + 28)
    f === :d24 && return Ptr{UInt32}(x + 32)
    f === :d23 && return Ptr{UInt32}(x + 36)
    f === :d22 && return Ptr{UInt32}(x + 40)
    f === :d21 && return Ptr{UInt32}(x + 44)
    f === :d20 && return Ptr{UInt32}(x + 48)
    f === :d19 && return Ptr{UInt32}(x + 52)
    f === :d18 && return Ptr{UInt32}(x + 56)
    f === :d17 && return Ptr{UInt32}(x + 60)
    f === :d16 && return Ptr{UInt32}(x + 64)
    f === :d15 && return Ptr{UInt32}(x + 68)
    f === :d14 && return Ptr{UInt32}(x + 72)
    f === :d13 && return Ptr{UInt32}(x + 76)
    f === :d12 && return Ptr{UInt32}(x + 80)
    f === :d11 && return Ptr{UInt32}(x + 84)
    f === :d10 && return Ptr{UInt32}(x + 88)
    f === :d9 && return Ptr{UInt32}(x + 92)
    f === :d8 && return Ptr{UInt32}(x + 96)
    f === :d7 && return Ptr{UInt32}(x + 100)
    f === :d6 && return Ptr{UInt32}(x + 104)
    f === :d5 && return Ptr{UInt32}(x + 108)
    f === :d4 && return Ptr{UInt32}(x + 112)
    f === :d3 && return Ptr{UInt32}(x + 116)
    f === :d2 && return Ptr{UInt32}(x + 120)
    f === :MSW && return Ptr{UInt32}(x + 124)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#355", f::Symbol)
    r = Ref{var"##Ctag#355"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#355"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#355"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#357"
    LSW::Int32
    d7::UInt32
    d6::UInt32
    d5::UInt32
    d4::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#357"}, f::Symbol)
    f === :LSW && return Ptr{Int32}(x + 0)
    f === :d7 && return Ptr{UInt32}(x + 4)
    f === :d6 && return Ptr{UInt32}(x + 8)
    f === :d5 && return Ptr{UInt32}(x + 12)
    f === :d4 && return Ptr{UInt32}(x + 16)
    f === :d3 && return Ptr{UInt32}(x + 20)
    f === :d2 && return Ptr{UInt32}(x + 24)
    f === :MSW && return Ptr{UInt32}(x + 28)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#357", f::Symbol)
    r = Ref{var"##Ctag#357"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#357"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#357"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


struct var"##Ctag#359"
    LSW::UInt32
    d15::UInt32
    d14::UInt32
    d13::UInt32
    d12::UInt32
    d11::UInt32
    d10::UInt32
    d9::UInt32
    d8::UInt32
    d7::UInt32
    d6::UInt32
    d5::UInt32
    d4::UInt32
    d3::UInt32
    d2::UInt32
    MSW::UInt32
end
function Base.getproperty(x::Ptr{var"##Ctag#359"}, f::Symbol)
    f === :LSW && return Ptr{UInt32}(x + 0)
    f === :d15 && return Ptr{UInt32}(x + 4)
    f === :d14 && return Ptr{UInt32}(x + 8)
    f === :d13 && return Ptr{UInt32}(x + 12)
    f === :d12 && return Ptr{UInt32}(x + 16)
    f === :d11 && return Ptr{UInt32}(x + 20)
    f === :d10 && return Ptr{UInt32}(x + 24)
    f === :d9 && return Ptr{UInt32}(x + 28)
    f === :d8 && return Ptr{UInt32}(x + 32)
    f === :d7 && return Ptr{UInt32}(x + 36)
    f === :d6 && return Ptr{UInt32}(x + 40)
    f === :d5 && return Ptr{UInt32}(x + 44)
    f === :d4 && return Ptr{UInt32}(x + 48)
    f === :d3 && return Ptr{UInt32}(x + 52)
    f === :d2 && return Ptr{UInt32}(x + 56)
    f === :MSW && return Ptr{UInt32}(x + 60)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#359", f::Symbol)
    r = Ref{var"##Ctag#359"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#359"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#359"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end


const vDSP_Version0 = 1126

const vDSP_Version1 = 100

const USE_NON_APPLE_STANDARD_DATATYPES = 1

# Skipping MacroDefinition: __VBASICOPS_INLINE_ATTR__ __attribute__ ( ( __always_inline__ , __nodebug__ ) )

# Skipping MacroDefinition: SPARSE_PUBLIC_INTERFACE __attribute__ ( ( overloadable ) )

const CBLAS_INDEX = Cint

# Skipping MacroDefinition: __SPARSE_ENUM_ATTR __attribute__ ( ( enum_extensibility ( open ) ) )

# Skipping MacroDefinition: __SPARSE_ENUM_ATTR_CLOSED __attribute__ ( ( enum_extensibility ( closed ) ) )

const _SPARSE_IMPLEMENTATION_TYPE = Float64

# Skipping MacroDefinition: _SPARSE_OLDSTYLE ( NAME ) NAME ## _double

# Skipping MacroDefinition: _SPARSE_VARIANT ( NAME ) NAME ## _Double

const _SPARSE_ATTRIBUTES = SparseAttributes_t

# Skipping MacroDefinition: BNNS_ENUM ( _name , _type , ... ) enum { __VA_ARGS__ } ; typedef _type _name

const BNNS_MAX_TENSOR_DIMENSION = 8

const QUADRATURE_INTEGRATE_QAG_WORKSPACE_PER_INTERVAL = 32

const QUADRATURE_INTEGRATE_QAGS_WORKSPACE_PER_INTERVAL = 152

end # module
