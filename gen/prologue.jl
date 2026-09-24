# Accelerate ships as a macOS *system framework*. Unlike the usual Clang.jl + JLL
# flow, there is no artifact to load — every Mac already has it. We ccall it by
# absolute path, mirroring the `libacc` const in the top-level AppleAccelerate module.
# NOTE: must stay in sync with `ACCELERATE_FRAMEWORK` in gen/generate.jl, which dlopens
# the same binary for its dead-symbol strip pass and verifies this line at generation time.
const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

# vecLibTypes.h defines vUInt32 as a 16-byte vector of unsigned int. Clang.jl omits
# its definition but still references it in the vBigNum union accessors. VecElement
# preserves the SIMD representation and 16-byte alignment; a plain UInt32 tuple does not.
const vUInt32 = NTuple{4, VecElement{UInt32}}

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

struct bnns_user_message_data_t
    data::Ptr{Cvoid}
    size::Csize_t
end

# bnns_graph_shape_t is NOT one of the `{ data, size }` handles above: it is
# `{ size_t rank; uint64_t *shape; }` — the count comes first. Clang.jl drops it for the same
# trailing-availability-attribute reason, so define it here with the header's field order.
struct bnns_graph_shape_t
    rank::Csize_t
    shape::Ptr{UInt64}
end

# bnns_graph_argument_t holds an anonymous union of three pointer types
# (BNNSTensor* / BNNSNDArrayDescriptor* / void*) followed by a size. Clang.jl drops structs
# with anonymous unions, so define it here. The union collapses to one pointer-sized field.
struct bnns_graph_argument_t
    data::Ptr{Cvoid}            # union { tensor; descriptor; data_ptr }
    data_ptr_size::Csize_t
end
