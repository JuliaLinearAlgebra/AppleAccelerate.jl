# ============================================================
# Idiomatic LinearAlgebra integration backed by Accelerate vDSP.
#
# AppleAccelerate already forwards Julia's BLAS/LAPACK to Accelerate at load time
# (see `load_accelerate`), so `*`, `mul!`, factorizations, etc. on standard
# arrays are already accelerated through libblastrampoline. This file adds an
# explicit, opt-in vDSP path (vDSP_mmul / vDSP_mtrans) for plain
# `Matrix{Float32}` / `Matrix{Float64}` operands, exposed through the
# AppleAccelerate-owned verbs `accelerate_mul!` and `accelerate_transpose!`.
#
# These are deliberately *new* functions owned by AppleAccelerate rather than
# methods added to `LinearAlgebra.mul!` / `transpose!` for `Matrix{<:BlasFloat}`:
# the latter would be type piracy (all argument types belong to Base/stdlib),
# would silently shadow Julia's well-tested BLAS `gemm`, and is flagged by Aqua.
# Users who want the direct vDSP kernel call `accelerate_mul!` explicitly; the
# ordinary `*` / `mul!` continue to use Accelerate via BLAS forwarding.
#
# Note on packaging: `LinearAlgebra` is a *hard* dependency of AppleAccelerate
# (the main module and `sparse.jl` reference `LinearAlgebra.BlasInt`,
# `LinearAlgebra.Factorization`, `LinearAlgebra.ldiv!`, etc. at load time, and
# BLAS/LAPACK forwarding in `__init__` requires `LinearAlgebra.BLAS`). Because of
# that it cannot also serve as an extension *trigger* — Julia/Pkg forbids a
# package from being both a used `[deps]` entry and a weakdep extension trigger
# (verified empirically: such a config breaks `Pkg.test` precompilation). So this
# LinearAlgebra integration lives in the main module rather than in a separate
# `ext/` file. See the PR description for the full rationale.
# ============================================================

const _VDSPFloat = Union{Float32,Float64}

"""
    accelerate_mul!(C, A, B) -> C

Compute the matrix product `C = A * B` using Accelerate's vDSP `vDSP_mmul`
kernel and store the result in `C`. Defined for `Matrix{Float32}` and
`Matrix{Float64}` arguments.

This is an explicit, opt-in alternative to `LinearAlgebra.mul!`. Unlike BLAS
`gemm`, the vDSP path does not support transposed or scaled operands; use
`LinearAlgebra.mul!` / `*` (also Accelerate-backed, via BLAS forwarding) for the
general case.

# Examples
```julia
using AppleAccelerate, LinearAlgebra
A = rand(4, 3); B = rand(3, 5); C = Matrix{Float64}(undef, 4, 5)
AppleAccelerate.accelerate_mul!(C, A, B)  # C ≈ A * B
```
"""
accelerate_mul!(C::Matrix{T}, A::Matrix{T}, B::Matrix{T}) where {T<:_VDSPFloat} =
    mmul!(C, A, B)

"""
    accelerate_mul(A, B) -> C

Allocating form of [`accelerate_mul!`](@ref): return `C = A * B` computed with
Accelerate's vDSP `vDSP_mmul` kernel, for `Matrix{Float32}` / `Matrix{Float64}`.
"""
accelerate_mul(A::Matrix{T}, B::Matrix{T}) where {T<:_VDSPFloat} = mmul(A, B)

"""
    accelerate_transpose!(C, A) -> C

Compute the matrix transpose `C = transpose(A)` using Accelerate's vDSP
`vDSP_mtrans` kernel and store the result in `C`. Defined for `Matrix{Float32}`
and `Matrix{Float64}` arguments. `C` must have size `reverse(size(A))`.
"""
function accelerate_transpose!(C::Matrix{T}, A::Matrix{T}) where {T<:_VDSPFloat}
    size(C) == reverse(size(A)) ||
        throw(DimensionMismatch("transpose destination has size $(size(C)), expected $(reverse(size(A)))"))
    return mtrans!(C, A)
end

"""
    accelerate_transpose(A) -> C

Allocating form of [`accelerate_transpose!`](@ref): return `C = transpose(A)`
computed with Accelerate's vDSP `vDSP_mtrans` kernel.
"""
accelerate_transpose(A::Matrix{T}) where {T<:_VDSPFloat} = mtrans(A)
