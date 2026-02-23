# API Reference

## BLAS/LAPACK

- [`AppleAccelerate.load_accelerate`](@ref)
- [`AppleAccelerate.get_macos_version`](@ref)
- [`AppleAccelerate.set_num_threads`](@ref)
- [`AppleAccelerate.get_num_threads`](@ref)

## Array Operations

- [`AppleAccelerate.@replaceBase`](@ref)
- Vector-vector: `vadd`, `vsub`, `vmul`, `vdiv` (and `!` variants)
- Vector-scalar: `vsadd`, `vssub`, `svsub`, `vsdiv`, `vsmul` (and `!` variants)
- Reductions: `maximum`, `minimum`, `mean`, `sum`, `findmax`, `findmin`, `meanmag`, `meansqr`, `summag`, `sumsqr`

## DSP & FFT

- `fft`, `ifft`, `bfft`, `plan_fft` — FFT (Float32 and Float64, 1D and 2D)
- [`AppleAccelerate.plan_dct`](@ref), [`AppleAccelerate.dct`](@ref) — DCT (Float32 only)
- `conv`, `conv!` — Convolution
- `xcorr`, `xcorr!` — Cross-correlation
- `biquadcreate`, `biquad` — Biquad IIR filtering (Float64 only)
- [`AppleAccelerate.blackman`](@ref), [`AppleAccelerate.hamming`](@ref), [`AppleAccelerate.hanning`](@ref), [`AppleAccelerate.hann`](@ref) — Window functions

## Sparse Linear Algebra

- [`AppleAccelerate.AASparseMatrix`](@ref) — Sparse matrix wrapper
- [`AppleAccelerate.AAFactorization`](@ref) — Lazy factorization
- [`AppleAccelerate.muladd!`](@ref) — Fused multiply-add
- `solve`, `solve!` — Sparse solvers
- `factor!` — Explicit factorization
