# API Reference

## BLAS/LAPACK

- [`AppleAccelerate.load_accelerate`](@ref)
- [`AppleAccelerate.get_macos_version`](@ref)
- [`AppleAccelerate.set_num_threads`](@ref)
- [`AppleAccelerate.get_num_threads`](@ref)

## Array Operations

- [`AppleAccelerate.@replaceBase`](@ref)
- Vector-vector: `vadd`, `vsub`, `vmul`, `vdiv` (and `!` variants)
- Vector-scalar: `vsadd`, `vssub`, `svsub`, `vsdiv`, `vsmul`, `svdiv` (and `!` variants)
- Reductions: `maximum`, `minimum`, `mean`, `sum`, `findmax`, `findmin`, `meanmag`, `meansqr`, `summag`, `sumsqr`, `dot`, `distancesq`
- Unary vDSP: `vneg`, `vnabs`, `vabs`, `vsq`, `vssq`, `vfrac`, `vreverse` (and `!` variants)
- Two-vector: `vmax`, `vmin`, `vmaxmg`, `vminmg`, `vdist`, `vtmerg` (and `!` variants)
- Compound arithmetic: `vam`, `vsbm`, `vaam`, `vsbsbm`, `vasbm`, `vpythg`, `vasm`, `vsbsm`, `vsma`, `vsmsa`, `vaddsub`, `venvlp` (and `!` variants)
- Clipping: `vclip`, `viclip`, `vthr`, `vthres`, `vcmprs` (and `!` variants)
- Type conversion: `vdouble`, `vsingle`
- Ramp generation: `vramp`, `vrampmul`
- Integration: `vrsum`, `vsimps`, `vtrapz`, `vswsum`, `vswmax` (and `!` variants)
- Interpolation: `vintb`, `vlint`, `vqint` (and `!` variants)
- Polynomial: `vpoly` (and `!` variant)
- Normalization: `vnormalize` (and `!` variant)
- Zero crossings: `nzcros`
- Decibel conversion: `vdbcon` (and `!` variant)

## DSP & FFT

- `plan_fft`, `fft`, `destroy_fftsetup` — FFT (Float32 and Float64)
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
