# API Reference

## Array Operations

- [`AppleAccelerate.@replaceBase`](@ref)
- Vector-vector: `vadd`, `vsub`, `vmul`, `vdiv` (and `!` variants)
- Vector-scalar: `vsadd`, `vssub`, `svsub`, `vsdiv`, [`svdiv`](@ref AppleAccelerate.svdiv), `vsmul` (and `!` variants)
- Unary vDSP: [`vneg`](@ref AppleAccelerate.vneg), [`vnabs`](@ref AppleAccelerate.vnabs), [`vabs`](@ref AppleAccelerate.vabs), [`vsq`](@ref AppleAccelerate.vsq), [`vssq`](@ref AppleAccelerate.vssq), [`vfrac`](@ref AppleAccelerate.vfrac), [`vreverse!`](@ref AppleAccelerate.vreverse!) (and `!` variants)
- Two-vector: [`vmax`](@ref AppleAccelerate.vmax), [`vmin`](@ref AppleAccelerate.vmin), [`vmaxmg`](@ref AppleAccelerate.vmaxmg), [`vminmg`](@ref AppleAccelerate.vminmg), [`vdist`](@ref AppleAccelerate.vdist), [`vtmerg`](@ref AppleAccelerate.vtmerg) (and `!` variants)
- Compound arithmetic: [`vam`](@ref AppleAccelerate.vam), [`vsbm`](@ref AppleAccelerate.vsbm), [`venvlp`](@ref AppleAccelerate.venvlp), [`vaam`](@ref AppleAccelerate.vaam), [`vsbsbm`](@ref AppleAccelerate.vsbsbm), [`vasbm`](@ref AppleAccelerate.vasbm), [`vpythg`](@ref AppleAccelerate.vpythg), [`vasm`](@ref AppleAccelerate.vasm), [`vsbsm`](@ref AppleAccelerate.vsbsm), [`vsma`](@ref AppleAccelerate.vsma), [`vsmsa`](@ref AppleAccelerate.vsmsa), [`vaddsub`](@ref AppleAccelerate.vaddsub) (and `!` variants)
- Reductions: `maximum`, `minimum`, `mean`, `sum`, `findmax`, `findmin`, `meanmag`, `meansqr`, `summag`, `sumsqr`, [`dot`](@ref AppleAccelerate.dot), [`distancesq`](@ref AppleAccelerate.distancesq)
- Clipping: [`vclip`](@ref AppleAccelerate.vclip), [`viclip`](@ref AppleAccelerate.viclip), [`vthr`](@ref AppleAccelerate.vthr), [`vthres`](@ref AppleAccelerate.vthres), [`vcmprs`](@ref AppleAccelerate.vcmprs) (and `!` variants)
- Type conversion: [`vdouble`](@ref AppleAccelerate.vdouble), [`vsingle`](@ref AppleAccelerate.vsingle)
- Ramp generation: [`vramp`](@ref AppleAccelerate.vramp), [`vrampmul`](@ref AppleAccelerate.vrampmul)
- Integration: [`vrsum`](@ref AppleAccelerate.vrsum), [`vsimps`](@ref AppleAccelerate.vsimps), [`vtrapz`](@ref AppleAccelerate.vtrapz), [`vswsum`](@ref AppleAccelerate.vswsum), [`vswmax`](@ref AppleAccelerate.vswmax) (and `!` variants)
- Interpolation: [`vintb`](@ref AppleAccelerate.vintb), [`vlint`](@ref AppleAccelerate.vlint), [`vqint`](@ref AppleAccelerate.vqint) (and `!` variants)
- Polynomial: [`vpoly`](@ref AppleAccelerate.vpoly) (and `!` variant)
- Normalization: [`vnormalize`](@ref AppleAccelerate.vnormalize) (and `!` variant)
- Zero crossings: [`nzcros`](@ref AppleAccelerate.nzcros)
- Decibel conversion: [`vdbcon`](@ref AppleAccelerate.vdbcon) (and `!` variant)

## Dense Linear Algebra

- [`AppleAccelerate.load_accelerate`](@ref)
- [`AppleAccelerate.get_macos_version`](@ref)
- [`AppleAccelerate.set_num_threads`](@ref)
- [`AppleAccelerate.get_num_threads`](@ref)

## Sparse Linear Algebra

- [`AppleAccelerate.AASparseMatrix`](@ref) — Sparse matrix wrapper
- [`AppleAccelerate.AAFactorization`](@ref) — Lazy factorization
- [`AppleAccelerate.muladd!`](@ref) — Fused multiply-add (`y += A*x`)
- `solve`, `solve!` — Sparse solvers
- `factor!` — Explicit factorization (Cholesky, LDLT, QR, CholeskyAtA)
- `factorize` — Create `AAFactorization` from `AASparseMatrix`
- `ldiv!` — In-place left division via `LinearAlgebra`
- `transpose` — Transpose (flag-based, no copy)
- `issymmetric`, `istriu`, `istril` — Structure queries

## Signal Processing

- `fft`, `ifft`, `bfft`, `fft!`, `ifft!`, `bfft!` — Complex FFT (Float32 and Float64, 1D and 2D)
- `rfft`, `irfft`, `brfft` — Real FFT (Float32 and Float64, 1D)
- `plan_fft`, `plan_rfft` — FFT plan creation
- [`AppleAccelerate.plan_dct`](@ref), [`AppleAccelerate.dct`](@ref) — DCT (Float32 only)
- `conv`, `conv!` — Convolution
- `xcorr`, `xcorr!` — Cross-correlation
- `biquadcreate`, `biquad` — Biquad IIR filtering (Float64 only)
- [`AppleAccelerate.blackman`](@ref), [`AppleAccelerate.hamming`](@ref), [`AppleAccelerate.hanning`](@ref), [`AppleAccelerate.hann`](@ref) — Window functions
