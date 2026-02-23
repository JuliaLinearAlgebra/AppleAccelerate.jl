# API Reference

## Module

```@docs
AppleAccelerate
```

## BLAS/LAPACK

```@docs
AppleAccelerate.load_accelerate
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
AppleAccelerate.get_macos_version
```

## Array Operations

```@docs
AppleAccelerate.@replaceBase
```

## DSP & FFT

```@docs
AppleAccelerate.plan_fft
AppleAccelerate.fft
AppleAccelerate.destroy_fftsetup
AppleAccelerate.FFTSetup
AppleAccelerate.plan_dct
AppleAccelerate.dct
AppleAccelerate.conv
AppleAccelerate.conv!
AppleAccelerate.xcorr
AppleAccelerate.xcorr!
AppleAccelerate.biquadcreate
AppleAccelerate.biquad
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
AppleAccelerate.hann
```

## Sparse Linear Algebra

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.AAFactorization
AppleAccelerate.factor!
AppleAccelerate.solve
AppleAccelerate.solve!
AppleAccelerate.muladd!
```
