# API Reference

## Array Operations — [vDSP](https://developer.apple.com/documentation/accelerate/vdsp) / [vecLib](https://developer.apple.com/documentation/accelerate/veclib)

### Unary vDSP Operations

```@docs
AppleAccelerate.vneg
AppleAccelerate.vnabs
AppleAccelerate.vabs
AppleAccelerate.vsq
AppleAccelerate.vssq
AppleAccelerate.vfrac
AppleAccelerate.vreverse!
AppleAccelerate.vreverse
```

### Vector Reductions

```@docs
AppleAccelerate.maximum
AppleAccelerate.minimum
AppleAccelerate.sum
AppleAccelerate.mean
AppleAccelerate.findmax
AppleAccelerate.findmin
AppleAccelerate.meanmag
AppleAccelerate.meansqr
AppleAccelerate.meanssqr
AppleAccelerate.summag
AppleAccelerate.sumsqr
AppleAccelerate.sumssqr
AppleAccelerate.dot
AppleAccelerate.distancesq
```

### Special Return Types

```@docs
AppleAccelerate.sincos
AppleAccelerate.cis
```

### Vector-Vector Arithmetic

```@docs
AppleAccelerate.vadd
AppleAccelerate.vadd!
AppleAccelerate.vsub
AppleAccelerate.vsub!
AppleAccelerate.vmul
AppleAccelerate.vmul!
AppleAccelerate.vdiv
AppleAccelerate.vdiv!
```

### Two-Vector Comparison & Distance

```@docs
AppleAccelerate.vmax
AppleAccelerate.vmin
AppleAccelerate.vmaxmg
AppleAccelerate.vminmg
AppleAccelerate.vdist
AppleAccelerate.vtmerg
```

### Vector-Scalar Operations

```@docs
AppleAccelerate.vsadd
AppleAccelerate.vsadd!
AppleAccelerate.vssub
AppleAccelerate.vssub!
AppleAccelerate.svsub
AppleAccelerate.svsub!
AppleAccelerate.vsmul
AppleAccelerate.vsmul!
AppleAccelerate.vsdiv
AppleAccelerate.vsdiv!
AppleAccelerate.svdiv
```

### Compound Arithmetic

```@docs
AppleAccelerate.vam
AppleAccelerate.vsbm
AppleAccelerate.vma
AppleAccelerate.vmsb
AppleAccelerate.venvlp
AppleAccelerate.vaam
AppleAccelerate.vsbsbm
AppleAccelerate.vasbm
AppleAccelerate.vmma
AppleAccelerate.vmmsb
AppleAccelerate.vpythg
AppleAccelerate.vasm
AppleAccelerate.vsbsm
AppleAccelerate.vsma
AppleAccelerate.vsmsa
AppleAccelerate.vmsa
AppleAccelerate.vsmsb
AppleAccelerate.vsmsma
AppleAccelerate.vaddsub
```

### Extra Reductions

```@docs
AppleAccelerate.rmsqv
AppleAccelerate.sve_svesq
AppleAccelerate.maxmgv
AppleAccelerate.minmgv
AppleAccelerate.maxmgvi
AppleAccelerate.minmgvi
```

### Clipping & Thresholding

```@docs
AppleAccelerate.vclip
AppleAccelerate.vclipc
AppleAccelerate.viclip
AppleAccelerate.vthr
AppleAccelerate.vthres
AppleAccelerate.vlim
AppleAccelerate.vthrsc
AppleAccelerate.vcmprs
```

### Type Conversion

```@docs
AppleAccelerate.vdouble
AppleAccelerate.vsingle
```

### Ramp Generation

```@docs
AppleAccelerate.vramp
AppleAccelerate.vrampmul
AppleAccelerate.vrampmul2
```

### Linear Average

```@docs
AppleAccelerate.vavlin
```

### Integration & Running Operations

```@docs
AppleAccelerate.vrsum
AppleAccelerate.vsimps
AppleAccelerate.vtrapz
AppleAccelerate.vswsum
AppleAccelerate.vswmax
```

### Interpolation

```@docs
AppleAccelerate.vintb
AppleAccelerate.vlint
AppleAccelerate.vqint
```

### Polynomial Evaluation

```@docs
AppleAccelerate.vpoly
```

### Normalization

```@docs
AppleAccelerate.vnormalize
```

### Zero Crossings

```@docs
AppleAccelerate.nzcros
```

### Decibel Conversion

```@docs
AppleAccelerate.vdbcon
```

### Vector Fill, Swap & Sort

```@docs
AppleAccelerate.vclr!
AppleAccelerate.vfill!
AppleAccelerate.vswap!
AppleAccelerate.vsort!
AppleAccelerate.vsorti
```

### Gathering & Indexing

```@docs
AppleAccelerate.vgathr
AppleAccelerate.vindex
AppleAccelerate.vgen
AppleAccelerate.vgenp
AppleAccelerate.vtabi
```

### Matrix Operations

```@docs
AppleAccelerate.mmul
AppleAccelerate.mtrans
AppleAccelerate.mmov
```

### Integer Operations

```@docs
AppleAccelerate.vaddi
AppleAccelerate.vabsi
AppleAccelerate.vfilli!
AppleAccelerate.veqvi
```

### Image Convolution

```@docs
AppleAccelerate.f3x3
AppleAccelerate.f5x5
AppleAccelerate.imgfir
```

### Format Conversion

```@docs
AppleAccelerate.ctoz
AppleAccelerate.ztoc
```

### Complex Array Operations

```@docs
AppleAccelerate.vconj
AppleAccelerate.vcopy
AppleAccelerate.vphase
AppleAccelerate.vmags
AppleAccelerate.vmagsa
AppleAccelerate.polar
AppleAccelerate.rect
AppleAccelerate.zvadd
AppleAccelerate.zvsub
AppleAccelerate.zvcmul
AppleAccelerate.zrvmul
AppleAccelerate.zrvdiv
AppleAccelerate.zrvadd
AppleAccelerate.zrvsub
AppleAccelerate.zvcma
AppleAccelerate.zvma
AppleAccelerate.zvsma
AppleAccelerate.zidotpr
AppleAccelerate.zrdotpr
AppleAccelerate.zvfill!
AppleAccelerate.zconv
AppleAccelerate.zmmul
```

## Dense Linear Algebra — [BLAS](https://developer.apple.com/documentation/accelerate/blas) / [LAPACK](https://developer.apple.com/documentation/accelerate/solving-systems-of-linear-equations-with-lapack)

```@docs
AppleAccelerate.load_accelerate
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
AppleAccelerate.get_macos_version
AppleAccelerate._read_macos_version
```

## Sparse Linear Algebra — [Sparse Solvers](https://developer.apple.com/documentation/accelerate/sparse_solvers)

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.AAFactorization
AppleAccelerate.muladd!
AppleAccelerate.factor!
AppleAccelerate.solve
AppleAccelerate.solve!
```

## Signal Processing — [vDSP](https://developer.apple.com/documentation/accelerate/vdsp)

### FFT — [Fast Fourier Transforms](https://developer.apple.com/documentation/accelerate/fast_fourier_transforms)

```@docs
AppleAccelerate.plan_fft
AppleAccelerate.fft
AppleAccelerate.ifft
AppleAccelerate.bfft
AppleAccelerate.fft!
AppleAccelerate.ifft!
AppleAccelerate.bfft!
```

### Real FFT

```@docs
AppleAccelerate.plan_rfft
AppleAccelerate.rfft
AppleAccelerate.irfft
AppleAccelerate.brfft
```

### DFT — [`vDSP_DFT_Execute`](https://developer.apple.com/documentation/accelerate/vdsp_dft_execute)

```@docs
AppleAccelerate.plan_dft
AppleAccelerate.dft
AppleAccelerate.idft
```

### DCT — [Discrete Cosine Transforms](https://developer.apple.com/documentation/accelerate/discrete_cosine_transforms)

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
AppleAccelerate.plan_destroy
```

### Convolution & Correlation — [`vDSP_conv`](https://developer.apple.com/documentation/accelerate/vdsp_conv)

```@docs
AppleAccelerate.conv
AppleAccelerate.conv!
AppleAccelerate.xcorr
AppleAccelerate.xcorr!
```

### Biquad Filtering — [`vDSP_biquad`](https://developer.apple.com/documentation/accelerate/vdsp_biquad)

```@docs
AppleAccelerate.biquadcreate
AppleAccelerate.biquad
AppleAccelerate.biquaddestroy
```

### Multi-Channel Biquad Filtering — [`vDSP_biquadm`](https://developer.apple.com/documentation/accelerate/vdsp_biquadm)

```@docs
AppleAccelerate.biquadm_create
AppleAccelerate.biquadm
```

### Recursive Filter

```@docs
AppleAccelerate.deq22
AppleAccelerate.deq22!
AppleAccelerate.desamp
AppleAccelerate.desamp!
AppleAccelerate.wiener
AppleAccelerate.wiener!
```

### Spectral Analysis

```@docs
AppleAccelerate.zaspec
AppleAccelerate.zaspec!
AppleAccelerate.zcspec
AppleAccelerate.zcspec!
AppleAccelerate.zcoher
AppleAccelerate.zcoher!
AppleAccelerate.ztrans
AppleAccelerate.ztrans!
```

### Window Functions — [Vector Generation](https://developer.apple.com/documentation/accelerate/vdsp/vector_generation)

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
```
