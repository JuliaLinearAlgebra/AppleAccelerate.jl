# API Reference

## Array Operations

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
AppleAccelerate.dot
AppleAccelerate.distancesq
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
AppleAccelerate.venvlp
AppleAccelerate.vaam
AppleAccelerate.vsbsbm
AppleAccelerate.vasbm
AppleAccelerate.vpythg
AppleAccelerate.vasm
AppleAccelerate.vsbsm
AppleAccelerate.vsma
AppleAccelerate.vsmsa
AppleAccelerate.vaddsub
```

### Clipping & Thresholding

```@docs
AppleAccelerate.vclip
AppleAccelerate.viclip
AppleAccelerate.vthr
AppleAccelerate.vthres
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

### Complex Array Operations

```@docs
AppleAccelerate.vconj
AppleAccelerate.vcopy
AppleAccelerate.vphase
AppleAccelerate.vmags
AppleAccelerate.vmagsa
AppleAccelerate.polar
AppleAccelerate.rect
```

## Dense Linear Algebra

```@docs
AppleAccelerate.load_accelerate
AppleAccelerate.set_num_threads
AppleAccelerate.get_num_threads
AppleAccelerate.get_macos_version
```

## Sparse Linear Algebra

```@docs
AppleAccelerate.AASparseMatrix
AppleAccelerate.AAFactorization
AppleAccelerate.muladd!
```

## Signal Processing

### DCT

```@docs
AppleAccelerate.plan_dct
AppleAccelerate.dct
AppleAccelerate.plan_destroy
```

### Spectral Analysis

```@docs
AppleAccelerate.zaspec
AppleAccelerate.zcspec
AppleAccelerate.zcoher
AppleAccelerate.ztrans
```

### Window Functions

```@docs
AppleAccelerate.blackman
AppleAccelerate.hamming
AppleAccelerate.hanning
```
