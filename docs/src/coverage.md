# Binding Coverage

How much of each Accelerate subframework has a hand-written idiomatic Julia API on top
of the generated `LibAccelerate` raw bindings. **Total** counts the C functions present
in the generated raw layer (which itself mirrors the in-scope headers); **Wrapped**
counts those referenced by name in the idiomatic layer.

!!! note
    Modules that construct C symbol names dynamically (e.g. `Symbol("vv", fa, suff)`)
    are undercounted by the static scan, so these numbers are a floor on real coverage.
    BLAS/LAPACK are intentionally excluded — they are forwarded via libblastrampoline,
    not ccall.

| Subframework | Wrapped | Total | Coverage |
|--------------|--------:|------:|---------:|
| vForce | 4 | 84 | 5% |
| vBigNum | 0 | 69 | 0% |
| vectorOps | 8 | 9 | 89% |
| vDSP | 164 | 461 | 36% |
| Sparse | 8 | 144 | 6% |
| BNNS | 0 | 139 | 0% |
| Quadrature | 1 | 1 | 100% |
| other | 0 | 1 | 0% |
| **Total** | **185** | **908** | **20%** |

