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
| vForce | 8 | 84 | 10% |
| vBigNum | 0 | 69 | 0% |
| vectorOps | 8 | 9 | 89% |
| vDSP | 164 | 460 | 36% |
| Sparse | 14 | 107 | 13% |
| BNNS | 5 | 136 | 4% |
| Quadrature | 1 | 1 | 100% |
| other | 0 | 1 | 0% |
| **Total** | **200** | **867** | **23%** |

