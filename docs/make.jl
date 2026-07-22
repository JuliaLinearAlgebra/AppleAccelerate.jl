using Documenter
using AppleAccelerate

makedocs(;
    sitename = "AppleAccelerate.jl",
    modules = [AppleAccelerate],
    pages = [
        "Introduction" => "index.md",
        "Array Operations (vDSP / vForce)" => "array.md",
        "Complex Array Operations (vDSP)" => "complex.md",
        "Dense Linear Algebra (BLAS / LAPACK)" => "blas.md",
        "Sparse Linear Algebra (libSparse)" => "sparse.md",
        "FFT & Transforms (vDSP)" => "fft.md",
        "Filtering & Spectral (vDSP)" => "filtering.md",
        "Neural Network Primitives (BNNS)" => "bnns.md",
        "Image Processing (vImage)" => "vimage.md",
        "Numerical Integration (Quadrature)" => "quadrature.md",
        "Benchmarks" => "benchmarks.md",
        "Architecture" => "extensions.md",
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    warnonly = [:missing_docs],
)

deploydocs(;
    repo = "github.com/JuliaLinearAlgebra/AppleAccelerate.jl.git",
    push_preview = true,
)
