using Documenter
using AppleAccelerate

makedocs(;
    sitename = "AppleAccelerate.jl",
    modules = [AppleAccelerate],
    pages = [
        "Introduction" => "index.md",
        "Array Operations" => "array.md",
        "Complex Array Operations" => "complex.md",
        "Dense Linear Algebra" => "blas.md",
        "Sparse Linear Algebra" => "sparse.md",
        "FFT & Transforms" => "fft.md",
        "Filtering & Spectral" => "filtering.md",
        "Neural Network Primitives (BNNS)" => "bnns.md",
        "Image Processing (vImage)" => "vimage.md",
        "Numerical Integration" => "quadrature.md",
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
