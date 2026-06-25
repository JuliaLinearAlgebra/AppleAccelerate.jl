using Documenter
using AppleAccelerate

makedocs(;
    sitename = "AppleAccelerate.jl",
    modules = [AppleAccelerate],
    pages = [
        "Introduction" => "index.md",
        "Array Operations" => "array.md",
        "Dense Linear Algebra" => "blas.md",
        "Sparse Linear Algebra" => "sparse.md",
        "FFT & Transforms" => "fft.md",
        "Filtering & Spectral" => "filtering.md",
        "Benchmarks" => "benchmarks.md",
        "Binding Coverage" => "coverage.md",
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
