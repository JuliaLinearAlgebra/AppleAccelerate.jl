using Documenter
using AppleAccelerate

makedocs(;
    sitename = "AppleAccelerate.jl",
    modules = [AppleAccelerate],
    pages = [
        "Home" => "index.md",
        "BLAS & LAPACK" => "blas.md",
        "Array Operations" => "array.md",
        "DSP & FFT" => "dsp.md",
        "Sparse Linear Algebra" => "sparse.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
)

deploydocs(;
    repo = "github.com/JuliaLinearAlgebra/AppleAccelerate.jl.git",
    push_preview = true,
)
