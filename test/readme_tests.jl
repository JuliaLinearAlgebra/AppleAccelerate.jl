# Every ```julia fence in README.md must run as written, each in a fresh module — the
# README promises "self-contained, copy-pasteable" blocks, and nothing else in CI reads it.
# (The manual's examples are checked separately, by Documenter, in the docs build.)

using SparseArrays  # so the README's own `using SparseArrays` resolves in the test env

function readme_julia_blocks(path)
    blocks = Pair{Int,String}[]
    start, buf = 0, String[]
    for (n, line) in enumerate(eachline(path))
        if start == 0
            startswith(line, "```julia") && (start = n; empty!(buf))
        elseif startswith(line, "```")
            push!(blocks, start => join(buf, "\n"))
            start = 0
        else
            push!(buf, line)
        end
    end
    start == 0 || error("unterminated ```julia fence at README.md:$start")
    return blocks
end

@testset "README examples" begin
    readme = joinpath(dirname(@__DIR__), "README.md")
    blocks = readme_julia_blocks(readme)
    # Installation snippets touch the package manager, not the API.
    runnable = filter(b -> !occursin("Pkg.add", last(b)), blocks)
    @test length(runnable) >= 8   # one per subsystem; a fence-syntax change must not skip them all
    for (line, code) in runnable
        @testset "README.md:$line" begin
            m = Module(:ReadmeBlock)
            @test (Base.include_string(m, code, "README.md:$line"); true)
        end
    end
end
