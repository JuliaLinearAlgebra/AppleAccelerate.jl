# Compare benchmark/benchmarks.jl between two git revisions and print a `judge` report.
#
# Usage (from the repository root):
#   julia --project=benchmark benchmark/compare.jl <baseline> [<target>]
#
#   julia --project=benchmark benchmark/compare.jl v0.8.0          # v0.8.0 vs HEAD
#   julia --project=benchmark benchmark/compare.jl v0.8.0 my-branch
#
# Each revision is exported with `git archive` into a temp directory (the working tree is
# not touched, and uncommitted changes are NOT measured — commit first) and benchmarked
# in its own Julia process against its own temp environment. The suite definition is
# always *this* checkout's benchmarks.jl, so a baseline that predates benchmark/ still
# works; benchmarks whose API the baseline lacks are skipped and reported as unmatched.
#
# Exit status is 1 if any benchmark regressed, so the script can gate a release checklist.
# Run it on an otherwise idle machine: Accelerate's SME/AMX co-processor is shared, and a
# busy machine produces false regressions. `TOLERANCE` (default 0.10) is the relative
# time change treated as noise.

using BenchmarkTools, Printf
using BenchmarkTools: leaves

const REPO = dirname(@__DIR__)
const SUITE_FILE = joinpath(@__DIR__, "benchmarks.jl")
const TOLERANCE = parse(Float64, get(ENV, "TOLERANCE", "0.10"))

function benchmark_revision(rev::AbstractString)
    sha = readchomp(`git -C $REPO rev-parse --verify $(rev * "^{commit}")`)
    src = mktempdir(); env = mktempdir()
    run(pipeline(`git -C $REPO archive $sha`, `tar -x -C $src`))
    out = joinpath(env, "results.json")
    script = """
        using Pkg
        Pkg.develop(PackageSpec(path = $(repr(src))); io = devnull)
        Pkg.add("BenchmarkTools"; io = devnull)
        using BenchmarkTools
        include($(repr(SUITE_FILE)))
        tune!(SUITE)
        BenchmarkTools.save($(repr(out)), run(SUITE; verbose = true))
        """
    @info "Benchmarking $rev ($(first(sha, 8)))"
    run(`$(Base.julia_cmd()) --project=$env --threads=1 -e $script`)
    return only(BenchmarkTools.load(out))
end

function main(args)
    if !(1 <= length(args) <= 2)
        println(stderr, "usage: compare.jl <baseline> [<target>]")
        return 2
    end
    baseline_rev, target_rev = args[1], get(args, 2, "HEAD")
    baseline = benchmark_revision(baseline_rev)
    target = benchmark_revision(target_rev)
    judgement = judge(minimum(target), minimum(baseline); time_tolerance = TOLERANCE)

    println("\n", "="^78)
    println("$target_rev vs $baseline_rev   (minimum times, tolerance ±$(round(Int, 100TOLERANCE))%)")
    println("="^78)
    rows = sort!(collect(leaves(judgement)); by = r -> join(string.(first(r)), " "))
    nregress = 0
    for (key, j) in rows
        mark = j.time === :regression ? "REGRESSION" : j.time === :improvement ? "improved" : ""
        nregress += j.time === :regression
        @printf("%-52s %7.2fx  %s\n", join(string.(key), " / "), j.ratio.time, mark)
    end
    nmax = max(length(leaves(target)), length(leaves(baseline)))
    length(rows) < nmax &&
        println("\n$(nmax - length(rows)) benchmark(s) exist in only one revision and were not compared.")
    nimprove = count(r -> last(r).time === :improvement, rows)
    println("\n$nregress regression(s), $nimprove improvement(s), $(length(rows)) compared.")
    return nregress == 0 ? 0 : 1
end

exit(main(ARGS))
