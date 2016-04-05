__precompile__()
module AppleAccelerate

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

if !isfile(libacc)
    error("Accelerate framework not found at $(libacc)")
end

include("Array.jl")
include("DSP.jl")
include("Util.jl")

end # module
