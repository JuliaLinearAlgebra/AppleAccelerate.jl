__precompile__()
module AppleAccelerate

const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"

include("Array.jl")
include("DSP.jl")
include("Util.jl")

end # module
