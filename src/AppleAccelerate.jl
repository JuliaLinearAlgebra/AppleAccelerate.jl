__precompile__()
module AppleAccelerate

using Libdl

try
    global const libacc = dlopen("/System/Library/Frameworks/Accelerate.framework/Accelerate")
catch
    error("Accelerate framework not found.")
end

get_fptr(s) = dlsym(libacc, s)

include("Array.jl")
include("DSP.jl")
include("Util.jl")

end # module
