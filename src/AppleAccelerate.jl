module AppleAccelerate

using Libdl

if Sys.isapple()

    try
        global const libacc = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
        let accel_lib = dlopen(libacc)
            dlclose(accel_lib)
        end
    catch
        error("Accelerate framework not found.")
    end


    include("Array.jl")
    include("DSP.jl")
    include("Util.jl")

end

end # module
