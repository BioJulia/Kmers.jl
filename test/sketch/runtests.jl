module TestSketches

include("../utils.jl")

using Base.Test
using Kmers.Sketches
import Kmers.Sketches: add!

include("basic.jl")
include("countmin.jl")

end # module TestSketches
