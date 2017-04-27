module Sketches

export CountMinSketch,
       savesketch,
       loadsketch,
       collisionrate

include("base.jl")
include("countmin.jl")
include("io.jl")

end
