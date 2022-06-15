module TestKmers

using Kmers
using BioSequences 
using Test

const GROUP = get(ENV, "GROUP", "All")

include("utils.jl")

if GROUP == "BioSequences" || GROUP == "All"
    include("biosequences_interface.jl")
    include("conversion.jl")
    include("comparisons.jl")
    include("length.jl")
    include("access.jl")
    include("random.jl")
    include("find.jl")
    include("print.jl")
    include("transformations.jl")
    include("mismatches.jl")
    include("debruijn_neighbors.jl")
    include("iteration.jl")
    include("translation.jl")
    #include("shuffle.jl")
end

end # module
