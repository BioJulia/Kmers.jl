using Kmers
import BioSequences: @dna_str, @aa_str, @rna_str
using Test

@testset "Mers" begin
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
    #include("mers/shuffle.jl")
end
