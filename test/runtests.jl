module TestKmers

using Kmers
using BioSequences
using Test

include("utils.jl")

@testset "BioSequences Interface" begin
    for A in [DNAAlphabet{2}, DNAAlphabet{4}, RNAAlphabet{2}, RNAAlphabet{4}, AminoAcidAlphabet]
        for K in (1, 9, 116)
            @test BioSequences.has_interface(
                BioSequence,
                Kmers.derive_type(Kmer{A, K}),
                rand(collect(A()), K),
                false,
            )
        end
    end
end

@testset "Construction" begin
end

# include("construction_and_conversion.jl")
# include("comparisons.jl")
# include("length.jl")
# include("access.jl")
# include("random.jl")
# include("find.jl")
# include("print.jl")
# include("transformations.jl")
# include("mismatches.jl")
# include("debruijn_neighbors.jl")
# include("iteration.jl")
# include("translation.jl")
#include("shuffle.jl")

end # module
