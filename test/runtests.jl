module TestKmers

using Kmers
import BioSequences: @dna_str, @aa_str, @rna_str, LongSequence
using Test

# Return a random DNA/RNA sequence of the given length.
function random_seq(n::Integer, nts, probs)
    cumprobs = cumsum(probs)
    x = Vector{Char}(undef, n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    return String(x)
end

function random_seq(::Type{A}, n::Integer) where {A<:Alphabet}
    # TODO: Resolve the use of symbols(A()).
    nts = symbols(A())
    probs = Vector{Float64}(undef, length(nts))
    fill!(probs, 1 / length(nts))
    return LongSequence{A}(random_seq(n, nts, probs))
end

function random_dna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'T', 'N'], probs)
end

function random_rna(n, probs=[0.24, 0.24, 0.24, 0.24, 0.04])
    return random_seq(n, ['A', 'C', 'G', 'U', 'N'], probs)
end

function random_dna_kmer(len)
    return random_dna(len, [0.25, 0.25, 0.25, 0.25])
end

function random_rna_kmer(len)
    return random_rna(len, [0.25, 0.25, 0.25, 0.25])
end

@testset "Kmers" begin
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

end # module